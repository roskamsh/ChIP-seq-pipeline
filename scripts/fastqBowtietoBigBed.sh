#!/bin/bash
#set -o errexit
#set -o pipefail


################################################################################
# Set default values for mapping parameters
################################################################################

nrMM="2"
mapLen="36"
gap="1000"
tempDir="/home/groups/CEDAR/roskamsh/test_data/temp"
USE=4
KEEPSEQ="0"
KEEPMM="0"
KEEPLOG="1"
FILTER="0"

################################################################################
# Requirements
################################################################################

# Programs:
# * fastx-toolkit
# * bowtie 0.12.9
# * Jim Kent's utils
# * bedtools 2.19.1

# Tools you will need for this analysis
module load fastx-toolkit/0.0.13
module load kent-ucsc/3.8f6f5e0a1cb75
module load bowtie/0.12.9
module load bedtools/2.19.1

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
echo >&2 "
$(basename $0) - Mapping FASTQ file to reference genome and returning BigBed file with coordinates of mapped reads (or fragments for paired end data)

USAGE: $(basename $0) -i <FASTQ file(s)> -n <outfile basename> -s <species> -l <read lengthn> -m <max mismatches> -o <output directory> [OPTIONS]
-i    FASTQ (gzipped) file(s) [required]
Use -i \"read_1.fastq.gz|read_2.fastq.gz\" to map paired-end gzipped fastq files
-n    Basename of the output file [required]
-a    Genome assembly (hg19, dm3, etc.) [required]
-l    Length to which sequenced reads will be trimmed before mapping (50 for human, 36 for fly) [default: $mapLen]
-m    Maximal number of mismatches for mapping (usually 3 for human, 2 for fly) [default: $nrMM]
-g    Maximal gap (in base pairs) between paired end reads [default: $gap]
-f     Filters [default: $FILTER]
        (0) Report all mapped reads
        (1) Report fragments collapsed on identical chromosome, start, end,
            and strand information
        (A) Report both all mapped reads AND unique fragments,
            i.e. report results from both filters 0 and 1
-r     Report mapping statistics (0/1) [default: $KEEPLOG]
-S     Report read/fragment sequences tags (0/1) [default: $KEEPSEQ]
-M     Report number of mismatches (column 5 of BED/BB file; 0/1) [default: $KEEPMM]
-u     Number of processes to use (will be overwritten by NSLOTS on SGE, see below) [default: $USE]
-o    Output directory [required]
-T    Directory to be used for writing temporary files [default: $tempDir]

NOTES:
Note that you need the bowtie index files for your genome assembly.
In addition the programs bowtie, fastx-toolkit and Jim Kent's utils must be installed on the machine.
"
exit 1
fi


################################################################################
# Parse input and error checking
################################################################################

while getopts "i:n:l:m:g:a:f:r:S:M:u:o:T" op
do
    case "$op" in
        i)  fastqFile="$OPTARG";;
        n)  sample="$OPTARG";;
        l)  mapLen="$OPTARG";;
        m)  nrMM="$OPTARG";;
        g)  gap="$OPTARG";;
        f)  FILTER="$OPTARG";;
        r)  KEEPLOG="$OPTARG";;
        u)  USE="$OPTARG";;
        a)  assembly="$OPTARG";;
        S)  KEEPSEQ="$OPTARG";;
        M)  KEEPMM="$OPTARG";;
        o)  outDir="$OPTARG";;
        T)  tempDir="$OPTARG";;
        \?) exit 1;;
    esac
done

if [ -z "$fastqFile" ]; then
    echo >&2 "ERROR: Specify input FASTQ file(s)!"
    exit 1
fi

if [ -z "$mapLen" ]; then
    echo >&2 "ERROR: Specify desired length of the read for mapping!"
    exit 1
fi

if [ -z "$sample" ]; then
    echo >&2 "ERROR: Specify basename of the output file!"
    exit 1
fi

if [ -z "$assembly" ]; then
    echo >&2 "ERROR: Specify genome assembly!"
    exit 1
fi

if [ -z "$outDir" ]; then
    echo >&2 "ERROR: Specify output directory!"
    exit 1
fi

if [ ! -d "$tempDir" ]; then
    mkdir $tempDir
fi


################################################################################
# Set index and chromosome sizes
################################################################################

INDEX="/home/groups/CEDAR/roskamsh/projects/Omics-QC-Pipeline-SE/data/genomes/GRCh38.primary_assembly.genome_bowtie1"
SIZES="/home/groups/CEDAR/roskamsh/projects/Omics-QC-Pipeline-SE/data/genomes/chrNameLength.txt"

[ -e "${INDEX}.1.ebwt" ] || echo >&2 "ERROR: No bowtie index files found for genome assembly ${assembly}!"
[ -e "${SIZES}" ] || echo >&2 "ERROR: No chromosome sizes file found for genome assembly ${assembly}!"

################################################################################
# Create tmp directory and start log file
################################################################################

# A temporary directory
mytemp=$(mktemp -d -p "${tempDir}")

# Generating temp log file to store bowtie stats
# logs
(
  echo -n "Start: "
  date
  echo "fastqBowtietoBigBed.sh"
  bowtie --version | head -n 1
  echo "$0 $@"
  echo
) > {mytemp}/LOG.log


################################################################################
# Set processor usage
# Either assigned by SGE or half of the available cores
# (at least 3 as we anyway need 3 later on)
################################################################################

if [ -z "$NSLOTS" ]; then
    USE=$(grep processor /proc/cpuinfo | wc -l | awk '{if($1>6){print int($1/2)}else{print 3}}')
else
    USE=$NSLOTS
fi


################################################################################
# Determine type of input - single end or paired end
################################################################################

b=( $(echo $fastqFile | tr "|" " ") )
if [ ${#b[@]} == 1 ]; then
    type="single"
elif [ ${#b[@]} == 2 ]; then
    type="paired"
    fastqFile1=${b[0]}
    fastqFile2=${b[1]}
fi

echo $type
# For NexSeq, you receive 4 fastq files, need to combine them prior to analysis
# This separation of se and pe would need to be done at the trimming step, extraction step


################################################################################
# Trim reads to specified length
################################################################################

if [ "$type" == "single" ]; then
    zcat $fastqFile | fastx_trimmer -Q33 -l $mapLen -o ${mytemp}/${sample}.trimmed.fq
elif [ "$type" == "paired" ]; then
    zcat $fastqFile1 | fastx_trimmer -Q33 -l $mapLen -o ${mytemp}/${sample}.read_1.trimmed.fq
    zcat $fastqFile2 | fastx_trimmer -Q33 -l $mapLen -o ${mytemp}/${sample}.read_2.trimmed.fq
fi


################################################################################
# Map with bowtie
################################################################################

if [ "$type" == "single" ]; then
    bowtie -p $USE -q -v $nrMM -m 1 --best --strata $INDEX ${mytemp}/${sample}.trimmed.fq > ${mytemp}/${sample}.reads.bowtie 2> ${outDir}${sample}.mapping.stat
    rm ${mytemp}/${sample}.trimmed.fq
elif [ "$type" == "paired" ]; then
    bowtie -p $USE -q -X $gap -v $nrMM -m 1 --best --strata $INDEX -1 ${mytemp}/${sample}.read_1.trimmed.fq -2 ${mytemp}/${sample}.read_2.trimmed.fq > ${mytemp}/${sample}.reads.bowtie 2> ${outDir}${sample}.mapping.stat
    rm ${mytemp}/${sample}.read_1.trimmed.fq ${mytemp}/${sample}.read_2.trimmed.fq
fi
#-q ensures fastq input/output


################################################################################
# Convert output to BED and sort
################################################################################

if [ "$type" == "single" ]; then

    cat ${mytemp}/${sample}.reads.bowtie | \
        awk -vFS="\t" -vOFS="\t" -vR=$mapLen '
        {
            print $3,$4,$4+R,".",0,$2
        }' | \
    sort -k1,1 -k2,2n -k6,6 -S 1G --parallel $(echo $USE | awk '{print $1-2}') > ${mytemp}/${sample}.reads.sorted.bed

elif [ "$type" == "paired" ]; then

    cat ${mytemp}/${sample}.reads.bowtie | \
    awk -vFS="\t" -vOFS="\t" '{
        if(NR%2==1)
        {
            split($1,T,"/")
            R[T[2]]=T[1]

            if (T[2]==1)
            {
                S="+"
            }
            else
            {
                S="-"
            }

            B=$4

        }
        else
        {
            print $3,B,$4+length($5),".",0,S

        }
    }' | \
    sort -k1,1 -k2,2n -k6,6 -S 1G --parallel $(echo $USE | awk '{print $1-2}') > ${mytemp}/${sample}.reads.sorted.bed

fi

rm ${mytemp}/${sample}.reads.bowtie

################################################################################
# Apply filters to collapse reads and/or remove information:
#
#   (0) Report all mapped reads
#   (1) Report fragments collapsed on identical chromosome, start, end,
#       and strand information
#   (A) Report both all mapped reads AND unique fragments,
#       i.e. report results from both filters 0 and 1
#
# Output is saved in $TMP/reads.filtered.bed
################################################################################

# Removing duplicate reads by defining if more than 1 read has the same chr, start and end, it is the same, and only keeping one of those in further analysis
(echo -n "Applying filters: "; date) >> $TMP/LOG.log

if [ $FILTER = "0" ]; then
    awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '{print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' ${mytemp}/${sample}.reads.sorted.bed > ${mytemp}/${sample}.reads.filtered.bed
elif [ $FILTER = "1" ]; then
    awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '(!x[$1" "$2" "$3" "$6]++){print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' ${mytemp}/${sample}.reads.sorted.bed > ${mytemp}/${sample}.reads.filtered.bed
elif [ $FILTER = "A" ]; then
# run filter 0
    awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '{print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' ${mytemp}/${sample}.reads.sorted.bed > ${mytemp}/${sample}.reads.filtered.0.bed
# run filter 1
    awk -vOFS='\t' -vS=$KEEPSEQ -vM=$KEEPMM '(!x[$1" "$2" "$3" "$6]++){print $1,$2,$3,S==1?$4:".",M==1?$5:"0",$6}' ${mytemp}/${sample}.reads.sorted.bed > ${mytemp}/${sample}.reads.filtered.1.bed
fi
#Keepseq = keep sequence of read
#KeepMM = keep the # of mismatches
# filter=0 line S==1?$4:"." means IF the first thing is not true, change the first thing to the second thing - i.e. is S is not equal to 1, change S=4
# filter=1 line !x[$1" "$2" "$3" "$6]++ means: do not keep reads with the same chr, start and name --> remove duplicates


# Clean up
rm ${mytemp}/${sample}.reads.sorted.bed

# Add if/else to collapse reads, like is done in bowtie_se.sh with FILTER=0,1,A option
# Don't collapse paired-end read fragments unless they have same sequence, start and end


################################################################################
# Convert BED to BigBed format and copy to output directory
################################################################################

if [ $FILTER = "A" ]; then

    (echo -n "Running bedToBigBed for filter A-0: "; date) >> ${mytemp}/LOG.log
    bedToBigBed ${mytemp}/${sample}.reads.filtered.0.bed $SIZES ${mytemp}/${sample}.all.bb 2>>${mytemp}/LOG.log
    echo >> ${mytemp}/LOG.log

    (echo -n "Running bedToBigBed for filter A-2: "; date) >> ${mytemp}/LOG.log
    bedToBigBed ${mytemp}/${sample}.reads.filtered.1.bed $SIZES ${mytemp}/${sample}.bb 2>>${mytemp}/LOG.log

    rm ${mytemp}/${sample}.reads.filtered.{0,1}.bed

else
(echo -n "Running bedToBigBed: "; date) >> ${mytemp}/LOG.log
    bedToBigBed ${mytemp}/${sample}.reads.filtered.bed $SIZES ${mytemp}/${sample}.reads.final.bb
    rm ${mytemp}/${sample}.reads.filtered.bed
fi

## Moving final product bigBed file to your output directory
if [ "$type" == "single" && "$FILTER" == "A"]; then
    mv ${mytemp}/${sample}.all.bb ${outDir}/${sample}.mapped.reads.all.bb
    mv ${mytemp}/${sample}.bb ${outDir}/${sample}.mapped.reads.bb
elif [ "$type" == "single" && ("$FILTER" == "0" || "$FILTER" == "1")]; then
    mv ${mytemp}/${sample}.reads.final.bb ${outDir}/${sample}.mapped.reads.bb
elif [ "$type" == "paired" && "$FILTER" == "A"]; then
    mv ${mytemp}/${sample}.all.bb ${outDir}/${sample}.mapped.fragments.all.bb
    mv ${mytemp}/${sample}.bb ${outDir}/${sample}.mapped.fragments.bb
elif [ "$type" == "paired" && ("$FILTER" == "0" || "$FILTER" == "1")]; then
    mv ${mytemp}/${sample}.reads.final.bb ${outDir}/${sample}.mapped.fragments.bb
fi

## Remove temporary bigbed files
if [ "$FILTER" == "A"]; then
    rm ${mytemp}/${sample}.all.bb
    rm ${mytemp}/${sample}.bb
else
    rm ${mytemp}/${sample}.reads.final.bb

(echo -n "End: "
date ) >> $TMP/LOG.log

if [ $KEEPLOG = "1" ]; then
    mv ${mytemp}/LOG.log ${sample}.log
fi
# Storing all stats/log information from bedToBigBed in logfile named the same as our outfile


################################################################################
# Clean temp folder and exit
################################################################################

rm -rf ${mytemp}

exit 0




