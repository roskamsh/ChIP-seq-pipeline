#!/bin/bash
set -o errexit
set -o pipefail

################################################################################
# Requirements
################################################################################

# Programs:
# * Jim Kent's utils (http://hgwdev.cse.ucsc.edu/~kent/src/)
# * bedtools (http://code.google.com/p/bedtools/)


# Files:
# * /pathTo/${ASSEMBLY}.chrom.sizes

################################################################################
# Set default values
################################################################################

ASSEMBLY="hg38" # Genome assembly
OUTFILE=""     # Outfile, if not supplied ${input%.bb}.bw
EXT="0"        # Fragment length of the experiment
EXT2="0"       # Addititional extension
NORM="1"       # Normalize read coverge to 1 million mapped reads

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
  echo >&2 "
$(basename $0) - Compute a BigWig file from a BigBed file of mapped reads

USAGE: $(basename $0) -i <BigBed input file> [OPTIONS]
 -i     Input file (BigBed) [required]
 -o     Output file (BigWig) [default: <input file>.bw]
 -g     Genome assembly (e.g. dm3, hg19) [default: $ASSEMBLY]
 -e     Extend mapped read coordinates to total length of n (off: 0) [default: $EXT]
        Set to 0 for paired-end data
        Set to fragment size (e.g. 600) for single-end data to extend fragment coordinates
 -E     Additional extension [default: $EXT2]
        Enables the artificial extension of reads.
 -n     Normalize the coverage to 1 million mapped reads (0/1) [default: $NORM]

NOTES:
The program requires a file with the chromosome sizes of your assembly stored 
here /groups/stark/genomes/chrom/ASSEMBLY.chrom.sizes.

The program uses up to 4 CPU cores, and about 2GB of RAM for a 100MB BigBed file.
"
  exit 1
fi

################################################################################
# Parse input
################################################################################

while getopts "i:o:g:e:n:E:" o
do
  case "$o" in
    i) INFILE="$OPTARG";;
    o) OUTFILE="$OPTARG";;
    g) ASSEMBLY="$OPTARG";;
    e) EXT="$OPTARG";;
    E) EXT2="$OPTARG";;
    n) NORM="$OPTARG";;
   \?) exit 1;;
  esac
done

################################################################################
# Set chromosome size file
################################################################################

# Set $INDEX and $SIZES
INDEX=/home/groups/CEDAR/tools/indices/bowtie/${ASSEMBLY}/${ASSEMBLY}
SIZES=/home/groups/CEDAR/tools/genomes/chrom/${ASSEMBLY}.chrom.sizes


# Throw error message if chromosome sizes file does not exist
if [ ! -e "$SIZES" ]; then
  echo >&2 "ERROR: No chromosome size file found for genome assembly ${ASSEMBLY}!"
  exit 1
fi

################################################################################
# Run main program
################################################################################

# Get number of total fragments
N=$(bigBedInfo $INFILE | awk '(/^itemCount/){gsub(/,/,"",$2);print $2}')

# Get mapped read length and compute extension paramter E
# Set to zero if EXT should be zero
E=$(bigBedToBed $INFILE stdout -maxItems=1 | awk -v E=$EXT 'NR==1{if(E==0){print "0"}else{print E-($3-$2)}}')

# Error if extension lower than 0 
if [ "$E" -lt "0" ]; then
  echo >&2 "ERROR: Total fragment length should be larger than read length!" 
  exit 1
fi

# Set outfile to infile.bw or to the user given name
if [ "$OUTFILE" = "" ]; then
  OUTFILE=${INFILE%.bb}.bw
else
  OUTFILE=$OUTFILE
fi

# Compute coverage and get BedGraph file
# Convert BedGraph file to BigBed
# If required extend the reads prior to computing the coverage
# EXT2 is and additional extension parameter
/home/groups/CEDAR/tools/kentUtils/bigBedToBed $INFILE stdout | \
  if [ "$E" = "0" -a "$EXT2" = "0" ]; then
    cat  
  else
    awk -vE=$E -vE2=$EXT2 -vC=$SIZES -vOFS="\t" '
      BEGIN {while(getline<C){chr[$1]=$2}} 
      { 
        if($6=="+"){$3=$3+E+E2}else{$2=$2-E-E2}
        if($2<0){$2=0}
        if($3>chr[$1]){$3=chr[$1]}
        print $0
      }'  
  fi | \
  genomeCoverageBed -i stdin -bg -g $SIZES | \
    if [ "$NORM" = "1" ]; then
      awk -vN=$N -vOFS="\t" '{print $1,$2,$3,1e6*$4/N}'
    else
      awk -vN=$N -vOFS="\t" '{print $1,$2,$3,$4}'
    fi | \
  /home/groups/CEDAR/tools/kentUtils/wigToBigWig stdin $SIZES $OUTFILE

# Exit
exit 0

