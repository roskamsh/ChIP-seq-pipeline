#!/bin/bash
#set -o errexit
#set -o pipefail

################################################################################
# Help
################################################################################

if [ $# -eq 0 ]; then
echo >&2 "
$(basename $0) - Calling broad or narrow peaks, dependig on the input.

USAGE: $(basename $0) -t <treatment bed> -c <control bed> -n <outfile basename> -p <peak call> -a <genome assembly> [OPTIONS]
-t    Treatment bed file                                                                     [required]
-c    Control bed file                                                                       [required]
-n    Basename of the output file                                                            [required]
-p    Type of peak calling (narrow or broad)                                                 [required]
-a    Genome assembly (hg38, hg19, dm3, etc.)                                                [required]
"
exit 1
fi


################################################################################
# Parse input and error checking
################################################################################

while getopts "t:c:n:p:a:" op
do
    case "$op" in
        t)  treated="$OPTARG";;
        c)  control="$OPTARG";;
        n)  sample="$OPTARG";;
        p)  peak="$OPTARG";;
        a)  assembly="$OPTARG";;
        \?) exit 1;;
    esac
done

if [ -z "$treated" ]; then
    echo >&2 "ERROR: Specify treatment bed file!"
    exit 1
fi
if [ -z "$control" ]; then
    echo >&2 "ERROR: Specify control bed file!"
    exit 1
fi
if [ -z "$sample" ]; then
    echo >&2 "ERROR: Specify basename of output file!"
    exit 1
fi
if [ -z "$peak" ]; then
    echo >&2 "ERROR: Specify the type of peak calling!"
    exit 1
fi
if [ -z "$assembly" ]; then
    echo >&2 "ERROR: Specify genome assembly!"
    exit 1
fi

################################################################################
# Create tmp directory and start log file
################################################################################

# Grab current working directory to input into peak calling commands
WDIR=$PWD

# A temporary directory
TMP=$(mktemp -d)
cd $TMP

# Generating temp log file to store bowtie stats
# logs
(
  echo -n "Start: "
  date
  echo "callPeaks.sh"
) > $TMP/LOG.log


##############################################################################
# Determine required parameters for peak calling
##############################################################################

if [ "$assembly" == "hg38" ] || [ "$assembly" == "hg19" ]; then
    gsize="2.7e9"
elif [ "$assembly" == "mm9" ] || [ "$assembly" == "mm10" ]; then
    gsize="1.7e9"
elif [ "$assembly" == "ce6" ] || [ "$assembly" == "ce10" ]; then
    gsize="9e7"
elif [ "$assembly" == "dm3" ] || [ "$assembly" == "dm6" ]; then
    gsize="1.2e9"
fi

inputDir=$(dirname -- "$treated")
basenameTreat=$(basename -- "$treated")
basenameControl=$(basename -- "$control")

if [ "$peak" == "narrow" ]; then
    outDir="results/macs2/${sample}"
elif [ "$peak" == "broad" ]; then
    outDir="results/SICER/${sample}"
fi

###############################################################################
# Call peaks
###############################################################################

if [ "$peak" == "narrow" ]; then
    macs2 callpeak -t ${WDIR}/${treated} -c ${WDIR}/${control} -n $sample --outdir $TMP -g $gsize &> $TMP/LOG.log
elif [ "$peak" == "broad" ]; then
    SICER.sh ${WDIR}/${inputDir} $basenameTreat $basenameControl $TMP $assembly 1 200 150 0.8 600 1e-8 &> $TMP/LOG.log
fi

(echo -n "End: "
date ) >> $TMP/LOG.log

##############################################################################
# Move output files to output directory and remove TMP directory
#############################################################################
mv $TMP/LOG.log ${TMP}/${sample}.log
mv $TMP/* ${WDIR}/${outDir}
rm -rf $TMP

# Exit
exit 0
