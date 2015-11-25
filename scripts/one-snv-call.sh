#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles

module load python
module load gcc/4.8.1
module use /.mounts/labs/simpsonlab/modules/
module load openblas
module load python-packages/2

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ ! -f "$1" ] || [ ! -f "$2" ] || [ ! -f "$3" ]
then
    echo "$0 - performs call validation based on readcounts files"
    echo "Usage: $0 INPUT.VCF NORMAL_READCOUNTS.RC TUMOUR_READCOUNTS.RC OUTPUT.VCF"
    echo "$0 $1 $2 $3 $4"
    exit 
fi

INPUT_VCF=$1
NORMAL_RC=$2
TUMOUR_RC=$3
OUTPUT_VCF=$4

if [ ! -f "$OUTPUT_VCF" ]
then
    ./scripts/snv_readcounts.py ${INPUT_VCF} ${NORMAL_RC} ${TUMOUR_RC} | \
        ./scripts/snv_indel_call.py --error 0.01 --callthreshold 0.02 --strandbias -1 --germlineprob 0.01 --mindepth 25 \
        | sort -k1,1 -k2,2n  > ${OUTPUT_VCF}
fi
