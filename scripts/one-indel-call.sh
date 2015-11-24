#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles

module load python
module load gcc/4.8.1
module use /.mounts/labs/simpsonlab/modules/
module load openblas
module load python-packages/2

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ ! -f "$1" ] 
then
    echo "$0 - performs call validation based on readcounts files"
    echo "Usage: $0 VALIDATION_VCF_WITH_COUNTS.VCF OUTPUT.VCF"
    echo "$0 $1 $2"
    exit 
fi

INPUT_VCF=$1
OUTPUT_VCF=$2

if [ ! -f "$OUTPUT_VCF" ]
then
    ./scripts/snv_indel_call.py --error 0.01 --callthreshold 0.02 --mindepth 25 --strandbias -1 \
        -i ${INPUT_VCF} \
        | sort -k1,1 -k2,2n | sed -e 's/^\([1-9XYh]\)/chr\1/' > ${OUTPUT_VCF}
fi
