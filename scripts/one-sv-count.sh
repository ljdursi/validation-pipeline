#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

module load python
module load python-packages/2
module load gcc/4.8.1
module load 

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3"] || [ -z "$4" ] \
    || [ ! -f "$1" ] || [ ! -f "$2" ]  || [ ! -f "$3" ]
then
    echo "$0 - performs read counts in support of a set of structural variant calls"
    echo "Usage: $0 INPUT.VCF NORMAL.BAM TUMOR.BAM OUTPUT.VCF"
    echo "invocation was: $0 $1 $2 $3 $4"
    exit 
fi

INPUT_VCF=$1
NORMAL_BAM=$2
TUMOUR_BAM=$3
OUTPUT_VCF=$4

if [ ! -f ${OUTPUT_VCF} ]
then
    ./scripts/count-sv-support.py ${INPUT_VCF} ${NORMAL_BAM} ${TUMOUR_BAM} \
        | ./scripts/filter_by_variant_len.py -m 300 -o ${OUTPUT_VCF}
fi
