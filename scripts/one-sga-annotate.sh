#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load parallel
module load gcc/4.8.1

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa.gz
SGABIN=/.mounts/labs/simpsonlab/users/jsimpson/code/sga/src/build/SGA/sga 

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] \
    || [ ! -f "$1" ] || [ ! -f "$2" ] || [ ! -f "$3" ]
then
    echo "$0 - performs call annotation based on bam files"
    echo "Usage: $0 INPUT.VCF NORMAL.BAM TUMOUR.BAM OUTPUT.VCF [REFERENCE=$DEFREF]"
    echo "invocation was: $0 $1 $2 $3 $4 $5"
    exit 
fi

INPUT_VCF=$1
NORMAL_BAM=$2
TUMOUR_BAM=$3
OUTPUT_VCF=$4
REFERENCE=${5-$DEFREF}

if [ ! -f "$OUTPUT_VCF" ]
then
    ${SGABIN} somatic-variant-filters --annotate-only -t 4 \
        --tumor $TUMOR_BAM --normal $NORMAL_BAM --reference $REFERENCE $INPUT_VCF > $OUTPUT_VCF
fi
