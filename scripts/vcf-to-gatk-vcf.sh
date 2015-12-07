#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ ! -f "$1" ] 
then
    echo "$0 - converts VCF to something GATK will ingest"
    echo "Usage: $0 INPUT.VCF OUTPUT.VCF"
    echo "invocation was: $0 $1 $2"
    exit 
fi

INPUT_VCF=$1
OUTPUT_VCF=$2

if [ ! -f ${OUTPUT_VCF} ]
then
    grep "^#" ${INPUT_VCF} > ${OUTPUT_VCF}
    grep "^[c]*[h]*[r]*[0-9][0-9]*" ${INPUT_VCF} | sed -e 's/^chr//' | sort -k1,1n -k2,2n >> ${OUTPUT_VCF}
    grep "^[c]*[h]*[r]*[XY]" ${INPUT_VCF} | sed -e 's/^chr//' | sort -k1,1 -k2,2n >> ${OUTPUT_VCF}
    grep "^MT" ${INPUT_VCF} | sort -k1,1 -k2,2n >> ${OUTPUT_VCF}
    grep "^GL" ${INPUT_VCF} | sort -k1,1 -k2,2n >> ${OUTPUT_VCF}
    grep "^hs" ${INPUT_VCF} | sort -k1,1 -k2,2n >> ${OUTPUT_VCF}
fi
