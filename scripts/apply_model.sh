#!/bin/bash -l
#
# Calls an R script to apply a model given a model file,
# an input, an output, and optionally a threshold

module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module load R/3.2.5

readonly DEFTHRESH=0.5

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ ! -f "$1" ] || [ ! -f "$2" ]
then
    echo "$0 - Apply an ensemble model (provided) to the VCF"
    echo "Usage: $0 modelfile input.vcf output.vcf [thresh=${DEFTHRESH}]"
    exit 
fi

readonly MODEL=$1
readonly INPUTVCF=$2
readonly OUTPUTVCF=$3
readonly THRESH=${4:-${DEFTHRESH}}

Rscript --vanilla scripts/filter_calls_by_model.R $MODEL $INPUTVCF $OUTPUTVCF $THRESH \
    && mv ${OUTPUTVCF}.bgz ${OUTPUTVCF}.gz \
    && mv ${OUTPUTVCF}.bgz.tbi ${OUTPUTVCF}.gz.tbi
