#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles

module load bcftools/1.1

module load python
module load gcc/4.8.1
module use /.mounts/labs/simpsonlab/modules/
module load openblas
module load python-packages/2

if [ $# -eq 0 ] || [ -z "$1"] || [ -z "$2"] || [ -z "$3"] || [ ! -f "$1" ] || [ ! -f "$2" ] 
then
    echo "$0 - performs variant calling on an mpileup file"
    echo "Usage: $0 INPUT.MPILEUP TARGETS.GZ OUTPUT.VCF"
    exit 
fi

MPILEUPFILE=$1
TARGETS=$2
VCFFILE=$3

bcftools call -c --output-type v --targets-file ${TARGETS} ${MPILEUPFILE} | \
    sed -e 's/",Version="[0-9.]*">/">/' | \
    ./scripts/genotypes.py --error 0.01 --callthreshold 0.02 --strandbias 0.1 --mindepth 10 \
    > $VCFFILE
