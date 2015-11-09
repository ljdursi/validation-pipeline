#!/bin/bash 
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

module load samtools/1.2
module load bcftools/1.1

DEFREF=/.mounts/labs/simpsonlab/data/references/hs37d5.fa 

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]
then
    echo "$0 - performs variant calling on two BAMS"
    echo "Usage: $0 NORMAL.BAM TUMOUR.BAM REGIONS.BED TARGETS.BED [REFERENCE=$DEFREF]
    fastqdir bamdir [njobs=8] [ref=$DEFREF]"
    exit 
fi

NORMAL=$1
TUMOUR=$2
REGIONS=$3
TARGETS=$4
REFERENCE=${5-$DEFREF}

samtools mpileup -C50 --positions $REGIONS -gf $REFERENCE $1 $2 \
    | bcftools call --multiallelic-caller --variants-only --targets-file $TARGETS --> ${3}.vcf

