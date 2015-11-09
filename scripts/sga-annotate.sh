#!/bin/bash
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load parallel
module load gcc/4.8.1

function sgaannotate {
    VCFDIR=$1
    SAMPLE=$2
    NORMALBAM=$3
    TUMORBAM=$4
    OUTVCFDIR=$5

    VCF=${VCFDIR}/${SAMPLE}.vcf
    ANNOTATEDVCF=${OUTVCFDIR}/${SAMPLE}-annotated.vcf

    REFERENCE=/.mounts/labs/simpsonlab/data/references/hs37d5.fa 
    SGABIN=/.mounts/labs/simpsonlab/users/jsimpson/code/sga/src/build/SGA/sga 
    echo "${SGABIN} somatic-variant-filters --annotate-only -t 1 --tumor $TUMORBAM --normal $NORMALBAM --reference $REFERENCE $VCF out to  $ANNOTATEDVCF"
    module load gcc/4.8.1
    ${SGABIN} somatic-variant-filters --annotate-only -t 1 \
            --tumor $TUMORBAM --normal $NORMALBAM --reference $REFERENCE $VCF > $ANNOTATEDVCF
}

export -f sgaannotate

if [ $# -eq 0 ]
then
    echo "$0 - annotates generates list of target variants for list of samples given"
    echo "Usage: $0 vcfindir vcfoutdir [samplesfile]"
    exit 
fi

if [ -z "$1" ]
then
    echo "No vcf indir supplied"
    exit
fi

if [ -z "$2" ]
then
    echo "No vcf out dir supplied"
    exit
fi

VCFDIR=$1
OUTVCFDIR=$2

BAMDIR=xfer.genome.wustl.edu/gxfer1/67722228712318/

SAMPLESFILE=${3-${BAMDIR}/sample-directory.txt}

JOBS=8

parallel -j $JOBS --colsep '\t' sgaannotate ${VCFDIR} {1} ${BAMDIR}{2} ${BAMDIR}{3} ${OUTVCFDIR} :::: $SAMPLESFILE
