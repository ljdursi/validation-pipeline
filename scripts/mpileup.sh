#!/bin/bash
module purge
module load samtools/1.2
module load bcftools/1.1
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load parallel

JOBS=8
BAMDIR=xfer.genome.wustl.edu/gxfer1/67722228712318/
SAMPLESFILE=${BAMDIR}/sample-directory.txt

function callmpileup {
    REFERENCE=/.mounts/labs/simpsonlab/data/references/hs37d5.fa 
    echo "samtools mpileup -gf $REFERENCE $1 $2 | bcftools call --multiallelic-caller --variants-only > ${3}.vcf"
    samtools mpileup -gf $REFERENCE $1 $2 | bcftools call --multiallelic-caller --variants-only > ${3}.vcf
}

export -f callmpileup

parallel --colsep '\t' callmpileup ${BAMDIR}{2} ${BAMDIR}{3} vcfs/{1} :::: $SAMPLESFILE
