#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

module load java/1.7.0_79
module load gatk/3.4.0

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ -z "$5" ] \
    || [ ! -f "$1" ] || [ ! -f "$2" ] || [ ! -f "$3" ] || [ ! -f "$4" ]
then
    echo "$0 - performs GATK Hapolotypecaller annotation"
    echo "Usage: $0 INPUT.VCF NORMAL.BAM TUMOUR.BAM INTERVALS.LIST OUTPUT.VCF [REFERENCE=$DEFREF]"
    echo "invocation was: $0 $1 $2 $3 $4 $5 $6"
    exit 
fi

INPUT_VCF=$1
NORMAL_BAM=$2
TUMOUR_BAM=$3
INTERVALS=$4
OUTPUT_VCF=$5
REFERENCE=${6-$DEFREF}

JMEM=8g

if [ ! -f "$OUTPUT_VCF" ]
then
    java -Xmx${JMEM} -jar ${GATKROOT}/GenomeAnalysisTK.jar \
        -R ${REFERENCE} \
        -T HaplotypeCaller \
        -I ${NORMAL_BAM} \
        -I ${TUMOUR_BAM} \
        --alleles ${INPUT_VCF} \
        -L ${INTERVALS} \
        --genotyping_mode GENOTYPE_GIVEN_ALLELES \
        --activeRegionSize 500 \
        --forceActive \
        -dt none \
        --minPruning 100 \
        -nct 4 \
        --output_mode EMIT_ALL_SITES \
        -o ${OUTPUT_VCF}
fi
