#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

module load java/1.7.0_79
module load gatk/3.4.0

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] \
    || [ ! -f "$1" ] || [ ! -f "$2" ] || [ ! -f "$3" ]
then
    echo "$0 - performs realignment around indels"
    echo "Usage: $0 INPUT.BAM GERMLINE_INDELS.VCF TARGETS OUTPUT.BAM [REFERENCE=$DEFREF]"
    exit 
fi


BAM=$1
VCF=$2
TARGETS=$3
OUTBAM=$4
REFERENCE=${5-$DEFREF}

JMEM=6g
if [ ! -f ${OUTBAM} ]
then
    java -Xmx${JMEM} -jar ${GATKROOT}/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R ${REFERENCE} \
        -targetIntervals ${TARGETS} \
        -I ${BAM} \
        -known ${VCF} \
        -o ${OUTBAM}
fi
