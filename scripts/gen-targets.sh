#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

module load java/1.7.0_79
module load gatk/3.4.0

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ ! -f "$1" ] 
then
    echo "$0 - generates realign targets around indels"
    echo "Usage: $0 GERMLINE_INDELS.VCF OUTPUT.TARGETS [REFERENCE=$DEFREF]"
    exit 
fi


VCF=$1
OUTTARGETS=$2
REFERENCE=${3-$DEFREF}

JMEM=6g
if [ ! -f ${OUTTARGETS} ]
then
    java -Xmx${JMEM} -jar ${GATKROOT}/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R ${REFERENCE} \
        --known ${VCF} \
        -o ${OUTTARGETS}
fi
