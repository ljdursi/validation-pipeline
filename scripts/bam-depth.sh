#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

module load java/1.7.0_79
module load gatk/3.4.0

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]
then
    echo "$0 - calculates depth of coverage in target regions for a given BAM file"
    echo "Usage: $0 FILE.BAM REGIONS.BED OUTPUTDIR [SUFFIX=''] [REFERENCE=$DEFREF]"
    exit 
fi

BAM=$1
BED=$2
OUTDIR=$3
SUFFIX=$4
REFERENCE=${5-$DEFREF}

base=$( basename $BAM .bam )
JMEM=5g
java -Xmx${JMEM} -Xms${JMEM} -jar ${GATKROOT}/GenomeAnalysisTK.jar \
                -T DepthOfCoverage \
                -R ${REFERENCE} \
                -I ${BAM} \
                -L ${BED} \
                -o ${OUTDIR}/${base}_target_coverage${SUFFIX} \
                -ct 1 -ct 10 -ct 20
