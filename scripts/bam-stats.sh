#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

module load java/1.7.0_79
module load samtools/1.2
module load bcftools/1.1
module load picard/1.92

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa.gz

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ]
then
    echo "$0 - performs samtools flagstat and picard CollectAlignmentSummaryMetrics stats"
    echo "     on a bam file"
    echo "Usage: $0 FILE.BAM OUTPUTDIR [REFERENCE=$DEFREF]"
    exit 
fi

BAM=$1
OUTDIR=$2
REFERENCE=${3-$DEFREF}


#samtools

base=$( basename $BAM .bam )
echo "samtools flagstat $BAM > ${OUTDIR}/${base}_flagstat_metrics.txt"
samtools flagstat $BAM > ${OUTDIR}/${base}_flagstat_metrics.txt

#picard
JMEM=4g
java -Xmx${JMEM} -Xms${JMEM}  -jar ${PICARDROOT}/CollectAlignmentSummaryMetrics.jar \
                METRIC_ACCUMULATION_LEVEL=ALL_READS \
                INPUT=${BAM} OUTPUT=${OUTDIR}/${base}_picard_alignment_metrics.txt \
                REFERENCE_SEQUENCE=${REFERENCE} \
                VALIDATION_STRINGENCY=LENIENT
