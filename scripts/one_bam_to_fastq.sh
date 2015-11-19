#!/bin/bash -l 
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load samtools/1.2
module load java/1.8.0_40
module load picard/1.92
module load parallel

DEFNJOBS=8

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ ! -f $1 ] 
then
    echo "$0 - converts BWA ALN aligned BAMS to FASTQ"
    echo "     BAMS are separated by RG so that they can be recombined with RG information "
    echo "Usage: $0 bamindir fastqoutdir [njobs=${DEFNJOBS}]"
    echo "invocation was: $0 $1 $2 $3"
    exit 
fi

BAMFILE=$1
FASTQDIR=$2
NJOBS=${3-8}
MEM=4g

BASE=${FASTQDIR}/$( basename $BAMFILE .bam )
if [ ! -d ${BASE} ]
then
    echo "Converting ${BAMFILE} to fastq into ${FASTQDIR}/${BASE}"
    mkdir ${BASE}
    samtools view -H $BAMFILE > ${BASE}/bam-header.txt
    echo "java -Xmx${MEM} -Xms${MEM} -jar ${PICARDROOT}/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT \
            OUTPUT_PER_RG=true OUTPUT_DIR=${BASE} INPUT=${BAMFILE}"
    java -Xmx${MEM} -Xms${MEM} -jar ${PICARDROOT}/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT \
            OUTPUT_PER_RG=true OUTPUT_DIR=${BASE} INPUT=${BAMFILE} &&
    for file in ${BASE}/*.fastq; do gzip $file; done
fi  
