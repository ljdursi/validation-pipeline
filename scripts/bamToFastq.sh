#!/bin/bash
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load parallel

function picardbamtofastq {
    FASTQDIR=$1
    BAMFILE=$2
    MEM=${3-4g}

    BASE=${FASTQDIR}/$( basename $BAMFILE .bam )

    if [ ! -d ${BASE} ]
    then
        module load samtools/1.2
        module load java/1.8.0_40
        module load picard/1.92
        mkdir ${BASE}
        samtools view -H $BAMFILE > ${BASE}/bam-header.txt
        echo "java -Xmx${MEM} -Xms${MEM} -jar ${PICARDROOT}/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT \
                OUTPUT_PER_RG=true OUTPUT_DIR=${BASE} INPUT=${BAMFILE}"
        java -Xmx${MEM} -Xms${MEM} -jar ${PICARDROOT}/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT \
                OUTPUT_PER_RG=true OUTPUT_DIR=${BASE} INPUT=${BAMFILE} &&
            for file in ${BASE}/*.fastq; do gzip $file; done
    fi  
}

export -f picardbamtofastq

if [ $# -eq 0 ]
then
    echo "$0 - converts BWA ALN aligned BAMS to FASTQ"
    echo "Usage: $0 bamindir fastqoutdir [njobs=8]"
    exit 
fi

if [ -z "$1" ]
then
    echo "No BAM input dir supplied"
    exit
fi

if [ -z "$2" ]
then
    echo "No FASTQ output dir supplied"
    exit
fi

BAMDIR=$1
FASTQDIR=$2

JOBS=${3-8}

find ${BAMDIR} -name *.bam -print | parallel -j $JOBS picardbamtofastq ${FASTQDIR} {} 4g
