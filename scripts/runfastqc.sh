#!/bin/bash
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load parallel

function dofastqc {
    PAIR1=$1
    FASTQCDIR=$2

    PAIR2=${PAIR1/_1/_2}
    SAMPLE=$( dirname $PAIR1 )
    SAMPLE=$( basename $SAMPLE )
    module load java/1.8.0_40
    module load fastqc/0.11.2
    module load R/3.2.1

    mkdir -p ${FASTQCDIR}/${SAMPLE}
    echo "running: fastqc -o $FASTQCDIR/${SAMPLE} $PAIR1 $PAIR2"
    fastqc -o $FASTQCDIR/${SAMPLE} $PAIR1 $PAIR2
}

export -f dofastqc

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ]
then
    echo "$0 - runs fastqc on a directory of fastqs"
    echo "Usage: $0 fastqdir fastqcdir [njobs=8]"
    exit 
fi

FASTQDIR=$1
FASTQCDIR=$2
JOBS=${3-8}

mkdir -p $FASTQCDIR

find ${FASTQDIR} -name "*_1.fastq*" | parallel -j $JOBS dofastqc {} ${FASTQCDIR} 
