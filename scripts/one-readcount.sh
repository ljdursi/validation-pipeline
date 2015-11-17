#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

module load gcc/4.8.1
module load bam-readcount

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3"] || [ ! -f "$1" ] || [ ! -f "$2" ] 
then
    echo "$0 - performs read counts against a BAM at target locations"
    echo "Usage: $0 INPUT.BAM TARGETS.BED OUTPUT.RC [REFERENCE=$DEFREF]"
    echo "$0 $1 $2 $3"
    exit 
fi

BAM=$1
TARGETS=$2
OUTPUT=$3
REFERENCE=${4-$DEFREF}

if [ ! -f ${OUTPUT} ]
then
    bam-readcount --reference-fasta $REFERENCE --site-list <( awk '{ printf "%s\t%d\t%d\n",$1,$2+1,$3+1 }' $TARGETS ) \
        --max-count 8000 --max-warnings 0 $BAM > ${OUTPUT}
fi
