#!/bin/bash
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load parallel

function bwamem {
    rg=$1
    reference=$2
    pair1=$3
    pair2=$4
    outbam=$5

    module load bwa/0.7.12
    module load samtools/0.1.19
    bwa mem -T 0 -R "${rg}" $reference $pair1 $pair2 | samtools view -bS - > ${outbam}
}

function alignfastqs {
    SAMPLEBASEDIR=$1
    FASTQDIR=$2
    BAMDIR=$3
    REFERENCE=$4

    if [ "$FASTQDIR" != "$SAMPLEBASEDIR" ]
    then
        BAMSUBDIR=${BAMDIR}/${SAMPLEBASEDIR#$FASTQDIR}

        HEADERFILE=${SAMPLEBASEDIR}/bam-header.txt
        if [ ! -d ${BAMSUBDIR} ]
        then
            mkdir -p ${BAMSUBDIR}
        fi    

        for pair1 in ${SAMPLEBASEDIR}/*_1.fastq*
        do
            pair2=${pair1/_1/_2}
            stripgz=${pair1%.gz}
            base=${stripgz%_1.fastq}
            base=$(basename $base)

            if [ ! -f ${BAMSUBDIR}/${base}.bam ]
            then
                rg=$( grep $base $HEADERFILE )
                echo "bwamem \"$rg\" $REFERENCE $pair1 $pair2 ${BAMSUBDIR}/${base}.bam "
                bwamem "$rg" $REFERENCE $pair1 $pair2 ${BAMSUBDIR}/${base}.bam 
            fi

        done
        module load picard/1.92
        mem=4g
        inputs=$( ls ${BAMSUBDIR}/*.bam | xargs -n1 -I{} echo "INPUT="{} )
        java -Xmx${mem} -Xms${mem} -jar ${PICARDROOT}/MergeSamFiles.jar \
                $inputs OUTPUT=${BAMSUBDIR}.bam VALIDATION_STRINGENCY=LENIENT MERGE_SEQUENCE_DICTIONARIES=true

        module load samtools/0.1.19
        samtools index ${BAMSUBDIR}.bam
    fi
}

export -f alignfastqs
export -f bwamem

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa.gz

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ]
then
    echo "$0 - Aligns FASTQ files"
    echo "Usage: $0 fastqdir bamdir [njobs=8] [ref=$DEFREF]"
    exit 
fi

FASTQDIR=$1
BAMDIR=$2
JOBS=${3-8}
REFERENCE=${4-$DEFREF}

find ${FASTQDIR} -type d -print | parallel -j $JOBS alignfastqs {} ${FASTQDIR} ${BAMDIR} ${REFERENCE}
