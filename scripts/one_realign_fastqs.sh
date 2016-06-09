#!/bin/bash -l
#
# realigns the FASTQ files that are in directory FASTQDIR
# (eg, bams/sample-id/*.fastq) into individual BAMs using
# bwa mem, uses Picard to merge the BAMs, and then samtools
# index to index it.

module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load samtools/0.1.19
module load parallel

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa.gz

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ]
then
    echo "$0 - Aligns FASTQ files"
    echo "Usage: $0 fastqdir output.bam [njobs=8] [ref=$DEFREF]"
    exit 
fi

FASTQDIR=$1
OUTPUTBAM=$2
JOBS=${3:-4}
REFERENCE=${4:-$DEFREF}

HEADERFILE=${FASTQDIR}/bam-header.txt

##
## BWA MEM align a particular FASTQ pair
## This is in a function so we can call gnu parallel over it
##
function bwamem {
    rg=$1
    reference=$2
    pair1=$3
    pair2=$4
    outbam=$5

    module load bwa/0.7.12
    >&2 echo "invocation is: rg=$1, ref=$2 p1=$3, p2=$4, outbam=$5"
    >&2 echo "bwa mem -T 0 -R _${rg}_ $reference $pair1 $pair2 | samtools view -bS - > ${outbam}"
    bwa mem -T 0 -R "${rg}" $reference $pair1 $pair2 | samtools view -bS - > ${outbam}
}

function call_bwamem {
    firstfastq=$1
    headerfile=$2
    bamdir=$3
    reference=$4

    pair=${firstfastq//.gz/}
    pair=${pair//_1.fastq/}
    base=$( basename $pair )

    rg=$( grep $base $headerfile )

    >&2 echo "invocation is: ffq=$1, header=$2 bamdir=$3, ref=$4"
    >&2 echo "bwamem _${rg}_ $reference ${pair}_1.fastq* ${pair}_2.fastq* ${bamdir}/${base}.bam"
    bwamem "$rg" $reference ${pair}_1.fastq* ${pair}_2.fastq* ${bamdir}/${base}.bam
}   

export -f bwamem
export -f call_bwamem

if [ ! -f "${OUTPUTBAM}" ]
then
    SUBBAMS_DIR=${OUTPUTBAM%%.bam}
    mkdir -p ${SUBBAMS_DIR}

    ls ${FASTQDIR}/*_1.fastq* | parallel --ungroup -j $JOBS call_bwamem {} ${HEADERFILE} ${SUBBAMS_DIR} ${REFERENCE}

    module load picard/1.92
    mem=14g
    inputs=$( ls ${SUBBAMS_DIR}/*.bam | xargs -n1 -I{} echo "INPUT="{} )
    java -Xmx${mem} -Xms${mem} -jar ${PICARDROOT}/MergeSamFiles.jar \
            $inputs OUTPUT=${OUTPUTBAM} VALIDATION_STRINGENCY=LENIENT MERGE_SEQUENCE_DICTIONARIES=true

    module load samtools/0.1.19
    samtools index ${OUTPUTBAM}
    ln -s $( basename $OUTPUTBAM .bam ).bai ${OUTPUTBAM}.bai
fi
