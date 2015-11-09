#!/bin/bash 
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/
module load bedtools/2.24.0

if [ $# -eq 0 ]
then
    echo "$0 - filters vcf files by bed"
    echo "Usage: $0 vcffile bedfile"
    exit 
fi

if [ -z "$1" ]
then
    echo "No vcffile supplied"
    exit
fi

if [ -z "$2" ]
then
    echo "No bedfile supplied"
    exit
fi

VCF=$1
BED=$2


zgrep "^#" $VCF 
bedtools intersect -a $VCF -b $BED -wa 
