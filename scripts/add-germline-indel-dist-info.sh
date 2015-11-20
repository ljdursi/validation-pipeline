#!/bin/bash -l
module purge 
module load bedtools/2.21.0

GENOME=/oicr/data/pancanxfer/validation/beds/hg19.genome

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ ! -f "$2" ] 
then
    echo "$0 - annotates somatic calls with distance to nearest germline indel"
    echo "Usage: $0 CALLS.VCF GERMLINE-INDELS.VCF"
    echo "invocation was: $0 $1 $2"
    exit 
fi

INPUT_VCF=$1
GERMLINE_VCF=$2

grep -v "^#" ${GERMLINE_VCF} | sed -e 's/^\([^c]\)/chr\1/' | sort -k1,1 -k2,2n | uniq >> $$.vcf
bedtools closest -header -d -a ${INPUT_VCF} \
    -b $$.vcf \
    | awk '/^[^#]/{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tGermIndelDist=%s;%s\n", $1, $2, $3, $4, $5, $6, $7, $20, $8 } /^#/{print $0}' \
    || sed -e 's/;;/;/'
rm $$.vcf
