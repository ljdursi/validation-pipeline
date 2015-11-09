#!/bin/bash -x
module purge
module load bcftools/1.1
module load tabix/0.2.6

if [ $# -eq 0 ]
then
    echo "$0 - generates list of target variants for list of samples given"
    echo "Usage: $0 samplesfile vcfdir outputfile"
    exit 
fi

if [ -z "$1" ]
then
    echo "No samplesfile supplied"
    exit
fi

if [ -z "$2" ]
then
    echo "No Vcf master dir supplied"
    exit
fi

if [ -z "$3" ]
then
    echo "No output file given"
    exit
fi

SAMPLESFILE=$1
VCFDIR=$2
OUT=$3

# get an array from the SAMPLESFILE, and prepend each entry with the directory,
# and end with *.vcf
samples=( $( cut -f 1 $SAMPLESFILE ) )
vcfs=( "${samples[@]/#/${VCFDIR}}" )
snvs=( "${vcfs[@]/%/.snv_mnv.selected.vcf}" )
indels=( "${vcfs[@]/%/.indel.selected.vcf}" )

cat > ${OUT} <<EOF
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
cat ${snvs[@]} ${indels[@]} | grep -v "^#" | sort -k1,1 -k2,2n >> ${OUT}
bgzip ${OUT}
tabix -p vcf ${OUT}.gz
