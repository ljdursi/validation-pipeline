#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

DEFREF=/oicr/data/pancanxfer/validation/reference/bwa-0.6.2/genome.fa

if [ $# -eq 0 ] || [ -z "$1" ] 
then
    echo "$0 - generates IGV batch file for visualizing variants on a VCF"
    echo "Usage: $0 variants.vcf"
    exit 
fi

VCF=$1

donor=$( basename $VCF .vcf )
donor=$( echo $donor | cut -f 1 -d \. )

normal=$( grep $donor metadata/sample-directory.txt | awk '{print $2}')
tumour=$( grep $donor metadata/sample-directory.txt | awk '{print $3}')

if [ -z "$normal" ] || [ -z "$tumour" ]
then
    exit
fi

cat << EOF
new
genome hg19
load ../remappedbams/${normal}
load ../remappedbams/${tumour}
snapshotDirectory ../plots/${donor}
EOF

grep -v "^#" $VCF | while read -r line || [[ -n "$line" ]]; do
    items=(${line//	/ })

    chrom=${items[0]}
    pos=${items[1]}
    ref=${items[3]}
    alt=${items[4]}

    echo "goto ${chrom}:${pos}"
    echo "collapse"
    echo "snapshot ${chrom}_${pos}_${ref}_${alt}.png"
done 

echo "exit"
