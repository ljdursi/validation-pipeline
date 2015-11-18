#!/bin/bash -l
module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module use /.mounts/labs/simpsonlab/modules/

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
wgsnormal=$( grep $donor metadata/sample_metadata.csv | cut -f 10 -d , )
wgstumour=$( grep $donor metadata/sample_metadata.csv | cut -f 17 -d , )

if [ -z "$normal" ] || [ -z "$tumour" ] 
then
    exit
fi

cat << EOF
new
genome hg19
load ../remappedbams/${normal}
load ../remappedbams/${tumour}
load ../wgs-bams/${wgsnormal}
load ../wgs-bams/${wgstumour}
snapshotDirectory ../plots/${donor}
EOF

grep -v "^#" $VCF | while read -r line || [[ -n "$line" ]]; do
    items=(${line//	/ })

    chrom=${items[0]}
    pos=${items[1]}
    ref=${items[3]}
    alt=${items[4]}

    if (( ${#ref} == ${#alt} ))
    then
        vartype="${ref}_${alt}"
        max=${ref}
    elif (( ${#ref} > ${#alt} )) 
    then
        vartype="del"
        max=$ref 
    else
        vartype="ins"
        max=$alt 
    fi

    startp=$(( $pos - 10 ))
    endp=$(( $pos + ${#max} + 10 ))

    echo "goto ${chrom}:${startp}-${endp}"
    echo "collapse"
    echo "snapshot ${chrom}_${pos}_${vartype}.png"
done 

echo "exit"
