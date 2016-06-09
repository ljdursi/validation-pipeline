#!/bin/bash 
module purge 
module load bcftools/1.1
module load tabix/0.2.6

REFERENCE=/.mounts/labs/simpsonlab/data/references/hs37d5.fa 

for vcffile in "$@"
do
    if [ ! -f ${vcffile}.tbi ]; then
        bgzip ${vcffile}
        tabix -p vcf ${vcffile}.gz
        vcffile=${vcffile}.gz
    fi
    out="${vcffile%.*.*}".norm.vcf 
    bcftools norm --check-ref x -m- -f $REFERENCE -o $out $vcffile 
    bgzip $out
    tabix -p vcf ${out}.gz
done
