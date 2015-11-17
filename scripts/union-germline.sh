#!/bin/bash

module purge 
module unuse /oicr/local/boutroslab/Modules/modulefiles
module load gcc/4.8.1
module load tabix
module use /.mounts/labs/simpsonlab/modules/
module load vcflib

DIR=germline-calls
cases=$( ls ${DIR}/*.vcf.gz | cut -f 1 -d . | sed -e "s#${DIR}/##" | sort | uniq  )
for pair in ${cases}
do
    vcf=${DIR}/union/${pair}.vcf
    echo ${pair}
    echo "##fileformat=VCFv4.1" > $vcf
    echo '##INFO=<ID=GERMLINE,Number=0,Type=Flag,Description="Indicates if record is a germline mutation">' >> ${vcf}
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">' >> ${vcf}
    echo '##SAMPLE=<ID=CONTROL,SampleName=control_NA,Individual=NA,Description="Control">' >> ${vcf}
    echo '##SAMPLE=<ID=TUMOR,SampleName=tumor_NA,Individual=NA,Description="Tumor">' >> ${vcf}
    echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CONTROL	TUMOR" >> ${vcf}
    zcat ${DIR}/${pair}*vcf.gz \
        | awk '$1 !~ /^#/ && (length($4) != length($5)) && ($7 ~ /\./ || $7 ~ /PASS/) { chr=$1; if (chr ~/[0-9XY]/) {chr="chr"chr}; printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tGERMLINE\tGL\t0/1\t0/1\n", $1, $2, $3, $4, $5, $6, $7}'  > tmp.$$
    grep "^#" tmp.$$ >> ${vcf}
    grep "^[1-9]" tmp.$$ | sort -k1,1n -k2,2n | uniq >> ${vcf}
    grep "^[XY]" tmp.$$ | sort -k1,1 -k2,2n | uniq >> ${vcf}
    grep "^MT" tmp.$$ | sort -k1,1 -k2,2n | uniq >> ${vcf}
    grep "^hs37d5" tmp.$$ | sort -k1,1 -k2,2n | uniq >> ${vcf}
    rm tmp.$$
    bgzip ${vcf}
    tabix -p vcf ${vcf}.gz
    zcat ${vcf}.gz > ${vcf}
done
