#!/bin/bash

module load gcc/4.8.1
module load vcflib
module load bcftools/1.1
module load bedtools/2.24.0
module load tabix

readonly ARRAYFILE=metadata/ValidationSamples.csv
readonly VARIANTDIR=selected-variants/somatic
readonly SLOP=50

function f_vcf_sort {
    grep -v '^#' \
        | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n
}

for array in 1 2 3 4
#for array in 2
do
    for vartype in snv_mnv indel
    do
        arrayvcf=beds/Array_${array}_${vartype}.vcf

        cat beds/header.txt > ${arrayvcf}
        vcfcat ${VARIANTDIR}/../somatic_array_${array}/*.${vartype}.selected.vcf \
            | grep -v "^#" \
            | sed -e 's/^chr//' \
            | f_vcf_sort \
            >> ${arrayvcf}
        bgzip -f ${arrayvcf}
        tabix -p vcf ${arrayvcf}.gz

        arraybed=${arrayvcf%.vcf}.bed
        zcat ${arrayvcf}.gz | vcf2bed.py > ${arraybed}

        inflatedbed=${arraybed%.bed}_inflated.bed
        bedtools slop -i ${arraybed} -b $SLOP -g beds/hg19.genome | bedtools merge -i - > ${inflatedbed}
        bgzip -f ${inflatedbed}
        zcat ${inflatedbed}.gz > ${inflatedbed}
        tabix -p bed ${inflatedbed}.gz

        arraytargets=${arrayvcf%.vcf}.targets
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${arrayvcf}.gz | bgzip -f -c > ${arraytargets}.gz
        tabix -p bed ${arraytargets}.gz

        mkdir -p beds/bysample
        for file in ${VARIANTDIR}/../somatic_array_${array}/*.${vartype}.selected.vcf
        do
            base=beds/bysample/$( basename $file .selected.vcf )

            filebed=${base}.bed
            cat ${file} | vcf2bed.py | sed -e 's/^chr//' > ${filebed}

            fileinflatedbed=${base}_inflated.bed
            bedtools slop -i ${filebed} -b $SLOP -g beds/hg19.genome | bedtools merge -i - | sed -e 's/^chr//' > ${fileinflatedbed}
            bgzip -f ${fileinflatedbed}
            zcat ${fileinflatedbed}.gz > ${fileinflatedbed}
            tabix -p bed ${fileinflatedbed}.gz

            bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' <( cat ${file} | sed -e 's/^chr//' ) | bgzip -f -c > ${base}_targets.gz
            zcat ${base}_targets.gz > ${base}_targets
            tabix -p vcf ${base}_targets.gz
        done
    done
done
