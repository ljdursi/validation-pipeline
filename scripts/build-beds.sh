#!/bin/bash

module load vcflib
module load bcftools/1.1
module load bedtools/2.24.0
module load tabix

ARRAYFILE=metadata/ValidationSamples.csv
VARIANTDIR=selected-variants/somatic
SLOP=50

#for array in 1 2 3 4
for array in 2
do
    for vartype in snv_mnv indel
    do
        arrayvcf=beds/Array_${array}_${vartype}.vcf

        vcfcat ${VARIANTDIR}/../somatic_array_${array}/*.${vartype}.selected.vcf | vcfsort | sed -e 's/^chr//' > ${arrayvcf}
        bgzip ${arrayvcf}
        tabix -p vcf ${arrayvcf}.gz

        arraybed=${arrayvcf%.vcf}.bed
        zcat ${arrayvcf}.gz | vcf2bed.py > ${arraybed}

        inflatedbed=${arraybed%.bed}_inflated.bed
        bedtools slop -i ${arraybed} -b $SLOP -g beds/hg19.genome | bedtools merge -i - > ${inflatedbed}
        bgzip ${inflatedbed}
        tabix -p bed ${inflatedbed}

        arraytargets=${arrayvcf%.vcf}.targets
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${arrayvcf}.gz > ${arraytargets}
    done
done
