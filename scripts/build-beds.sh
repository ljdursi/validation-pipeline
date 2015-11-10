#!/bin/bash

module load gcc/4.8.1
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

        vcfcat ${VARIANTDIR}/../somatic_array_${array}/*.${vartype}.selected.vcf | vcfsort | sed -e 's/^chr//' | grep -v 'Callers=smufin;' > ${arrayvcf}
        bgzip ${arrayvcf}
        tabix -p vcf ${arrayvcf}.gz

        arraybed=${arrayvcf%.vcf}.bed
        zcat ${arrayvcf}.gz | vcf2bed.py > ${arraybed}

        inflatedbed=${arraybed%.bed}_inflated.bed
        bedtools slop -i ${arraybed} -b $SLOP -g beds/hg19.genome | bedtools merge -i - > ${inflatedbed}
        bgzip ${inflatedbed}
        tabix -p bed ${inflatedbed}.gz

        arraytargets=${arrayvcf%.vcf}.targets
        bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' ${arrayvcf}.gz | bgzip -c > ${arraytargets}.gz
        tabix -p vcf ${arraytargets}.gz

        mkdir -p beds/bysample
        for file in ${VARIANTDIR}/../somatic_array_${array}/*.${vartype}.selected.vcf
        do
            base=beds/bysample/$( basename $file .selected.vcf )

            filebed=${base}.bed
            cat ${file} | vcf2bed.py > ${filebed}

            fileinflatedbed=${base}_inflated.bed
            bedtools slop -i ${filebed} -b $SLOP -g beds/hg19.genome | bedtools merge -i - | sed -e 's/^chr//' > ${fileinflatedbed}
            bgzip ${fileinflatedbed}
            tabix -p bed ${fileinflatedbed}.gz

            bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' <( cat ${file} | grep -v 'Callers=smufin;' | sed -e 's/^chr//' ) | bgzip -c > ${base}_targets.gz
            tabix -p vcf ${base}_targets.gz
        done
    done
done
