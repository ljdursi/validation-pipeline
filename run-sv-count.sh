#!/bin/bash

DEFBAMDIR=germ-indel-realigned
DEFREADCOUNTDIR=readcounts/germline-realigned

BAMDIR=${BAMDIR-$DEFBAMDIR}
READCOUNTDIR=${READCOUNTDIR-$DEFREADCOUNTDIR}

while read -r line || [[ -n "$line" ]]; do
    items=(${line//	/ })

    normal=${items[1]}
    tumour=${items[2]}
    donor=${items[0]}

    if [ -f ${BAMDIR}/${normal} ] && [ -f ${BAMDIR}/${tumour} ] 
    then
        for vartype in sv
        do
            if [ -f selected-variants/somatic/${donor}.${vartype}.selected.vcf ]
            then
                qsub -cwd -e logs -o logs -l h_vmem=2g -pe smp 2 \
                    ./scripts/one-sv-count.sh \
                    selected-variants/somatic/${donor}.${vartype}.selected.vcf \
                    ${BAMDIR}/${normal} ${BAMDIR}/${tumour} \
                    ${READCOUNTDIR}/${donor}.${vartype}.vcf
            fi
        done
    fi
done < metadata/sample-directory.txt
