#!/bin/bash

module load bedtools/2.24.0

function pairstobed {
    SAMPLE=$1
    BEDFILE=$2

    grep ${SAMPLE} morris-pairs.txt | awk '{print $3, $4}' | \
        tr '_	' ' ' | awk '{printf "%s\t%d\t%d\n", $1, $2, $4}' | \
        bedtools sort | bedtools slop -b 20 -g ../beds/hg19.genome | bedtools merge -i - > $BEDFILE
}

pairstobed b02b4bba-6e66-44fb-a48f-38c309aaaac5 ce8fd3aa-31a0-441f-ae4a-8e3114942f16.snv_mnv.evolution.bed
pairstobed f848b66f-bd9e-4fba-afd4-eb58848d1ef4 b92ab845-7c4c-4498-88c2-75c2cb770b62.snv_mnv.evolution.bed
