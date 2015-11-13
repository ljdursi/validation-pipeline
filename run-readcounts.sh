#!/bin/bash

READCOUNTDIR=readcounts
mkdir -p $READCOUNTDIR

for file in remappedbams/*.bam
do
    base=$( basename $file .bam )
    donor=$( grep $base metadata/sample-directory.txt | awk '{print $1}' )

    for vartype in snv_mnv indel
    do
        qsub -cwd -e logs -o logs -l h_vmem=4g ./scripts/one-readcount.sh $file beds/bysample/${donor}.${vartype}.bed ${READCOUNTDIR}/${base}.${vartype}.rc
    done
done 
