#!/bin/bash

DEFREADCOUNTSDIR=readcounts/germline-realigned
DEFBAMDIR=germ-indel-realigned

## can be overridden at commandline, eg BAMDIR=remappedbams READCOUNTSDIR=readcounts/not-realigned ./run-readcounts.sh
READCOUNTSDIR=${READCOUNTSDIR:-${DEFREADCOUNTSDIR}}
BAMDIR=${BAMDIR:-${DEFBAMDIR}}

echo "bams: ${BAMDIR}"
echo "readcounts: ${READCOUNTSDIR}"

mkdir -p $READCOUNTSDIR

for file in ${BAMDIR}/*.bam
do
    base=$( basename $file .bam )
    donor=$( grep $base metadata/sample-directory.txt | awk '{print $1}' )

    for vartype in snv_mnv 
    do
        if [ ! -f ${READCOUNTSDIR}/${base}.${vartype}.rc ]
        then
            qsub -cwd -e logs -o logs -l h_vmem=4g ./scripts/one-readcount.sh $file beds/bysample/${donor}.${vartype}.bed ${READCOUNTSDIR}/${base}.${vartype}.rc
        fi
    done
done 
