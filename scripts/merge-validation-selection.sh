#!/bin/bash

if [ $# -eq 0 ] || [ ! -f "$1" ] || [ ! -f "$2" ] 
then
    echo "$0 - merges validation and selection vcfs"
    echo "Usage: $0 validation.vcf selection.vcf"
    exit 
fi

VALIDATION=$1
SELECTION=$2

join -j 1\
    <( grep -v "^#" ${VALIDATION} | awk '{id=$1"_"$2"_"$3"_"$4; printf "%s\t%s\n", id, $0}' | sort -k 1,1 ) \
    <( grep -v "^#" ${SELECTION} |  awk '{id=$1"_"$2"_"$3"_"$4; printf "%s\t%s\n", id, $0}' | sort -k 1,1 ) \
    | awk '{printf "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s;%s\n", $2, $3, $4, $5, $6, $7, $8, $9, $17}'
