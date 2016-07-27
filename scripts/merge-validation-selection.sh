#!/bin/bash

if [ $# -eq 0 ] || [ -z "$1" ] || [ ! -f "$2" ] 
then
    echo "$0 - merges validation and selection vcfs"
    echo "Usage: $0 validation.vcf selection.vcf"
    exit 
fi

VALIDATION=$1
SELECTION=$2

join -j 1\
    <( grep -v "^#" ${VALIDATION} | sed -e 's/^\([1-9XYh]\)/chr\1/' | sed -e 's/Callers=[^;]*;//' | awk '{id=$1"_"$2"_"$4"_"$5; printf "%s\t%s\n", id, $0}' | sort -k 1,1 | uniq ) \
    <( grep -v "^#" ${SELECTION} | sed -e 's/^\([1-9XYh]\)/chr\1/' | awk '{id=$1"_"$2"_"$4"_"$5; printf "%s\t%s\n", id, $0}' | sort -k 1,1 | uniq ) \
    | awk '{printf "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s;%s\n", $2, $3, $4, $5, $6, $7, $8, $9, $17}'
