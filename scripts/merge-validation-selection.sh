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
    <( grep -v "^#" ${VALIDATION} | awk '$8 !~ /^DP=[0-9];/ && $10 !~ /^\.\/\./ && $11 !~ /^\.\/\./ {id=$1"_"$2; if ($5 == "." || $11 ~ /0\/0:/) { filter="NOTSEEN";} else if ($10 !~ /0\/0:/) { filter = "GERMLINE" } else filter="PASS"; printf "%s\t%s\t%s\n", id, filter, $0}' | sort -k 1,1 ) \
    <( grep -v "^#" ${SELECTION} | sed -e 's/^chr//' | awk '{id=$1"_"$2; printf "%s\t%s\n", id, $0}' | sort -k 1,1 ) \
    | awk '{printf "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s;%s\t%s\t%s\t%s\n", $3, $4, $5, $6, $18, $8, $2, $10, $21, $11, $12, $13}'
