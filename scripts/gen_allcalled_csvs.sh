#!/bin/bash

declare -A CASES
CASES[1]=$( cut -f 1,5 metadata/ValidationSamples.csv | tr -d \" | awk '$1 == 1{print $2}' )
CASES[2]=$( cut -f 1,5 metadata/ValidationSamples.csv | tr -d \" | awk '$1 == 2{print $2}' )
CASES[3]=$( cut -f 1,5 metadata/ValidationSamples.csv | tr -d \" | awk '$1 == 3{print $2}' )
CASES[4]=$( cut -f 1,5 metadata/ValidationSamples.csv | tr -d \" | awk '$1 == 4{print $2}' )

ALLCALLVCFDIR=/.mounts/labs/simpsonlab/users/jdursi/pcawg-validation-63/maskrepeats/master_norepeats/
VALIDATIONVCFDIR=combined/germline-realigned
#VALIDATIONVCFDIR=combined

for array in 1 2 3 4
do
    ALLCALLDIR=${ALLCALLVCFDIR}
    if [ $array == 2 ]
    then
        ALLCALLDIR=/.mounts/labs/simpsonlab/users/jdursi/pcawg-validation-63/master_nohally/
    fi
    for vartype in snv_mnv indel sv
    do
        echo "Array ${array} ${vartype}"
        files=$( for case in ${CASES[$array]}; do ls ${ALLCALLVCFDIR}/${case}.${vartype}.master.vcf 2> /dev/null ; done )
        if [ ! -z "${files}" ]
        then
            echo " all calls"
            ./scripts/combined_vcf_to_csv.py -t ${vartype} $files -o array${array}_allcalls_${vartype}.csv
        fi
        files=$( for case in ${CASES[$array]}; do ls ${VALIDATIONVCFDIR}/${case}.${vartype}.vcf 2> /dev/null; done )
        if [ ! -z "${files}" ]
        then
            echo " validated calls"
            ./scripts/combined_vcf_to_csv.py -t ${vartype} $files -o array${array}_${vartype}.csv
        fi
    done
done
