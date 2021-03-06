#!/bin/bash

module load python/2.7.2

declare -A CASES
CASES[1]=$( cut -f 1,5 metadata/ValidationSamples.csv | tr -d \" | awk '$1 == 1{print $2}' )
CASES[2]=$( cut -f 1,5 metadata/ValidationSamples.csv | tr -d \" | awk '$1 == 2{print $2}' )
CASES[3]=$( cut -f 1,5 metadata/ValidationSamples.csv | tr -d \" | awk '$1 == 3{print $2}' )
CASES[4]=$( cut -f 1,5 metadata/ValidationSamples.csv | tr -d \" | awk '$1 == 4{print $2}' )

#DEFALLCALLDIR=/.mounts/labs/simpsonlab/users/jdursi/pcawg-validation-63/maskrepeats/master_norepeats/
#DEFCOMBINEDIR=combined/germline-realigned
DEFALLCALLDIR=./newmasters/annotated
DEFCOMBINEDIR=./combined/newmasters
DEFOUTDIR=./csvs/newmasters
DEFEXCISE=""

COMBINEDIR=${COMBINEDIR-${DEFCOMBINEDIR}}
ALLCALLDIR=${ALLCALLDIR-${DEFALLCALLDIR}}
OUTDIR=${OUTDIR-${DEFOUTDIR}}
EXCISE=${EXCISE-${DEFEXCISE}}

mkdir -p ${OUTDIR}

for array in 1 2 3 4
do
    for vartype in snv_mnv indel sv
    do
        echo "Array ${array} ${vartype}"
        files=$( for case in ${CASES[$array]}; do ls "${ALLCALLDIR}/${case}.${vartype}.master.vcf" 2> /dev/null ; done )
        if [ ! -z "${files}" ]
        then
            echo " all calls"
            ./scripts/combined_vcf_to_csv.py ${EXCISE} -t ${vartype} $files -o "${OUTDIR}/array${array}_allcalls_${vartype}.csv"
        fi
        files=$( for case in ${CASES[$array]}; do ls "${COMBINEDIR}/${case}.${vartype}.vcf" 2> /dev/null; done )
        if [ ! -z "${files}" ]
        then
            echo " validated calls"
            ./scripts/combined_vcf_to_csv.py ${EXCISE} -t ${vartype} $files -o "${OUTDIR}/array${array}_${vartype}.csv"
        fi
    done
done
