#!/bin/bash -l
module purge 

if [[ $# -eq 0 ]] || [[ ! -f "$1" ]] || [[ -z "$2" ]] || [[ ! -d "$2" ]] || [[ -z "$3" ]] || [[ ! -d "$3" ]]
then
    >&2 echo "$0 - merges a validation VCF and the union VCF"
    >&2 echo "Usage: $0 FILE.VCF UNIONDIR OUTPUTDIR" 
    >&2 echo "invocation: $0 $1 $2 $3"
    exit 
fi

readonly FILE=$1
readonly WGSDIR=$2
readonly COMBINEDIR=$3

readonly BASE=$( basename $FILE .vcf )
#readonly DONOR=$( echo $BASE | cut -f 1 -d . )
readonly MASTERFILE=${WGSDIR}/${BASE}.master.vcf

if [[ ! -f ${MASTERFILE} ]]
then
    >&2> echo "$0: WGS file ${MASTERFILE} not found"
    exit
fi

./scripts/merge-validation-selection.sh ${FILE} ${MASTERFILE} > ${COMBINEDIR}/${BASE}.vcf
