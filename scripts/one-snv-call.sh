#!/bin/bash -l
module purge 
module use /.mounts/labs/simpsonlab/modules/
module load python-packages/2

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ ! -f "$1" ] || [ ! -f "$2" ] || [ ! -f "$3" ]
then
    echo "$0 - performs call validation based on readcounts files"
    echo "Usage: $0 INPUT.VCF NORMAL_READCOUNTS.RC TUMOUR_READCOUNTS.RC OUTPUT.VCF"
    echo "$0 $1 $2 $3 $4"
    exit 
fi

INPUT_VCF=$1
NORMAL_RC=$2
TUMOUR_RC=$3
OUTPUT_VCF=$4

if [ ! -f "$OUTPUT_VCF" ]
then
    cat > ${OUTPUT_VCF} <<EOF
##fileformat=VCFv4.1
##INFO=<ID=Callers,Number=.,Type=String,Description="Callers that made this call">
##INFO=<ID=NormalEvidenceReads,Number=2,Type=Integer,Description="Number of reads in normal sample supporting the alt allele (forward,backward)">
##INFO=<ID=TumourEvidenceReads,Number=2,Type=Integer,Description="Number of reads in tumour sample supporting the alt allele (forward,backward)">
##INFO=<ID=TumourReads,Number=1,Type=Integer,Description="Total number of reads in tumour sample">
##INFO=<ID=NormalReads,Number=1,Type=Integer,Description="Total number of reads in normal sample">
##INFO=<ID=Validation_status,Number=.,Type=String,Description="Validation status, as info field">
##FILTER=<ID=LOWDEPTH,Description="Insufficient depth to validate the call">
##FILTER=<ID=NOTSEEN,Description="Variant not seen in Tumour">
##FILTER=<ID=STRANDBIAS,Description="Too much strand bias in Tumour sample to validate the call">
##FILTER=<ID=NORMALEVIDENCE,Description="Evidence for presence in normal sample">
##FILTER=<ID=GERMLINE,Description="Germline het or hom">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
    ./scripts/snv_readcounts.py ${INPUT_VCF} ${NORMAL_RC} ${TUMOUR_RC} | \
        ./scripts/snv_indel_call.py --error 0.01 --callthreshold 0.02 --strandbias 0.001 --germlineprob 0.01 --mindepth 25 \
        | sed -e 's/;;/;/' \
        | grep -v "^#" \
        | sort -k1,1 -k2,2n  >> ${OUTPUT_VCF}
fi
