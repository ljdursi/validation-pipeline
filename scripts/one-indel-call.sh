#!/bin/bash -l
module purge 
module use /.mounts/labs/simpsonlab/modules/
module load python-packages/2

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ ! -f "$1" ] 
then
    echo "$0 - performs call validation based on readcounts files"
    echo "Usage: $0 VALIDATION_VCF_WITH_COUNTS.VCF OUTPUT.VCF"
    echo "$0 $1 $2"
    exit 
fi

INPUT_VCF=$1
OUTPUT_VCF=$2

if [ ! -f "$OUTPUT_VCF" ]
then
    cat > ${OUTPUT_VCF} <<EOF
##fileformat=VCFv4.1
##INFO=<ID=Callers,Number=.,Type=String,Description="Callers that made this call">
##INFO=<ID=NormalEvidenceReads,Number=1,Type=Integer,Description="Number of reads in normal sample supporting the alt allele">
##INFO=<ID=TumourEvidenceReads,Number=1,Type=Integer,Description="Number of reads in tumour sample supporting the alt allele">
##INFO=<ID=TumourReads,Number=1,Type=Integer,Description="Total number of reads in tumour sample">
##INFO=<ID=NormalReads,Number=1,Type=Integer,Description="Total number of reads in normal sample">
##INFO=<ID=Validation_status,Number=.,Type=String,Description="Validation status, as info field">
##INFO=<ID=RepeatRefCount,Number=1,Type=Integer,Description="Reference repeat unit count at this location">
##INFO=<ID=RepeatUnit,Number=1,Type=String,Description="Reference repeat unit">
##INFO=<ID=3pContext,Number=1,Type=String,Description="Genomic context at 3-end">
##INFO=<ID=5pContext,Number=1,Type=String,Description="Genomic context at 3-end">
##FILTER=<ID=LOWSUPPORT,Description="Insufficient depth to validate the call">
##FILTER=<ID=NOTSEEN,Description="Variant not seen in Tumour">
##FILTER=<ID=STRANDBIAS,Description="Too much strand bias in Tumour sample to validate the call">
##FILTER=<ID=NORMALEVIDENCE,Description="Evidence for presence in normal sample">
##FILTER=<ID=GERMLINE,Description="Germline het or hom">
EOF
    ./scripts/snv_indel_call.py --error 0.01 --callthreshold 0.02 --mindepth 25 --strandbias -1 --germlineprob 0.05 \
        -i ${INPUT_VCF} \
        | sed -e 's/;;/;/' \
        | grep -v "^##" \
        | sort -k1,1 -k2,2n | sed -e 's/^\([1-9XYh]\)/chr\1/' >> ${OUTPUT_VCF}
fi
