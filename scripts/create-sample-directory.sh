#!/bin/bash 
#
# creates a sample directory - a table of (case, control.bam, tumour.bam) -
# from the data given in metadata and the cghub query
#

QUERYFILENAME=${1-"ingest-bams/xfer-cghub/query4"}
SAMPLECSV=${2-"metadata/sample_metadata.csv"}

SAMPLETOBAM="sample-to-bam.txt"

join -j 1 \
    <( egrep  "(filename)|( sample_id)" ${QUERYFILENAME} \
        | awk '{ printf "%s",$3; if (n == 2) { print ""; n=0 } else { printf "\t"; n++ }}' \
        | awk '{print $3, $1}' | sort ) \
    <( cut -f 4,11,16 -d , $SAMPLECSV | tr , " " \
        | awk '{printf "%s %s CONTROL\n%s %s TUMOUR\n", $1, $3, $2, $3}' | sort ) \
    > ${SAMPLETOBAM}

join -j 1 \
    <( cat ${SAMPLETOBAM} | grep CONTROL | awk '{print $3, $2}' | sort ) \
    <( cat ${SAMPLETOBAM} | grep TUMOUR | awk '{print $3, $2}' | sort ) 

rm ${SAMPLETOBAM}
