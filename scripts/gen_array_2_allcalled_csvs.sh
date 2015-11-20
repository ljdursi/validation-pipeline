#!/bin/bash

# this has to be done manually with a different set of files for now,
# because everyone else used a different set of files, which had been
# masked for repeats

CASES="1eb36468-c560-4041-b2d9-89c464425a6e\
       1f1a065b-1458-4846-99b3-7370bbf7b367\
       211554b9-3e23-4b17-a0b6-e773495457a9\
       7867d1aa-c1b4-41fc-9ff1-c45467ce0ad9\
       7d23a9b6-f59d-4e93-9987-48e740a60159\
       7e39feb6-c04e-47a3-b5a6-b7b59a7fc013\
       a6e8dd23-c8a5-445a-ae4b-b9f92ed6a73e\
       b4c03bfa-fc41-4568-9006-0a2b2ba56ddf\
       ce8fd3aa-31a0-441f-ae4a-8e3114942f16\
       d080db6b-583b-46fe-9e2b-b70069ebe960\
       f740d082-ae6f-4c26-82c0-43ef49680d1d"

VCFDIR=/.mounts/labs/simpsonlab/users/jdursi/pcawg-validation-63/master_nohally/
for vartype in snv_mnv indel sv
do
    files=$( for case in $CASES; do ls ${VCFDIR}/${case}.${vartype}.master.vcf; done )
    ./scripts/combined_vcf_to_csv.py -t ${vartype} $files -o array2_allcalls_${vartype}.csv
done
