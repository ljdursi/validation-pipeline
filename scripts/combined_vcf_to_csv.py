#!/usr/bin/env python
"""
VCF filter; adjust genotypes based on DP4 (4-depth) scores
"""
from __future__ import print_function
import argparse
import vcf
import sys
import numpy
import scipy.stats
import collections
from os import listdir
from os.path import isfile, join, splitext

def parse_combined():
    snv_callers=['adiscan', 'broad_mutect', 'dkfz', 'lohcomplete', 'mda_hgsc_gatk_muse', 
                 'oicr_bl', 'oicr_sga', 'sanger', 'smufin', 'wustl']
    indel_callers=['broad_mutect', 'crg_clindel', 'dkfz', 'novobreak', 'oicr_sga', 
                   'sanger', 'smufin', 'wustl']
    sv_callers=['broad_merged', 'destruct', 'embl_delly', 'novobreak', 'sanger', 'smufin']

    parser = argparse.ArgumentParser( description='Set genotypes based on DP4 scores')
    parser.add_argument('directory')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-t', '--variant_type', type=str, choices=["snv_mnv","indel","sv"], default="snv_mnv")
    args = parser.parse_args()

    if args.variant_type == "indel":
        caller_names = indel_callers
    elif args.variant_type == "snv_mnv":
        caller_names = snv_callers
    else:
        caller_names = sv_callers

    headerfields = ['chrom', 'pos', 'sample', 'status', 'ref', 'alt', 'val_tvaf', 'val_nvaf', 'wgs_tvaf', 'wgs_nvaf', 'repeat_count', 'varlen', 'indel_dist' ]
    headerfields = headerfields + caller_names
    print(','.join(headerfields), file=args.output)

    # Get all of the files in the directory for this variant type
    files = [ f for f in listdir(args.directory) if isfile(args.directory + "/" + f) and f.find(args.variant_type) != -1 ]

    for f in files:
        sample_id = f.split(".")[0]
        for line in open(args.directory + "/" + f):
            if line[0] == '#':
                continue
            itemdict = collections.defaultdict(str)
            items = line.split()
            callers_found = False
            validation_tumour_depth, validation_tumour_var_depth = None, None
            validation_normal_depth, validation_normal_var_depth = None, None
            itemdict['chrom'], itemdict['pos'] = items[0], items[1]
            itemdict['status'], itemdict['ref'], itemdict['alt'] = items[6], items[3], items[4]
            itemdict['sample'] = sample_id
            itemdict['indel_dist'] = "NA"
            if 'GERMLINE' in itemdict['status']:
                itemdict['status'] = 'GERMLINE'
            if 'NOTSEEN' in itemdict['status']:
                itemdict['status'] = 'NOTSEEN'
            for field in items[7].split(';'):
                splitfield = field.split('=')
                if len(splitfield) < 2:
                    continue
                key, val = splitfield[0], splitfield[1]
                if key == 'Callers' and not callers_found:
                    callers_found = True
                    for caller in caller_names:
                        if caller in field:
                            itemdict[caller] = '1'
                        else:
                            itemdict[caller] = '0'
                if key == 'RepeatRefCount':
                    itemdict['repeat_count'] = val
                if key == 'RepeatRefCount':
                    itemdict['repeat_count'] = val
                if key == 'GermIndelDist':
                    itemdict['indel_dist'] = val
                if key == 'NormalVAF':
                    itemdict['wgs_nvaf'] = val
                if key == 'TumorVAF':
                    itemdict['wgs_tvaf'] = val
                if key == 'NormalEvidenceReads':
                    validation_normal_var_depth = sum([int(x) for x in val.split(',')])
                if key == 'NormalReads':
                    validation_normal_depth = int(val)
                if key == 'TumourEvidenceReads':
                    validation_tumour_var_depth = sum([int(x) for x in val.split(',')])
                if key == 'TumourReads':
                    validation_tumour_depth = int(val)
                    
            itemdict['val_tvaf'] = str(1.0*validation_tumour_var_depth/(validation_tumour_depth+0.00001))
            itemdict['val_nvaf'] = str(1.0*validation_normal_var_depth/(validation_normal_depth+0.00001))
            print(','.join([itemdict[h] for h in headerfields]), file=args.output)

if __name__ == "__main__":
    parse_combined()

