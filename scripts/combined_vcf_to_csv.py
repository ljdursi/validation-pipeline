#!/usr/bin/env python
"""
Take combined VCF files and translate them to CSV for analysis
"""
from __future__ import print_function
import argparse
import sys
import collections
from os import listdir
from os.path import isfile, isdir, basename

def build_filelist(filenames, match):
    """ If given a directory, descend into directory and generate 
        list of files that match "match" (but do not descend into subdirectories)."""
    file_list = []
    for filename in filenames:
        if isfile(filename) and filename.find(match) != -1:
            file_list.append(filename)
        elif isdir(filename):
            for f in listdir(filename):
                pathname = os.path.join(filename, f)
                if isfile(pathname) and filename.find(match) != -1:
                    file_list.append(pathname)
    return file_list


def parse_combined():
    snv_callers=['adiscan', 'broad_mutect', 'dkfz', 'lohcomplete', 'mda_hgsc_gatk_muse', 
                 'oicr_bl', 'oicr_sga', 'sanger', 'smufin', 'wustl']
    indel_callers=['broad_mutect', 'crg_clindel', 'dkfz', 'novobreak', 'oicr_sga', 
                   'sanger', 'smufin', 'wustl']
    sv_callers=['broad_merged', 'destruct', 'embl_delly', 'novobreak', 'sanger', 'smufin']

    parser = argparse.ArgumentParser( description='Set genotypes based on DP4 scores')
    parser.add_argument('inputs', nargs="+")
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
    files = build_filelist(args.inputs, args.variant_type)

    for f in files:
        sample_id = basename(f).split(".")[0]
        #for line in open(args.directory + "/" + f):
        for line in open(f, 'r'):
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
            
            if not validation_tumour_var_depth is None and not validation_tumour_depth is None:
                itemdict['val_tvaf'] = str(1.0*validation_tumour_var_depth/(validation_tumour_depth+0.00001))
            if not validation_normal_var_depth is None and not validation_normal_depth is None:
                itemdict['val_nvaf'] = str(1.0*validation_normal_var_depth/(validation_normal_depth+0.00001))
            print(','.join([itemdict[h] for h in headerfields]), file=args.output)

if __name__ == "__main__":
    parse_combined()

