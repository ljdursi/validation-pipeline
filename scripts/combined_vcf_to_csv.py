#!/usr/bin/env python
"""
Take combined VCF files and translate them to CSV for analysis
"""
from __future__ import print_function
import argparse
import sys
import collections
from os import listdir
import os.path
from os.path import isfile, isdir, basename

def build_filelist(filenames, match):
    """If given a directory, descend into directory and generate
       list of files that match "match" (but do not descend into
       subdirectories)."""
    file_list = []
    for filename in filenames:
        if isfile(filename) and filename.find(match) != -1:
            file_list.append(filename)
        elif isdir(filename):
            for onefile in listdir(filename):
                pathname = os.path.join(filename, onefile)
                if isfile(pathname) and filename.find(match) != -1:
                    file_list.append(pathname)
    return file_list


def parse_combined(files, output, callers, excise=False):
    """
    Parse a list of combined VCFs into a single CSV
    """
    headerfields = ['chrom', 'pos', 'sample', 'status', 'ref', 'alt',
                    'val_tvaf', 'val_nvaf', 'val_tdepth', 'val_ndepth']

    # some columns require special handling: most just involve
    # copying the field, perhaps applying some transformation.  This
    # diectionary allows us to automate the simple cases:

    def sum_ints(comma_delimited_string):
        """ Parse a comma-delimited list of integers and return the sum"""
        return sum([int(x) for x in comma_delimited_string.split(',')])

    def binarize(number):
        """ Return 0 if the string is integer 0, or 1 otherwise """
        return 0 if int(number) == 0 else 1

    def int_or_NA(maybe_int):
        """ Return an int if it can be converted as such, or "NA" """
        if maybe_int == "":
            return "NA"
        return int(maybe_int)

    def float_or_NA(maybe_float):
        """ Return an float if it can be converted as such, or "NA" """
        try:
            return float(maybe_float)
        except Exception as e:
            return "NA"


    fielddict = {'varlen':('varlen', int_or_NA),
                 'repeat_masker':('repeat_masker', binarize),
                 'cosmic':('cosmic', binarize),
                 'dbsnp':('dbsnp', binarize),
                 'gencode_prioritized':('gencode', str),
                 'wgs_NormalReads':('wgs_ndepth', int_or_NA),
                 'wgs_TumorReads':('wgs_tdepth', int_or_NA),
                 'wgs_TumorAvgVarBaseQ':('wgs_tvar_avgbaseq', float_or_NA),
                 'wgs_TumorEvidenceReads':('wgs_tvardepth', sum_ints),
                 'wgs_TumorVarDepth':('wgs_tvardepth', int_or_NA),
                 'wgs_NormalEvidenceReads':('wgs_nvardepth', sum_ints),
                 'wgs_NormalVarDepth':('wgs_nvardepth', int_or_NA),
                 'wgs_TumorAvgVarPosn':('wgs_tvar_avgbaseposn', float_or_NA),
                 'GermIndelDist':('indel_dist', int_or_NA),
                 'wgs_NormalVAF':('wgs_nvaf', float_or_NA),
                 'wgs_TumorVAF':('wgs_tvaf', float_or_NA),
                 'wgs_RepeatRefCount':('repeat_count', int_or_NA),
                 'muse_feature':('muse_feature', int_or_NA)
                 }

    headerfields = headerfields + [x for x, _ in fielddict.values()] +  callers
    print(','.join(headerfields), file=output)

    # Get all of the files in the directory for this variant type
    for combined_file in files:
        sample_id = basename(combined_file).split(".")[0]
        for line in open(combined_file, 'r'):
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

            # Split the info fields
            for field in items[7].split(';'):
                splitfield = field.split('=')
                if len(splitfield) != 2:
                    continue
                key, val = splitfield[0], splitfield[1]

                # handle the simple cases first
                if key in fielddict:
                    label, fun = fielddict[key]
                    try:
                        itemdict[label] = fun(val)
                    except Exception as e:
                        print("Failed converting "+str(key)+" = "+str(val), file=sys.stderr)
                        print("line = <"+line+"> ", file=sys.stderr)
                        print(e, file=sys.stderr)

                # Deal with the callers: only do this once
                if key == 'Callers' and not callers_found:
                    callers_found = True
                    for caller in callers:
                        if caller in field:
                            itemdict[caller] = '1'
                        else:
                            itemdict[caller] = '0'

                if key == 'NormalReads':
                    validation_normal_depth = int(val)
                    itemdict['val_ndepth'] = validation_normal_depth
                if key == 'NormalEvidenceReads':
                    validation_normal_var_depth = sum_ints(val)
                if key == 'TumourReads':
                    validation_tumour_depth = int(val)
                    itemdict['val_tdepth'] = validation_tumour_depth
                if key == 'TumourEvidenceReads':
                    validation_tumour_var_depth = sum_ints(val)

            if not validation_tumour_var_depth is None and not validation_tumour_depth is None:
                itemdict['val_tvaf'] = str(1.0*validation_tumour_var_depth/(validation_tumour_depth+0.00001))
            if not validation_normal_var_depth is None and not validation_normal_depth is None:
                itemdict['val_nvaf'] = str(1.0*validation_normal_var_depth/(validation_normal_depth+0.00001))

            if excise:
                itemdict['chrom'] = '1'
                itemdict['pos'] = '1'
            print(','.join([str(itemdict[h]) for h in headerfields]), file=output)

def main():
    """ Main program: parse arguments, call parse_combined """
    snv_callers = ['adiscan', 'broad_mutect', 'dkfz', 'lohcomplete',
                   'mda_hgsc_gatk_muse', 'oicr_bl', 'oicr_sga', 'sanger',
                   'smufin', 'wustl']
    indel_callers = ['broad_mutect', 'crg_clindel', 'dkfz', 'novobreak',
                     'oicr_sga', 'sanger', 'smufin', 'wustl']
    sv_callers = ['broad_merged', 'destruct', 'embl_delly', 'novobreak',
                  'oicr_bl', 'sanger', 'smufin', 'wustl']

    parser = argparse.ArgumentParser()
    parser.add_argument('inputs', nargs="+")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-x', '--excise', action="store_true")
    parser.add_argument('-t', '--variant_type', type=str,
                        choices=["snv_mnv", "indel", "sv"], default="snv_mnv")
    args = parser.parse_args()

    if args.variant_type == "indel":
        caller_names = indel_callers
    elif args.variant_type == "snv_mnv":
        caller_names = snv_callers
    else:
        caller_names = sv_callers

    files = build_filelist(args.inputs, args.variant_type)
    parse_combined(files, args.output, caller_names, args.excise)
    return 0

if __name__ == "__main__":
    sys.exit(main())

