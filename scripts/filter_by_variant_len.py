#!/usr/bin/env python
"""
VCF filter; annotate variants with their variant length (indels and some SVs), 
    and filter variants which have a absolute-value variant length below some value
"""
from __future__ import print_function
import argparse
import vcf
import sys

def variant_length(record):
    """
    Returns:    Length of variants for an indel or sv
                Returns max abs value of variant for multi-allelic calls
                Returns None for no length - eg, translocation
    """
    ref = record.REF
    lens = []
    for alt in record.ALT:
        if isinstance(alt, vcf.model._Substitution):
            lens.append(len(alt)-len(ref))
        elif isinstance(alt, vcf.model._Breakend):
            if alt.chr == record.CHROM:
                lens.append(record.POS - alt.pos)
            
    if len(lens) == 0:
        return None
    return sorted(lens, key=lambda x:-abs(x))[0]

def svlen_filter():
    parser = argparse.ArgumentParser( description='Set genotypes based on DP4 scores')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-m', '--min_abs_len', type=int, default=-1, help='minimum variant length; -1 to not filter')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(args.input)
    vcf_writer = vcf.Writer(args.output, vcf_reader)

    for record in vcf_reader:
        varlen = variant_length(record)
        if varlen is not None:
            if args.min_abs_len > -1 and abs(varlen) < args.min_abs_len:
                continue
            record.INFO['varlen'] = [varlen]
        vcf_writer.write_record(record)

if __name__ == "__main__":
    svlen_filter()
