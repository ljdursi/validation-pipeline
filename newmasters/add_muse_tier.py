#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf
import sys

def tier_db(muse_vcf_filename):
    reader = vcf.Reader(filename=muse_vcf_filename)
    tier_dict = {}

    for record in reader:
        for alt in record.ALT:
            key = (record.CHROM, record.POS, alt)
            tier_dict[key] = record.FILTER
            print(key, tier_dict[key])

    return tier_dict

def main(): 
    parser = argparse.ArgumentParser(description='Add MUSE tier level to VCF')
    parser.add_argument('input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-m', '--muse', type=str, help='Muse input file', required=True)
    args = parser.parse_args()

    tiers = tier_db(args.muse)
    return 0

if __name__ == "__main__":
    sys.exit(main())
