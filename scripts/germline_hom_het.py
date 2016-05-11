#!/usr/bin/env python
"""
Routines for reading in paired bam-readcount files, 
and comparing to a somatic VCF
"""
from __future__ import print_function
import argparse
import vcf
import sys
import numpy
import scipy.stats

def main():
    """
    Driver program - Read in a VCF file and normal/tumour read counts
    for each base at each position, and output read counts at each call
    """
    parser = argparse.ArgumentParser(description='Search validation data for germline homs/hets')
    parser.add_argument('vcffile', nargs='+', help='Vcf file(s) to check')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'),
                        default=sys.stdout, help='Output VCF (default: stdout)')
    parser.add_argument('-e', '--errorrate', type=float, default=0.02,
                        help='Error rate')
    parser.add_argument('-a', '--alpha', type=float, default=0.05,
                        help='prob threshold for calling hets/homs')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(filename=args.vcffile[0])
    vcf_writer = vcf.Writer(args.output, vcf_reader)

    for filename in args.vcffile:
        vcf_reader = vcf.Reader(filename=filename)

        for record in vcf_reader:
            if not 'NormalReads' in record.INFO or \
                    not 'NormalEvidenceReads' in record.INFO:
                continue

            if 'LOWDEPTH' in record.FILTER:
                continue

            normdepth = int(record.INFO['NormalReads'][0])
            normevidence = sum([int(nr) for nr in record.INFO['NormalEvidenceReads']])

            impl_vaf = normevidence*1./normdepth
            p_het = scipy.stats.binom_test(normevidence, normdepth, max(args.errorrate, impl_vaf)) 
            p_hom = scipy.stats.binom_test(normevidence, normdepth, max(.75, impl_vaf)) 

            if p_het > 1.-args.alpha:
                record.INFO['HET'] = 1.-p_het
            if p_hom > 1.-args.alpha:
                record.INFO['HOM'] = 1.-p_hom
            vcf_writer.write_record(record)

        return 0

if __name__ == "__main__":
    main()
