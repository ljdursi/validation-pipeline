#!/usr/bin/env python
"""
VCF filter; adjust genotypes based on DP4 (4-depth) scores
"""
from __future__ import print_function
import argparse
import vcf
import sys
import scipy.stats

def call_from_dp4(dp4, error_rate, prob_threshold, strand_bias):
    totdepth = sum(dp4)
    refdepth = sum(dp4[0:2])
    altdepth = sum(dp4[2:4])

    if min(dp4[2:3]) < strand_bias*altdepth:
        return False
    
    prob = 1.-scipy.stats.binom.cdf(altdepth, totdepth, error_rate)
    return prob < prob_threshold


def filter_genotypes():
    parser = argparse.ArgumentParser( description='Set genotypes based on DP4 scores')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-e', '--error', type=float, default=0.01, help='Error rate')
    parser.add_argument('-t', '--callthreshold', type=float, default=0.02, help='Max prob to call')
    parser.add_argument('-s', '--strandbias', type=float, default=0.10, help='minimum strand ratio')
    parser.add_argument('-m', '--mindepth', type=int, default=10, help='minimum total depth')
    args = parser.parse_args()

    vcf_reader = vcf.Reader(args.input)
    vcf_writer = vcf.Writer(args.output, vcf_reader)

    for record in vcf_reader:
        if record.INFO['DP'] < args.mindepth:
            continue

        for sample in record.samples:
            varcall = call_from_dp4(sample.data.DP4, args.error, args.callthreshold, args.strandbias)
            if varcall and sample.data.GT == '0/0':
                sample.data = sample.data._replace(GT='0/1')
        vcf_writer.write_record(record)

if __name__ == "__main__":
    filter_genotypes()
