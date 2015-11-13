#!/usr/bin/env python
"""
VCF filter; adjust genotypes based on DP4 (4-depth) scores
"""
from __future__ import print_function
import argparse
import vcf
import sys
import scipy.stats

def call_from_depths(totdepth, evidence, error_rate, prob_threshold, strand_bias):
    forward, reverse = evidence
    altdepth = forward+reverse

    if min([forward,reverse]) <= strand_bias*(forward+reverse):
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
        normal_reads = int(record.INFO['NormalReads'][0])
        tumour_reads = int(record.INFO['TumourReads'][0])
        normal_evidence = [int(x) for x in record.INFO['NormalEvidenceReads']]
        tumour_evidence = [int(x) for x in record.INFO['TumourEvidenceReads']]

        if min(normal_reads, tumour_reads) < args.mindepth:
            continue

        if call_from_depths(normal_reads, normal_evidence, args.error, args.callthreshold, args.strandbias):
            record.FILTER = ['GERMLINE']
        elif call_from_depths(tumour_reads, normal_evidence, args.error, args.callthreshold, args.strandbias):
            record.FILTER = ['PASS']
        else:
            record.FILTER = ['PASS']
        vcf_writer.write_record(record)

if __name__ == "__main__":
    filter_genotypes()
