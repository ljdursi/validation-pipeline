#!/usr/bin/env python
"""
VCF filter; make calls based on supporting depth
"""
from __future__ import print_function
import argparse
import vcf
import sys
import numpy
import scipy.stats

def call_from_depths(totdepth, evidence, error_rate, prob_threshold):
    """
    Null hypothesis: counts are consistent with a binomial probabilty of error_rate
                     
    Returns:    True if we can reject the null phyothesis w/ probability 1-prob_threshold
                False otherwise
    """
    altdepth = sum(evidence)
    prob = 1.-scipy.stats.binom.cdf(altdepth, totdepth, error_rate)
    return prob < prob_threshold

def reject_from_strandbias(totdepth, evidence, strand_bias):
    if len(evidence) == 1:
        return False  # can't reject; don't have strand information

    forward, reverse = evidence
    if min([forward,reverse]) <= strand_bias*(forward+reverse):
        return True
    return False

def germline_evidence(n_depth, n_evidence, t_depth, t_evidence, alpha):
    """
    Null hypothesis: germline results are consistent with being drawn from same
                     binomial distribution as tumour.

    Returns: True if we cannot reject the null w/ probability 1-alpha
             False otherwise

    Use binomial proportion test, eg http://itl.nist.gov/div898/software/dataplot/refman1/auxillar/binotest.htm
    """
    phat_tumour = sum(t_evidence)*1./t_depth
    phat_normal = sum(n_evidence)*1./n_depth
    if phat_normal == 0.:
        return False

    if phat_normal >= phat_tumour/2:
        return True

    # consistent w/ homozygous or het?  True
    if scipy.stats.binom_test(sum(n_evidence), n_depth, 1.) >= alpha:
        return True
    if scipy.stats.binom_test(sum(n_evidence), n_depth, 0.5) >= alpha:
        return True

    # consistent w/ phat_tumour, w/in factor of 2?  True
    p1 = scipy.stats.binom_test(sum(n_evidence), n_depth, phat_tumour) 
    p2 = scipy.stats.binom_test(sum(n_evidence), n_depth, phat_tumour/2)
    p = max(p1, p2)

    if p >= alpha:
        return True
    return False
#    phat = sum(t_evidence + n_evidence)*1./(t_depth + n_depth)
#
#    z = (phat_normal-phat_tumour)/numpy.sqrt(phat * (1-phat) * (1./(n_depth+.0001) + 1./(t_depth+.0001)) + 0.0001)
#    if scipy.stats.norm.cdf(z) < alpha:
#        return False
#    return True


def filter_genotypes():
    parser = argparse.ArgumentParser( description='Set genotypes based on DP4 scores')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-e', '--error', type=float, default=0.01, help='Error rate')
    parser.add_argument('-t', '--callthreshold', type=float, default=0.02, help='Max prob to call')
    parser.add_argument('-s', '--strandbias', type=float, default=0.10, help='minimum strand ratio')
    parser.add_argument('-m', '--mindepth', type=int, default=10, help='minimum total depth')
    parser.add_argument('-g', '--germlineprob', type=float, default=0.02, help='Maximum prob of germline')
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

        record.FILTER = []
        if not call_from_depths(tumour_reads, tumour_evidence, args.error, args.callthreshold):
            record.FILTER += ['NOTSEEN']
        if reject_from_strandbias(tumour_reads, tumour_evidence, args.strandbias):
            record.FILTER += ['STRANDBIAS']
        if not(normal_reads == 0 or tumour_reads == 0) and germline_evidence(normal_reads, normal_evidence, tumour_reads, tumour_evidence, args.germlineprob):
            record.FILTER += ['GERMLINE']
        if len(record.FILTER) == 0:
            record.FILTER = ['PASS']
        vcf_writer.write_record(record)

if __name__ == "__main__":
    filter_genotypes()
