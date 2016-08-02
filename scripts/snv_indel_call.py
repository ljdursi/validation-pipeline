#!/usr/bin/env python
"""
VCF filter; make calls based on supporting depth
"""
from __future__ import print_function
import argparse
import vcf
import vcf.parser
import sys
import numpy
import scipy.stats

def call_from_depths(totdepth, evidence, error_rate, prob_threshold):
    """
    Null hypothesis: counts are consistent with a binomial probabilty of error_rate
                     
    Returns:    True if we can reject the null hypothesis w/ probability 1-prob_threshold
                False otherwise
    """
    altdepth = sum(evidence)
    prob = 1.-scipy.stats.binom.cdf(altdepth, totdepth, error_rate)
    return prob < prob_threshold

def reject_from_strandbias(totdepth, evidence, prob_threshold, mindepth=30):
    if len(evidence) == 1:
        return False  # can't reject; don't have strand information

    if sum(evidence) < mindepth:
        return False

#    forward, reverse = evidence
#    if min([forward,reverse]) <= strand_bias*(forward+reverse):
#        return True
    prob = 1.-scipy.stats.binom.cdf(max(evidence), sum(evidence), 0.66)
    return prob < prob_threshold

def germline_hom_het(n_depth, n_evidence, t_depth, t_evidence, alpha):
    # consistent w/ homozygous or het?  True
    if scipy.stats.binom_test(sum(n_evidence), n_depth, 1.) >= alpha:
        return True
    if scipy.stats.binom_test(sum(n_evidence), n_depth, .95) >= alpha:
        return True
    if scipy.stats.binom_test(sum(n_evidence), n_depth, 0.5) >= alpha:
        return True
    if scipy.stats.binom_test(sum(n_evidence), n_depth, 0.45) >= alpha:
        return True
    return False


def reject_from_normal_evidence(n_depth, n_evidence, t_depth, t_evidence, alpha):
    """
    Null hypothesis: germline results are consistent with being drawn from same
                     binomial distribution as tumour.

    Returns: True if we cannot reject the null w/ probability 1-alpha
             False otherwise

    Use fisher exact test to determine if these are consistent
    """
    phat_tumour = sum(t_evidence)*1./t_depth
    phat_normal = sum(n_evidence)*1./n_depth
    if phat_normal == 0.:
        return False

    if phat_normal >= phat_tumour/2:
        return True

    oddsratio, prob = scipy.stats.fisher_exact([[sum(t_evidence)/2, sum(n_evidence)], [t_depth-sum(t_evidence)/2, n_depth-sum(n_evidence)]])
    if prob < alpha:
        return False
    else:
        return True
#    phat = sum(t_evidence + n_evidence)*1./(t_depth + n_depth)
#
#    z = (phat_normal-phat_tumour)/numpy.sqrt(phat * (1-phat) * (1./n_depth + 1./t_depth))
#    if scipy.stats.norm.cdf(z) < alpha:
#        return False
    return True


def filter_calls():
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
    vcf_reader.infos['Validation_status'] = vcf.parser._Info(id='Validation_status', num='.', type='String',
                                                             desc='Status from validation data',
                                                             source=None, version=None)
    vcf_writer = vcf.Writer(args.output, vcf_reader)

    for record in vcf_reader:
        normal_reads = int(record.INFO['NormalReads'][0])
        tumour_reads = int(record.INFO['TumourReads'][0])
        normal_evidence = [int(x) for x in record.INFO['NormalEvidenceReads']]
        tumour_evidence = [int(x) for x in record.INFO['TumourEvidenceReads']]

        if min(normal_reads, tumour_reads) < args.mindepth:
            record.FILTER = ['LOWDEPTH']
            record.INFO['Validation_status'] = 'LOWDEPTH'
            vcf_writer.write_record(record)
            continue

        record.FILTER = []
        if not call_from_depths(tumour_reads, tumour_evidence, args.error, args.callthreshold) or sum(tumour_evidence) < 7:
            record.FILTER += ['NOTSEEN']
        if (tumour_reads > 0) > 0 and reject_from_strandbias(tumour_reads, tumour_evidence, args.strandbias):
            record.FILTER += ['STRANDBIAS']
        if not(normal_reads == 0 or tumour_reads == 0) and germline_hom_het(normal_reads, normal_evidence, tumour_reads, tumour_evidence, args.germlineprob):
            record.FILTER += ['GERMLINE']
        elif call_from_depths(normal_reads, normal_evidence, args.germlineprob, args.callthreshold):
            record.FILTER += ['NORMALEVIDENCE']
        elif not(normal_reads == 0 or tumour_reads == 0) and reject_from_normal_evidence(normal_reads, normal_evidence, tumour_reads, tumour_evidence, args.germlineprob):
            record.FILTER += ['NORMALEVIDENCE']
        if len(record.FILTER) == 0:
            record.FILTER = ['PASS']
        record.INFO['Validation_status'] = ','.join(record.FILTER)
        vcf_writer.write_record(record)

if __name__ == "__main__":
    filter_calls()
