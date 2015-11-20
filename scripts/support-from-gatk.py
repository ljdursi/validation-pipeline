#!/usr/bin/env python
"""
Count the support for structural variation calls in a pair of bam files
"""
import vcf
import pysam
import argparse
import sys

def get_sample_name(bam_fn):
    bam_fh = pysam.AlignmentFile(bam_fn, 'rb')
    header = bam_fh.header
    sample_name = ""
    for l in header["PG"]:
        if "SM" in l:
            if sample_name != "":
                assert(l["SM"] == sample_name)
            else:
                sample_name = l["SM"]

    assert(sample_name != "")
    return(sample_name)

# Main
parser = argparse.ArgumentParser(description='Annotate calls with after running GATK HaplotypeCaller')
parser.add_argument('original_vcf', help='Name of vcf file to validate')
parser.add_argument('gatk_vcf', help='Name of vcf file annotated by GATK')
parser.add_argument('normal_bam', help='Name of bam file for the normal genome')
parser.add_argument('tumour_bam', help='Name of bam file for the tumour genome')
parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, 
                        help='Output VCF (default: stdout)')
args = parser.parse_args()

# Read the original VCF file into a dict
selected_sites = dict()
vcf_reader = vcf.Reader(open(args.original_vcf))

for record in vcf_reader:
    key = record.CHROM + str(record.POS)
    selected_sites[key] = record

# Get the sample name for the tumor and normal from the bam headers
sample_names = dict()
sample_names[get_sample_name(args.normal_bam)] = "Normal"
sample_names[get_sample_name(args.tumour_bam)] = "Tumour"

vcf_writer = vcf.Writer(args.output, vcf_reader)

# Annotate each site
gatk_vcf_reader = vcf.Reader(open(args.gatk_vcf))

for gatk_record in gatk_vcf_reader:
    key = gatk_record.CHROM + str(gatk_record.POS)
    if key not in selected_sites:
        continue
    original_record = selected_sites[key]

    for sample_name, sample_type in sample_names.items():

        # defaults
        evidence_str = "0,0"
        total_str = "0"

        call = gatk_record.genotype(sample_name).data
        if call.GT != None:
            assert(len(call.SAC) == 4)
            evidence_str = ",".join(str(x) for x in call.SAC[2:4])
            total_str = str(call.AD[0] + call.AD[1])
    
        original_record.INFO[sample_type + 'EvidenceReads'] = evidence_str
        original_record.INFO[sample_type + 'TotalReads'] = total_str
        
for record in selected_sites:
    vcf_writer.write_record(selected_sites[record])
