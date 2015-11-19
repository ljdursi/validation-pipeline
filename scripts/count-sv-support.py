#!/usr/bin/env python
"""
Count the support for structural variation calls in a pair of bam files
"""
import pysam
import vcf
import argparse
import sys

MAX_DISTANCE=1000000000
WINDOW=300
CLIP_WINDOW=5


# returns the reference positions at which the read is softclipped
def get_softclip_positions(read):
    out = list()
    ap = read.get_aligned_pairs(matches_only=True)
    if len(ap) == 0:
        return out

    # Check for softclip at the start of the read
    if ap[0][0] != 0:
        clip_start_ref = ap[0][1]
        clip_length = ap[0][0]
        #print '5prime clip', clip_start_ref, clip_length
        out.append(clip_start_ref)
           
    # Check for softclip at the end of the read
    if ap[-1][0] != read.query_length - 1:
        clip_start_ref = ap[-1][1]
        clip_length = read.query_length - ap[-1][0]
        #print '3prime clip', clip_start_ref, clip_length
        out.append(clip_start_ref)

    return out

# count the number of reads that are softclipped near the provided breakpoint
def count_clipped_reads_at_breakpoint(bam_fh, bp_chr, bp_pos):
    evidence_count = 0
    
    for read in bam_fh.fetch(bp_chr, bp_pos - WINDOW, bp_pos + WINDOW):

        if read.is_secondary or read.is_supplementary:
            continue

        for clip_pos in get_softclip_positions(read):
            if abs(clip_pos - bp_pos) < CLIP_WINDOW:
                evidence_count += 1

    return evidence_count

# count the number of paired reads that are roughly concordant with the SV call
def count_paired_reads_for_breakends(bam_fh, record):

    # the coordinates of the first breakend
    chr_1 = record.CHROM
    pos_1 = record.POS

    assert(len(record.ALT) == 1)
    
    # the coordinates of the second breakend
    alt = record.ALT[0]
    chr_2 = alt.chr
    pos_2 = alt.pos
 
    pe_evidence_count = 0
    total_pairs_count = 0
    
    # Establish the the orientation of the breakpoint
    bp_opp_orientation = alt.orientation != alt.remoteOrientation
    
    # do not double count a read and its pair if they are both within the window 
    observed_reads = dict()

    for read in bam_fh.fetch(chr_1, pos_1 - WINDOW, pos_1 + WINDOW):
        
        # skip non-primary alignments
        if read.is_secondary or read.is_supplementary:
            continue

        # only count each pair once
        if read.qname in observed_reads:
            continue

        # calculate distance from the pair to the second breakend if its on the same chr
        mate_dist_to_bp_2 = MAX_DISTANCE
        if bam_fh.references[read.next_reference_id] == chr_2:
            mate_dist_to_bp_2 = pos_2 - read.next_reference_start

        # check orientation
        read_opp_orientation = read.is_reverse != read.mate_is_reverse

        # We consider a read pair to be "evidence" if these conditions are met
        # - it is within the validation pulldown region of the second half of the breakend
        # - it is in the expected orientation
        if abs(mate_dist_to_bp_2) < WINDOW and read_opp_orientation == bp_opp_orientation:
            pe_evidence_count += 1
        
        total_pairs_count += 1
        observed_reads[read.qname] = 1
    
    return pe_evidence_count, total_pairs_count

# Main
def main():
    parser = argparse.ArgumentParser(description='Validate somatic SV calls against normal/tumor BAMs')
    parser.add_argument('vcffile', help='Name of vcf file to validate')
    parser.add_argument('normal_bam', help='Name of normal bam file')
    parser.add_argument('tumour_bam', help='Name of tumour bam file')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, 
                            help='Output VCF (default: stdout)')
    args = parser.parse_args()

    normal_bam_fh = pysam.AlignmentFile(args.normal_bam, 'rb')
    tumour_bam_fh = pysam.AlignmentFile(args.tumour_bam, 'rb')
    vcf_reader = vcf.Reader(open(args.vcffile))

    bam_fh_list = [ normal_bam_fh, tumour_bam_fh ]
    name_list = [ "Normal", "Tumour" ]

    # constants
    vcf_out_fh = vcf.Writer(args.output, vcf_reader)

    for record in vcf_reader:

        # Iterate over the reads mapped to the first half of the breakend
        for bam_fh, name in zip(bam_fh_list, name_list):
            spanned_evidence_count, total_pair_count = count_paired_reads_for_breakends(bam_fh, record)
            bp_1_reads = count_clipped_reads_at_breakpoint(bam_fh, record.CHROM, record.POS)
            bp_2_reads = count_clipped_reads_at_breakpoint(bam_fh, record.ALT[0].chr, record.ALT[0].pos)

            record.INFO[name + 'EvidenceReads'] = spanned_evidence_count
            record.INFO[name + 'Reads'] = total_pair_count
            record.INFO[name + 'Bp1ClipEvidence'] = bp_1_reads
            record.INFO[name + 'Bp2ClipEvidence'] = bp_2_reads

        # pysam puts an empty key if there is a single INFO tag ending in a delimiter. get rid of it
        del record.INFO['']
        vcf_out_fh.write_record(record)

if __name__ == "__main__":
    main()
