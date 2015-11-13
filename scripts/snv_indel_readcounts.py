#!/usr/bin/env python
"""
Routines for reading in paired bam-readcount files, 
and comparing to a somatic VCF
"""
import argparse
import vcf
import sys

def main():
    """
    Driver program - Read in a VCF file and normal/tumour read counts for each base at each position,
    and output read counts at each call
    """
    parser = argparse.ArgumentParser(description='Validate somatic calls against bamcounts of normal/tumor BAMs')
    parser.add_argument('vcffile', help='Name of vcf file to validate')
    parser.add_argument('normal_bamcounts', help='Name of normal bamcounts file')
    parser.add_argument('tumour_bamcounts', help='Name of tumour bamcounts file')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, 
                        help='Output VCF (default: stdout)')
    args = parser.parse_args()

    normal_rc = BamReadcountFile(args.normal_bamcounts)
    tumour_rc = BamReadcountFile(args.tumour_bamcounts)

    vcf_reader = vcf.Reader(open(args.vcffile,'r'))
    vcf_writer = vcf.Writer(args.output, vcf_reader)

    for record in vcf_reader:
        alt = str(record.ALT[0])
        record.INFO['NormalReads'] = normal_rc.reads(record.CHROM, record.POS, alt)
        record.INFO['NormalEvidenceReads'] = ','.join([str(x) for x in normal_rc.evidence_reads(record.CHROM, record.POS, alt)])
        record.INFO['TumourReads'] = tumour_rc.reads(record.CHROM, record.POS, alt)
        record.INFO['TumourEvidenceReads'] = ','.join([str(x) for x in tumour_rc.evidence_reads(record.CHROM, record.POS, alt)])

        vcf_writer.write_record(record)

    return 0

class BamReadcountFile(object):
    """
    Contains the data for a bam-readcount file
    """
    @classmethod
    def chrom_pos_tuple(cls, chrom, pos=None):
        if pos is None:
            return chrom
        else:
            return chrom, pos

    def __init__(self, filename):
        self.__nullentry = ReadCountEntry()
        self.__entries = {}
        with open(filename, 'r') as f:
            for line in f:
                entry = ReadCountEntry(line)
                self.__entries[(entry.chrom, entry.pos)] = entry

    def __getitem__(self, chrom_pos):
        " Returns the read count entry at this position "
        if not chrom_pos in self.__entries:
            return self.__nullentry
        return self.__entries[chrom_pos]

    def __contains__(self, chrom_pos):
        " Returns the read count entry at this position "
        return chrom_pos in self.__entries

    def depth(self, chrom, pos=None):
        " Returns the read depth "
        chrom_pos = self.chrom_pos_tuple(chrom, pos)
        return self[chrom_pos].depth

    def dp8(self, chrom, pos=None):
        " Returns the 8-depth (ACGT forward, ACGT reverse) at this position "
        chrom_pos = self.chrom_pos_tuple(chrom, pos)
        return self.__entries[chrom_pos].dp8

    def dp8_str(self, chrom, pos=None):
        " Returns the 8-depth (ACGT forward, ACGT reverse) at this position "
        chrom_pos = self.chrom_pos_tuple(chrom, pos)
        return ''.join( [ str(x) for x in self.dp8 ] )

    def evidence_reads(self, chrom, pos, alt):
        """
        Minimum number of reads supporing 'alt' at chrom, pos
        """
        npairs = [self[(chrom, pos+off)][base] for off, base in enumerate(alt)]
        forward = min([n[0] for n in npairs])
        reverse = min([n[1] for n in npairs])
        return forward, reverse

    def reads(self, chrom, pos, alt):
        """
        Minimum number of reads covering chrom,pos:pos+len(alt) 
        """
        nreads = [self.depth(chrom, p) for p in range(pos, pos+len(alt))]
        return min(nreads)


class ReadCountEntry(object):
    """
    Handles a position in a bam-readcount file
    """
    def __init__(self, line_string=""):
        self.__depths = {}
        if line_string == "":
            for base in ['A', 'C', 'G', 'T']:
                self.__depths[base] = BaseDepth()
            self.__chrom = '0'
            self.__pos = 0
            self.__ref = 'N'
            self.__totdepth = 0
        else:
            fields = line_string.split()
            if fields[0][0] != 'c':
                fields[0] = 'chr'+fields[0]
            self.__chrom = fields[0]
            self.__pos = int(fields[1])
            self.__ref = fields[2]
            self.__totdepth = fields[3]
            for entry in fields[4:]:
                self.__depths[entry[0]] = BaseDepth(entry)

    def __getitem__(self, base):
        " Returns the pair depth for a base at the current item "
        return self.__depths[base].depth_pair

    def __contains__(self, base):
        " Do we have data for base at this position? "
        return base in self.__depths[base]

    @property
    def depth(self):
        " Total depth at this position "
        return self.__totdepth

    @property
    def chrom(self):
        " Chromosome of this entry "
        return self.__chrom

    @property
    def pos(self):
        " Position for this entry "
        return self.__pos

    @property
    def dp8(self):
        " Returns the pair depth for a base at the current item "
        forwards = [ self.__depths[x].depth_pair[0] for x in ['A','C','G','T'] ]
        reverses = [ self.__depths[x].depth_pair[1] for x in ['A','C','G','T'] ]
        return tuple( forwards + reverses )

    @property
    def dp8_str(self):
        " Returns the pair depth for a base at the current item "
        counts = self.dp8
        return ','.join( [str(x) for x in counts] )

class BaseDepth(object):
    """
    Parses a field output from bam-readcount corresponding to base depth at a given location
    """
    def __init__(self, bam_readcounts_string="N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00"):
        fields = bam_readcounts_string.strip().split(':')
        self.__base = fields[0]
        self.__count = int(fields[1])
        self.__avg_map_qual = float(fields[2])
        self.__avg_base_qual = float(fields[3])
        self.__count_forward = int(fields[5])
        self.__count_reverse = int(fields[6])
        self.__avg_position = float(fields[7])
        self.__inputstr = bam_readcounts_string
        assert self.__count == self.__count_forward + self.__count_reverse

    @property
    def base(self):
        " Returns the base represented by this object "
        return self.__base

    @property
    def depth_pair(self):
        " Returns the forward and reverse depth for this base "
        return self.__count_forward, self.__count_reverse

    @property
    def depth(self):
        " Returns the total depth for this base "
        return self.__count_forward + self.__count_reverse

    @property
    def map_qual(self):
        " Returns the average map quality for reads supporting this base at this position"
        return self.__avg_map_qual

    @property
    def base_qual(self):
        " Returns the average base quality across reads with this base at this position"
        return self.__avg_base_qual

    @property
    def avg_frac_position(self):
        """
        Returns the average fractional position for this location across reads which support this base at this position
        """
        return self.__avg_position

    def __repr__(self):
        return "BaseDepth("+self.__inputstr+")"

if __name__ == "__main__":
    main()
