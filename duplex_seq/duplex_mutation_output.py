#!/usr/bin/env python

# USAGE: duplex_mutation_output.py
# CODED BY: John Letaw

from collections import OrderedDict
import argparse
import pysam

VERSION = '0.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input_bam', help='Input BAM')
    parser.add_argument('--input_bed', help='Input BED')
    parser.add_argument('--ref', help='Reference Genome')
    parser.add_argument('--outfile', help='Output Mutation Stats')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class DuplexBed():
    """
    1:19981652-19981652	1	19981652		C	T	tb133_NBL1	tb133_NBL1
    """
    def __init__(self, args):
        self.bed = {}
        self.fname = args.input_bed
        self.args = args
        self._parse_bed()

    def _parse_bed(self):
        """
        Arrange entries in the line.
        :return:
        """
        with open(self.fname, 'rU') as my_bed:
            for entry in my_bed:
                entry = entry.rstrip('\n').split()
                this_entry = BedEntry(entry, self.args.input_bam, self.args.ref)
                if this_entry.uniq_key not in self.bed:
                    self.bed[this_entry.uniq_key] = this_entry.uniq_val

class BedEntry():
    """
    Represent an entry in the BED file.
    """
    def __init__(self, bed_entry, inbam, ref_gen):
        self.bed_entry = bed_entry
        self.inbam = inbam
        self.ref_gen = ref_gen
        self.chrom = self.bed_entry[1]
        self.coord = self.bed_entry[2]
        self.ref = self.bed_entry[3]
        self.alt = self.bed_entry[4]
        self.design = self.bed_entry[5]
        self.gene = self.bed_entry[6]
        self.uniq_key = (self.chrom, self.coord, self.ref, self.alt)
        self.ival = self._create_ival_fmt()
        print(self._mpileup_grab())
        if self._mpileup_grab() != ['']:
            self.mpileup = MpileupData(self._mpileup_grab(), self.alt)
            self.uniq_val = [self.chrom, self.coord, self.ref, self.alt, self.design, self.gene, str(self.mpileup.depth),
                             str(self.mpileup.num_wt), str(self.mpileup.num_mut), str(self.mpileup.num_alt_non_mut),
                             str(self.mpileup.num_a), str(self.mpileup.num_t), str(self.mpileup.num_c),
                             str(self.mpileup.num_g), str(self.mpileup.num_n), str(self.mpileup.percent_wt),
                             str(self.mpileup.percent_mut), str(self.mpileup.percent_alt_non_mut),
                             str(self.mpileup.percent_a), str(self.mpileup.percent_t), str(self.mpileup.percent_g),
                             str(self.mpileup.percent_c), str(self.mpileup.percent_n)]
        else:
            print("HEREH")
            self.uniq_val = [self.chrom, self.coord, self.ref, self.alt, self.design, self.gene]

    def _create_ival_fmt(self):
        """
        Create a chrom:start-stop string, 0-based start, for samtools.
        :return:
        """
        start = str(int(self.coord) - 1)
        return ''.join([self.chrom, ':', self.coord, '-', self.coord])

    def _mpileup_grab(self):
        """
        Get information on this locus from samtools mpileup.
        :return:
        """
        return pysam.mpileup(self.inbam, "-f", self.ref_gen, "-r", self.ival).rstrip('\n').split('\t')


class MpileupData():
    """
    Everything derived from the mpileup string.
    where each line consists of chromosome, 1-based coordinate, reference base, the number of reads covering the
    site, read bases and base qualities. At the read base column, a dot stands for a match to the reference base
    on the forward strand, a comma for a match on the reverse strand, `ACGTN' for a mismatch on the forward strand
    and `acgtn' for a mismatch on the reverse strand. A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an
    insertion between this reference position and the next reference position. The length of the insertion is given
    by the integer in the pattern, followed by the inserted sequence.
    """
    def __init__(self, mpileup, mut):
        self.mpileup = mpileup
        self.mut = mut
        self.chrom = mpileup[0]
        self.coord = mpileup[1]
        self.ref = mpileup[2]
        self.depth = int(mpileup[3])
        self.quals = mpileup[4]
        self.seq = mpileup[5]
        self.num_n = self.seq.count('N')
        self.num_a = self.seq.count('A')
        self.num_c = self.seq.count('C')
        self.num_g = self.seq.count('G')
        self.num_t = self.seq.count('T')
        self.num_wt = self.seq.count('.') + self.seq.count(',')
        self.num_mut = self.seq.count(self.mut)
        self.num_alt_non_mut = self.num_a + self.num_c + self.num_g + self.num_t - self.num_mut
        self.percent_a = float((self.num_a / (self.depth - self.num_n + 1)) * 100)
        self.percent_c = float((self.num_c / (self.depth - self.num_n + 1)) * 100)
        self.percent_g = float((self.num_g / (self.depth - self.num_n + 1)) * 100)
        self.percent_t = float((self.num_t / (self.depth - self.num_n + 1)) * 100)
        self.percent_n = float((self.num_n / self.depth) * 100)
        self.percent_alt_non_mut = float((self.num_alt_non_mut / (self.depth - self.num_n + 1)) * 100)
        self.percent_wt = float((self.num_wt / (self.depth - self.num_n + 1)) * 100)
        self.percent_mut = float((self.num_mut / (self.depth - self.num_n + 1)) * 100)

def write_output(outfile, bed):
    """
    Write the mutation output.
    :return:
    """
    header = ['Chromosome', 'Position', 'Ref_Base', 'Mut_Base', 'Design_Target', 'GeneID', 'Depth', '#WT', '#Mut',
              '#Alt_non_mut_of_interest', '#A', '#T', '#G', '#C', '#N', 'Percent_WT', 'Percent_Mut',
              'Percent_Alt_non_mut_of_interest', 'Percent_A', 'Percent_T', 'Percent_G', 'Percent_C', 'Percent_N',
              'Sequence']
    with open(outfile, 'w') as to_write:
        to_write.write('\t'.join(header))
        to_write.write('\n')
        for v in sorted(bed.bed.values(), key=lambda x: (int(x[0]), int(x[1]))):
            to_write.write('\t'.join(v))
            to_write.write('\n')

def main():
    args = supply_args()
    my_bed = DuplexBed(args)
    write_output(args.outfile, my_bed)

if __name__ == "__main__":
    main()


