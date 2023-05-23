#!/usr/bin/env python

import argparse

VERSION = '0.0.1'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('raw_junc', help='Input chimeric junctions from STAR.')
    parser.add_argument('srch_junc', help='Input BEDPE file describing regions of interest.')
    parser.add_argument('bedpe_out', help='Output containing sites of interest.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class ChimJuncRec:
    """
    Set cotegory names for records.
    """
    def __init__(self, junc):
        self.chrom1 = junc[0]
        self.coord1 = junc[1]
        self.chrom2 = junc[3]
        self.coord2 = junc[4]
        self.rec = (self.chrom1, self.coord1, self.chrom2, self.coord2)


class ChimJuncFile:
    """
    HEADER:
    chr_donorA      brkpt_donorA    strand_donorA   chr_acceptorB   brkpt_acceptorB strand_acceptorB
    junction_type   repeat_left_lenA        repeat_right_lenB       read_name       start_alnA
    cigar_alnA      start_alnB      cigar_alnB      num_chim_aln    max_poss_aln_score      non_chim_aln_score
    this_chim_aln_score     bestall_chim_aln_score  PEmerged_bool   readgrp
    """
    def __init__(self, filename, regions=None):
        self.filename = filename
        self.regions = regions
        self.juncs = self._parse_cj()

    def _parse_cj(self) -> list:
        """
        Get relevant entries from raw STAR chim junctions.  Skip header, ignore entries that start with comment (#).
        :return:
        """
        juncs = []
        with open(self.filename, 'r') as myfile:
            next(myfile)
            for entry in myfile:
                if not entry.startswith('#'):
                    entry = ChimJuncRec(entry.rstrip('\n').split('\t'))
                    if self.regions:
                        if entry.chrom1 in self.regions.valid_chrs:
                            juncs.append(entry)
                    else:
                        juncs.append(entry)
            return juncs

    def filt_cj(self):
        """
        Get all of the junction records that lie within the BEDPE regions of interest.
        :return:
        """
        filt_juncs = {}
        for reg in self.regions.regions:
            if reg.name not in filt_juncs:
                filt_juncs[reg.name] = {}
            for junc in self.juncs:
                if reg.chrom1 == junc.chrom1:
                    if junc.coord1 > reg.start1 and junc.coord1 <= reg.end1:
                        if junc.coord2 > reg.start2 and junc.coord2 <= reg.end2:
                            if junc.rec not in filt_juncs[reg.name]:
                                filt_juncs[reg.name][junc.rec] = 1
                            else:
                                filt_juncs[reg.name][junc.rec] += 1
        return filt_juncs


class BedPeRec:
    """
    Set cotegory names for records.
    """
    def __init__(self, junc):
        self.chrom1 = junc[0]
        self.start1 = junc[1]
        self.end1 = junc[2]
        self.chrom2 = junc[3]
        self.start2 = junc[4]
        self.end2 = junc[5]
        self.name = junc[6]
        self.rec = '\t'.join(junc)


class BedPe:
    """
    https://bedtools.readthedocs.io/en/latest/content/general-usage.html
    chrom1
    start1 (0-based)
    end1
    chrom2
    start2 (0-based)
    end2
    name (.)
    score (.)
    strand1 (.)
    strand2 (.)
    """
    def __init__(self, filename):
        self.filename = filename
        self.regions = self._parse_bedpe()
        self.valid_chrs = [x.chrom1 for x in self.regions]

    def _parse_bedpe(self):
        """
        Get the BEDPE entries.
        :return:
        """
        regions = []
        with open(self.filename, 'r') as myfile:
            for entry in myfile:
                entry = entry.rstrip('\n').split('\t')
                regions.append(BedPeRec(entry))
        return regions


class Output:
    """
    Write the output.  Current headings:
    Region identifier
    Chrom1
    Start1
    Chrom2
    Start2
    Count
    """
    def __init__(self, filename, juncs):
        self.filename = filename
        self.header = ['Fusion', 'Chrom1', 'Coord1', 'Chrom2', 'Coord2', 'Count']
        self.juncs = juncs

    def write_me(self, thresh=10):
        with open(self.filename, 'w') as handle_out:
            handle_out.write('\t'.join(self.header))
            handle_out.write('\n')
            for fusion in self.juncs:
                for coords, count in self.juncs[fusion].items():
                    if count >= thresh:
                        to_write = [fusion, coords[0], coords[1], coords[2], coords[3], str(count)]
                        handle_out.write('\t'.join(to_write))
                        handle_out.write('\n')
        handle_out.close()


def main():
    args = supply_args()
    my_bed = BedPe(args.srch_junc)
    my_juncs = ChimJuncFile(args.raw_junc, my_bed).filt_cj()
    Output(args.bedpe_out, my_juncs).write_me()


if __name__ == "__main__":
    main()
