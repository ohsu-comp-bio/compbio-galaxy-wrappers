#!/usr/bin/env python

"""
Used on STAR chimeric junction output to try and find evidence of partial tandem duplication events.

0.0.1 - Initial commit
0.0.2 - Add gene name columns to output
0.1.0 - Allow for manual input of regions of interest
"""

import argparse
import json
import pysam
import re

VERSION = '0.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('raw_junc', help='Input chimeric junctions from STAR.')
    parser.add_argument('ref_fasta', help='Reference FASTA, FAI index should be in same location.')
    parser.add_argument('samp_met', help='Sample metrics JSON file with total_on_target_transcripts metric.')
    parser.add_argument('bedpe_out', help='Output containing sites of interest.')
    parser.add_argument('--srch_junc', help='Input BEDPE file describing regions of interest.')
    parser.add_argument('--manual_mode', action='store_true', help='Process manually entered coordinates.')
    parser.add_argument('--chr1', help='')
    parser.add_argument('--start1', type=int, help='')
    parser.add_argument('--end1', type=int, help='')
    parser.add_argument('--chr2', help='')
    parser.add_argument('--start2', type=int, help='')
    parser.add_argument('--end2', type=int, help='')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    if args.manual_mode:
        assert args.chr1 and args.start1 and args.end1 and args.chr2 and args.start2 and args.end2
    else:
        assert not (args.chr1 and args.start1 and args.end1 and args.chr2 and args.start2 and args.end2)

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

    def filt_cj(self, regions):
        """
        Get all of the junction records that lie within the BEDPE regions of interest.
        :return:
        """
        filt_juncs = {}
        for reg in regions:
            if reg.rec not in filt_juncs:
                filt_juncs[reg.rec] = {}
            for junc in self.juncs:
                if reg.chrom1 == junc.chrom1:
                    if int(junc.coord1) > int(reg.start1) and int(junc.coord1) <= int(reg.end1):
                        if int(junc.coord2) > int(reg.start2) and int(junc.coord2) <= int(reg.end2):
                            if junc.rec not in filt_juncs[reg.rec]:
                                filt_juncs[reg.rec][junc.rec] = 1
                            else:
                                filt_juncs[reg.rec][junc.rec] += 1
        return filt_juncs


class ManualRec:
    """
    Set cotegory names for records.
    """
    def __init__(self, junc):
        self.chrom1 = junc[0]
        self.start1 = junc[1] - 1
        self.end1 = junc[2]
        self.chrom2 = junc[3]
        self.start2 = junc[4] - 1
        self.end2 = junc[5]
        self.valid_chrs = self.chrom1
        self.rec = (self.chrom1, self.start1, self.end1, self.chrom2, self.start2, self.end2)


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
        self.refseq1 = junc[10]
        self.ccds1 = junc[11]
        self.exon1 = junc[12]
        self.refseq2 = junc[13]
        self.ccds2 = junc[14]
        self.exon2 = junc[15]
        self.gene1, self.gene2 = self._parse_gene_from_name(self.name)
        self.rec = (self.name, self.refseq1, self.ccds1, self.exon1, self.refseq2, self.ccds2, self.exon2,
                    self.gene1, self.gene2)
        self.ref_coords = self._ref_coord()

    @staticmethod
    def _parse_gene_from_name(fusion_name):
        """
        Get the HGNC gene names from the fusion name as defined in the input BEDPE.
        :param fusion_name:
        :return:
        """
        name_ptrn = re.compile(r'([A-Za-z0-9]*)[ei][0-9]{1,2}-([A-Za-z0-9]*)[ei][0-9]{1,2}')
        r = re.search(name_ptrn, fusion_name)
        gene1 = r.group(1)
        gene2 = r.group(2)
        return gene1, gene2

    def _ref_coord(self):
        """
        Figure out what the reference junction coordinates are for a given pair.  For instance, if the BEDPE says:
        chr11	118361909	118362034	chr11	118342375	118345031	KMT2Ae14-KMT2Ae3	.	.	.
        Reference would be:
        chr11	118362033	118362034	chr11	118342375	118342376	KMT2Ae14-KMT2Ae3-ref	.	.	.
        Output:
        KMT2Ae14-KMT2Ae3-ref	chr11	118362034	chr11	118342376	145
        :return:
        """
        return int(self.end1), int(self.start2)+1


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
        self.ref1_coords = []
        self.ref2_coords = []
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
                rec = BedPeRec(entry)
                regions.append(rec)
                self.ref1_coords.append(rec.ref_coords[0])
                self.ref2_coords.append(rec.ref_coords[1])
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
    def __init__(self, filename, juncs, bedpe, ref_fasta, tott):
        self.filename = filename
        self.header = ['Fusion', 'Gene1', 'Chrom1', 'Coord1', 'RefSeq1', 'CCDS1', 'Exon1', 'RefStatus1', 'Gene2',
                       'Chrom2', 'Coord2', 'RefSeq2', 'CCDS2', 'Exon2', 'RefStatus2', 'Count', 'NormCount',
                       'CombinedSeq']
        self.juncs = juncs
        self.bedpe = bedpe
        self.ref_fasta = ref_fasta
        self.tott = tott

    @staticmethod
    def _apply_ref(coord, ref_coord):
        """
        Apply additional label to denote whether the junction is reference or not.
        :return:
        """
        if int(coord) in ref_coord:
            return "REF"
        else:
            return "NON_REF"

    def _get_combined_seq(self, chrom1, coord1, chrom2, coord2, size=20):
        """
        Get the sequence from the reference FASTA and create a combined sequence.
        :return:
        """
        end1 = int(coord1) - 1
        start1 = end1 - size - 1
        start2 = int(coord2)
        end2 = start2 + size + 1
        seq1 = self.ref_fasta.fetch(reference=chrom1, start=start1, end=end1)
        seq2 = self.ref_fasta.fetch(reference=chrom2, start=start2, end=end2)
        return seq1 + seq2

    def _calc_on_target(self, *args):
        """
        calculate on-target cpm from junction & spanning frag count and sample level metrics
        :return:
        """
        num = 0.0
        for arg in args:
            num += float(arg)
        j_s_cpm = (num / self.tott) * 1e6
        return str(round(j_s_cpm, 3))

    def write_me(self, thresh=10):
        with open(self.filename, 'w') as handle_out:
            handle_out.write('\t'.join(self.header))
            handle_out.write('\n')
            for fusion in self.juncs:
                for coords, count in self.juncs[fusion].items():
                    if self.bedpe:
                        if count >= thresh:
                            to_write = [fusion[0], fusion[7], coords[0], coords[1], fusion[1], fusion[2], fusion[3],
                                        self._apply_ref(coords[1], self.bedpe.ref1_coords),
                                        fusion[8], coords[2], coords[3], fusion[4], fusion[5], fusion[6],
                                        self._apply_ref(coords[3], self.bedpe.ref2_coords),
                                        str(count), self._calc_on_target(count),
                                        self._get_combined_seq(coords[0], coords[1], coords[2], coords[3])]
                    else:
                        to_write = ['.', '.', coords[0], coords[1], '.', '.', '.',
                                    '.',
                                    '.', coords[2], coords[3], '.', '.', '.',
                                    '.',
                                    str(count), self._calc_on_target(count),
                                    self._get_combined_seq(coords[0], coords[1], coords[2], coords[3])]
                    handle_out.write('\t'.join(to_write))
                    handle_out.write('\n')

        handle_out.close()


class SampMetJson:
    def __init__(self, filename):
        self.filename = filename
        self.tott = self._get_sample_met()

    def _get_sample_met(self):
        with open(self.filename, 'r') as samp_metrics_fh:
            metrics = json.load(samp_metrics_fh)
            for entry in metrics['sampleRunMetrics']:
                if entry['metric'] == 'total_on_target_transcripts':
                    return float(entry['value'])
        return None


class ManualOutput:
    def __init__(self, filename):
        self.filename = filename

    def write_me(self, counts):
        with open(self.filename, 'w') as handle_out:
            for entry in counts:
                print(entry)
                for region, cnt in sorted(counts[entry].items(), key=lambda item: item[1], reverse=True):
                    print(region)
                    to_write = [region[0], region[1], region[2], region[3], str(cnt)]
                    if int(region[1]) == entry[1] or int(region[1]) == entry[2]:
                        to_write.append('REF')
                    else:
                        to_write.append('.')
                    if int(region[3]) == entry[4] or int(region[3]) == entry[5]:
                        to_write.append('REF')
                    else:
                        to_write.append('.')
                    handle_out.write('\t'.join(to_write))
                    handle_out.write('\n')


def main():
    args = supply_args()
    my_fasta = pysam.FastaFile(filename=args.ref_fasta)
    tott = SampMetJson(args.samp_met).tott
    if args.manual_mode:
        my_coord = ManualRec([args.chr1, args.start1, args.end1, args.chr2, args.start2, args.end2])
        my_juncs = ChimJuncFile(args.raw_junc, my_coord).filt_cj([my_coord])
        Output(args.bedpe_out, my_juncs, None, my_fasta, tott).write_me()
    else:
        my_bed = BedPe(args.srch_junc)
        chim_junc_file = ChimJuncFile(args.raw_junc, my_bed)
        my_juncs = chim_junc_file.filt_cj(chim_junc_file.regions.regions)
        Output(args.bedpe_out, my_juncs, my_bed, my_fasta, tott).write_me()
    my_fasta.close()


if __name__ == "__main__":
    main()
