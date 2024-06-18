#!/usr/bin/env python

"""
This tool attempts to determine whether the forward and reverse strand counts in a FreeBayes
VCF are biased, or not proportionate to each other.  Current implementation is utilizing Hoeffding's
Inequality to estimate upper and lower VAF bounds.  Will not assess 1/1 genotypes.

1.0.0 - Rewrite.
1.2.1 - Changed conf to 0.99999999999
1.3.0 - Don't apply StrandBias label when genotype is 1/1; make conf a parameter
1.3.1 - Removed vcftype and callers parameter and use VCF fields check to get strand counts.
1.3.2 - Include previous strand bias method.
1.3.2.4 - Exclude variants with alt strand counts below 10
"""

from __future__ import print_function
import argparse
import vcfpy
import numpy as np

from filter_strand_v0 import adj_alts


VERSION = '1.3.2.4'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('infile', help="Input VCF.")
    parser.add_argument('outfile', help="Output VCF.")
    parser.add_argument('--strandbias_flag', type=str, default='StrandBias', help='Strand bias flag.')
    parser.add_argument('--method', choices=['adjust_alts', 'hoeffding'], default='adjust_alts',
                        help="Method to use to access strand bias. Choices: adjust_alts, hoeffding.")
    parser.add_argument('--conf', type=float, default=0.9999, help="Confidence level for hoeffding calculation.")
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def get_strand_count(record):
    if 'SRF' in record.INFO:
        return [record.INFO['SRF'], record.INFO['SRR'], record.INFO['SAF'][0], record.INFO['SAR'][0]]
    elif 'SB' in record.calls[0].data and record.calls[0].data['SB']:
        return record.calls[0].data['SB']
    else:
        return None


class VcfReader:
    def __init__(self, filename):
        self.myvcf = vcfpy.Reader.from_path(filename)
        self.header = self.myvcf.header
        self.vrnts = self._read_vcf()
        self.myvcf.close()

    def _read_vcf(self):
        vrnts = []
        for vrnt in self.myvcf:
            vrnts.append(vrnt)
        return vrnts


class StrandOps:
    def __init__(self, info, conf):
        self.info = info
        self.conf = conf

    def assess_strand(self):
        """
        Check SAF, SAR, SRR, SRF fields to find out if we should call this StrandBias.
        :return:
        """
        if self.info:
            srf = self.info[0]
            srr = self.info[1]
            saf = self.info[2]
            sar = self.info[3]
            alt_dp = self._alt_dp(saf, sar)
            hoeff = self._hoeffding_t(alt_dp, self.conf)
            ref_vaf = self._strand_freq(srf, srr)
            alt_vaf = self._strand_freq(saf, sar)
            if ref_vaf:
                lower = max(ref_vaf - hoeff, 0)
                upper = min(ref_vaf + hoeff, 1)
            else:
                lower = 0
                upper = 1

            if lower <= alt_vaf <= upper:
                return True
            return False
        else:
            return False

    @staticmethod
    def _hoeffding_t(n, conf=0.95):
        """
        Find the predicted epsilon value based on alt read counts.
        This will need to be altered to allow for confidence intervals that aren't 0.95.
        :param n:
        :return:
        """
        # Strand count can sometimes be 0 in mutect2
        # https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected
        try:
            return np.sqrt((-1/(2*n))*(np.log((1-conf)/2)))
        except ZeroDivisionError:
            return 0

    @staticmethod
    def _set_epsilon(fwd, rev, fac=2):
        """
        Get the estimate of epsilon from the reference reads.
        :param fwd:
        :param rev:
        :return:
        """
        return (fwd / (fwd + rev)) / fac

    @staticmethod
    def _alt_dp(fwd, rev):
        """
        Get the total count of alternate allele reads.
        :return:
        """
        return fwd + rev

    @staticmethod
    def _strand_freq(fwd, rev):
        """
        Get the frequency of forward oriented reads.
        :return:
        """
        try:
            return fwd / (fwd + rev)
        except ZeroDivisionError:
            return 0


class FilterAdd:
    def __init__(self, filt):
        self.filt = filt

    def add_filt(self, text):
        """
        Add an annotation to the FILTER column.
        :return:
        """
        self.filt.append(text)


def main():
    args = supply_args()
    infile = VcfReader(args.infile)
    filter_annot = args.strandbias_flag
    # Add the header entry.
    header_add = vcfpy.OrderedDict()
    header_add['ID'] = filter_annot
    header_add['Description'] = f"Evidence of strand orientation bias at this locus according to bounds defined by {args.method} method."
    infile.header.add_filter_line(header_add)
    # Open writer using newly created header.
    writer = vcfpy.Writer.from_path(args.outfile, infile.header)

    for entry in infile.vrnts:
        strand_res = True
        info = get_strand_count(entry)
        if info:
            if entry.calls[0].gt_type != 2:
                if args.method == 'adjust_alts':
                    if info[2] >= 10 and info[3] >= 10 and info[0] > 0 and info[1] > 0:
                        strand_res = adj_alts(*info)
                elif args.method == 'hoeffding':
                    strand_res = StrandOps(info, args.conf).assess_strand()
                else:
                    raise Exception("Specify either adjust_alts or hoeffding methods")
        if not strand_res:
            filt = FilterAdd(entry.FILTER)
            filt.add_filt(text=filter_annot)
        writer.write_record(entry)

    writer.close()


if __name__ == "__main__":
    main()
