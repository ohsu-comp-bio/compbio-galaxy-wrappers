#!/usr/bin/env python

"""
This tool attempts to determine whether the forward and reverse strand counts in a FreeBayes
VCF are biased, or not proportionate to each other.  Current implementation is utilizing Hoeffding's
Inequality to estimate upper and lower VAF bounds.

1.0.0 - Rewrite.
"""

from __future__ import print_function
import argparse
import vcfpy
import numpy as np

VERSION = '1.2.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--infile', help='Input VCF.')
    parser.add_argument('--vcftype', help='Type of VCF.')
    parser.add_argument('--callers', help='Callers in a merged VCF.')
    parser.add_argument('--outfile', help='Output VCF.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def m2_strand_count(record):
    srf, srr, saf, sar = record.calls[0].data['SB']
    return [srf, srr, saf, sar]


def fb_strand_count(record):
    srf, srr, saf, sar, = record.INFO['SRF'], record.INFO['SRR'], record.INFO['SAF'][0], record.INFO['SAR'][0]
    return [srf, srr, saf, sar]


def get_strand_count(record, vcftype, callers):
    if vcftype == 'merged':
        callers = callers.split(' ')
        for flt in record.FILTER:
            if flt in callers:
                vcftype = flt[:2]
                break
    if vcftype == 'm2':
        return m2_strand_count(record)
    elif vcftype == 'fb':
        return fb_strand_count(record)
    else:
        pass
        # raise Exception('Unrecognized VCF type at {}:{}'.format(record.CHROM, record.POS))


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
    def __init__(self, info):
        self.info = info

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
            hoeff = self._hoeffding_t(alt_dp, conf=0.99)
            ref_vaf = self._strand_freq(srf, srr)
            alt_vaf = self._strand_freq(saf, sar)
            if ref_vaf:
                lower = max(ref_vaf - hoeff, 0)
                upper = min(ref_vaf + hoeff, 1)
            else:
                lower = 0
                upper = 1

            # if saf != 0 and sar != 0:
            #     print("SAF: {0}".format(saf))
            #     print("SAR: {0}".format(sar))
            #     print("SRF: {0}".format(srf))
            #     print("SRR: {0}".format(srr))
            #     print("Lower: {0}".format(lower))
            #     print("Upper: {0}".format(upper))
            #     print("Alt VAF: {0}".format(alt_vaf))

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
    filter_annot = 'StrandBias'
    conf = 0.95
    # Add the header entry.
    header_add = vcfpy.OrderedDict()
    header_add['ID'] = filter_annot
    header_add['Description'] = 'Evidence of strand orientation bias at this locus according to bounds defined ' \
                                'by Hoeffdings Inequality.'
    infile.header.add_filter_line(header_add)
    # Open writer using newly created header.
    writer = vcfpy.Writer.from_path(args.outfile, infile.header)

    for entry in infile.vrnts:
        info = get_strand_count(entry, args.vcftype, args.callers)
        if info:
            strand = StrandOps(info)
            strand_res = strand.assess_strand()
            if not strand_res:
                #print(entry.CHROM, entry.POS, entry.FILTER)
                filt = FilterAdd(entry.FILTER)
                filt.add_filt(text=filter_annot)
        writer.write_record(entry)

    writer.close()


if __name__ == "__main__":
    main()
