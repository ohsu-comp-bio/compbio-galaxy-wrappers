#!/usr/bin/env python

# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
from gatk_intervals import PerLocusRead
from gatk_intervals import ProbeQcRead
from vcf import VcfReader
from vcf import VcfWriter

import argparse
import sys

VERSION = '0.1.1'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--probeqc_before', type=ProbeQcRead,
                        help='ProbeQC before UMI deduplication')
    parser.add_argument('--probeqc_after', type=ProbeQcRead,
                        help='ProbeQC after UMI deduplication')
    parser.add_argument('--perlocus_before', type=PerLocusRead,
                        help='per locus before UMI deduplication')
    parser.add_argument('--perlocus_after', type=PerLocusRead,
                        help='per locus after UMI deduplication')
    parser.add_argument('invcf', type=VcfReader, help='Input VCF')
    parser.add_argument('outvcf', help='Output VCF')
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def get_avgd(record, umi_avgd):
    """
    Get the average depth from either the per locus or ProbeQC structures.
    :return:
    """
    for entry, avgd in umi_avgd.iteritems():
        if record[0] == entry[0] and record[1] == entry[1]:
            return '{:0.3f}'.format(avgd)

    return '.'


def main():

    args = supply_args()

    # For this tool, we will always set the header using below values.
    hfield = 'FORMAT'
    hid = 'PUMI'
    hnum = '1'
    htype = 'Float'
    hdesc = 'Percent UMI coverage for this position.'

    test_before = args.perlocus_before
    test_after = args.perlocus_after
    umi_avgd = test_after / test_before

    vcf_in = args.invcf
    vcf_in.add_header(hfield, hid, hnum, htype, hdesc)

    for record in vcf_in.vcf:
        chrom = record[0]
        coord = record[1]
        avgd = get_avgd(record, umi_avgd)
        vcf_in.add_format_field(chrom, coord, hid, avgd)

    VcfWriter(args.outvcf, vcf_in)


if __name__ == "__main__":
    main()
