#!/usr/bin/env python

"""
Hotspots inputs

Usage:
    hotspots_inputs.py ordered_test_hotspots.txt "Thrombocytosis Panel" --output_bed hotspots.bed --output_vcf hotspots.vcf

Details:
Hotspots inputs
"""

import argparse
from operator import attrgetter
from natsort import natsorted

VERSION = '1.0.0'


class Roi:
    def __init__(self, line, header):
        self.chrom = line[header.index('chrom')]
        self.start = int(line[header.index('pos')]) - 1
        self.end = line[header.index('pos')]
        self.pos = line[header.index('pos')]
        self.id = '.'
        self.ref = line[header.index('ref')]
        self.alt = line[header.index('alt')]
        self.qual = '.'
        self.filter = '.'
        self.info = '.'
        self.test = line[header.index('ordered test')]
        self.var_id = '{}:{}{}>{}'.format(self.chrom, self.pos, self.ref, self.alt)


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('input_file', help="Input file")
    parser.add_argument('ordered_test', help="Ordered test")
    parser.add_argument('--output_bed', help="Output BED.")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    roi_dict = {}
    with open(args.input_file, 'r') as infile:
        header = [col.lower() for col in next(infile).rstrip('\n').split('\t')]
        for line in infile:
            line = line.rstrip('\n').split('\t')
            roi = Roi(line, header)
            if roi.test == args.ordered_test:
                roi_dict[roi.var_id] = roi

    if args.output_bed:
        with open(args.output_bed, 'w') as bed_writer:
            if len(roi_dict) > 0:
                for roi in natsorted(roi_dict.values(), key=attrgetter('chrom', 'pos')):
                    bed_writer.write('{}\t{}\t{}\n'.format(roi.chrom, roi.start, roi.end))
            else:  # add fake pos if hotspots not in ordered tests hotspots file
                bed_writer.write('{}\t{}\t{}\n'.format('1', '1', '2'))

    if args.output_vcf:
        vcf_header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
        with open(args.output_vcf, 'w') as vcf_writer:
            vcf_writer.write('##fileformat=VCFv4.2\n')
            vcf_writer.write('\t'.join(vcf_header)+'\n')
            if len(roi_dict) > 0:
                for roi in natsorted(roi_dict.values(), key=attrgetter('chrom', 'pos')):
                    vcf_writer.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(roi.chrom, roi.pos, roi.id, roi.ref, roi.alt, roi.qual, roi.filter, roi.info))
            else:  # add fake pos if hotspots not in ordered tests hotspots file
                vcf_writer.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('1', '3', '.', 'T', 'C', '.', '.', '.'))

if __name__ == "__main__":
    main()
