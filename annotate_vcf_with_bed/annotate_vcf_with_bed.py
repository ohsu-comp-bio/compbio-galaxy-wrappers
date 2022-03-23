#!/usr/bin/env python

# DESCRIPTION: Based on a BED file, annotate VCFs in the FILTER column with a phrase of choice.

from bed import BedReader
from collections import OrderedDict
import argparse
import vcfpy


VERSION = '0.3.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('bed', help='Input BED')
    parser.add_argument('vcf', help='Input VCF')
    parser.add_argument('outfile', help='Output VCF')
    parser.add_argument('anno', help='Annotation to be added to FILTER column.')
    parser.add_argument('desc', help='Description of the annotation to be added to FILTER header entry.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def main():

    args = supply_args()
    vcf = vcfpy.Reader.from_path(args.vcf)
    # Add the header entry for the new FILTER.
    header_line = OrderedDict({"ID": args.anno,
                               "Description": args.desc})
    vcf.header.add_filter_line(header_line)
    writer = vcfpy.Writer.from_path(args.outfile, vcf.header)
    my_bed = BedReader(args.bed).split_coords()

    for rec in vcf:
        chrom = rec.CHROM
        pos = rec.POS
        if chrom in my_bed:
            if pos in my_bed[chrom]:
                rec.add_filter(args.anno)
        writer.write_record(rec)

    vcf.close()
    writer.close()


if __name__ == "__main__":
    main()
