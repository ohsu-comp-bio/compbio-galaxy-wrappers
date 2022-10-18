#!/usr/bin/env python

"""
Remove mrdf label from a merged VCF FILTER field of a record if variant is called by other callers.

Example usage: mrdf_label.py 'input.vcf' 'output.vcf'
"""

import argparse
import vcfpy

VERSION = '1.0.0'


def get_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', help="Input VCF")
    parser.add_argument('output_vcf', help="Output VCF")
    parser.add_argument('--label', help='Label to be removed from VCF FILTER field if conditions are met')
    parser.add_argument('--filter_condition', help='FILTER annotations corresponding to variant callers being used')
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)

    args = parser.parse_args()
    return args


def main():
    args = get_args()

    reader = vcfpy.Reader.from_path(args.input_vcf)
    header = reader.header
    writer = vcfpy.Writer.from_path(args.output_vcf, header=header)

    for record in reader:
        if args.label in record.FILTER and any(caller in record.FILTER for caller in args.callers.split(' ')):
            record.FILTER.remove(args.label)

        writer.write_record(record)


if __name__ == "__main__":
    main()
