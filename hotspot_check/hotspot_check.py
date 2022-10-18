#!/usr/bin/env python

"""
Check if hotspots are annotated as above or below background.

Example usage: hotspot_check.py 'input.vcf' 'output.vcf'

Output: JSON file of above and below background hotspot counts.
"""

import argparse
import vcfpy

VERSION = '1.0.0'


def get_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', help="Input VCF")
    parser.add_argument('output_metrics', help="Output JSON file containing hotspot metrics")
    parser.add_argument('--label', help='Hotspot label')
    parser.add_argument('--filter_condition', help='FILTER annotation to use for filtering')
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)

    args = parser.parse_args()
    return args


def main():
    args = get_args()

    reader = vcfpy.Reader.from_path(args.input_vcf)

    above = []
    below = []
    for record in reader:

        if args.label in record.FILTER:
            if any(filt in record.FILTER for filt in args.filter_condition.split(' ')):
                below.append(record)
            else:
                above.append(record)

    with open(args.output_metrics, 'w') as output_metrics:
        output_metrics.write("{{\"forced_calls_above\": {}, \"forced_calls_below\": {}}}".format(len(above), len(below)))


if __name__ == "__main__":
    main()
