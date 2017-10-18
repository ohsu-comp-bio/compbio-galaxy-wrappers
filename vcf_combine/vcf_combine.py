#!/usr/bin/env python

# DESCRIPTION:
# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
import argparse

VERSION = '0.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('remove_vcf', help='Remove variants from this VCF.')
    parser.add_argument('vcf_cands', help='Variants in this VCF will be '
                                          'removed from the other VCF if '
                                          'they exist.')
    parser.add_argument('output_vcf', help='Output VCF.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def build_cands(filename):
    """
    Create a list of variants that need to be removed from the other VCF.
    :return:
    """
    cands = []
    with open(filename, 'rU') as myfile:
        for line in myfile:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                variants = (line[0], line[1])
                cands.append(variants)
    return cands


def remove_cands(cands, infile, outfile):
    """
    Remove variants in cands from this file.
    :param filename:
    :return:
    """
    handle_out = open(outfile, 'w')
    with open(infile) as rem_file:
        for line in rem_file:
            if line.startswith('#'):
                handle_out.write(line)
            else:
                sline = line.rstrip('\n').split('\t')
                variant = (sline[0], sline[1])
                if variant not in cands:
                    handle_out.write(line)

    handle_out.close()


def main():

    args = supply_args()
    cands = build_cands(args.vcf_cands)
    remove_cands(cands, args.remove_vcf, args.output_vcf)


if __name__ == "__main__":
    main()
