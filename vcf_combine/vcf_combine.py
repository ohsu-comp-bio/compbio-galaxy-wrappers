#!/usr/bin/env python

# DESCRIPTION:
# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
from collections import OrderedDict
import argparse

VERSION = '0.2.0'


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
    parser.add_argument('output_vcf_base', help='Output VCF for source variants.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def build_cands(filename, outfile):
    """
    Create a list of variants that need to be removed from the other VCF.
    :return:
    """
    handle_out = open(outfile, 'w')
    cands = OrderedDict()
    with open(filename, 'rU') as myfile:
        for line in myfile:
            if not line.startswith('#'):
                sline = line.rstrip('\n').split('\t')
                variants = (sline[0], sline[1])
                cands[variants] = line
            else:
                handle_out.write(line)

    return cands, handle_out


def remove_cands(cands, infile, outfile, outfile_base):
    """
    Remove variants in cands from this file.
    :param filename:
    :return:
    """
    handle_out = open(outfile, 'w')
    annots = {}
    with open(infile) as rem_file:
        for line in rem_file:
            if line.startswith('#'):
                handle_out.write(line)
            else:
                sline = line.rstrip('\n').split('\t')
                filt = sline[6]
                variant = (sline[0], sline[1])
                if variant not in cands:
                    handle_out.write(line)
                else:
                    annots[variant] = filt

    handle_out.close()
    return annots


def add_filt(filename, annots, cands, handle_out):
    for variant in cands:
        if variant not in annots:
            handle_out.write(cands[variant])
        else:
            nline = cands[variant].rstrip('\n').split('\t')
            new_filt = ';'.join([nline[6], annots[variant]])
            to_write = nline[:6]
            to_write.append(new_filt)
            to_write.extend(nline[7:])
            handle_out.write('\t'.join(to_write))
            handle_out.write('\n')
    handle_out.close()


def main():

    args = supply_args()
    cands, handle_out = build_cands(args.vcf_cands, args.output_vcf_base)
    annots = remove_cands(cands, args.remove_vcf, args.output_vcf, args.output_vcf_base)
    add_filt(args.output_vcf_base, annots, cands, handle_out)

if __name__ == "__main__":
    main()
