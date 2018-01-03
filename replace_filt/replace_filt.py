#!/usr/bin/env python

# DESCRIPTION: Replace the FILTER field in a VCF, for all entries, or add something to it.
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
    parser.add_argument('input_vcf', help='Input VCF')
    parser.add_argument('output_vcf', help='Output VCF')
    parser.add_argument('anno', help='FILTER string to apply')
    parser.add_argument('size', type=int, help='Limit FILTER string to this size')

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def replace_filter(filt, anno):
    """
    Remove the dot if it's there, otherwise just append something with a
    semicolon.
    :return:
    """
    if filt == '.' or filt == 'PASS':
        return anno
    else:
        return ';'.join([filt, anno])


def check_filt_str(phrase, size=16):
    """
    Make sure the string isn't too long, and is only composed of alnum characters.
    :return:
    """
    if len(phrase) <= size and phrase.isalnum():
        return True
    return False


def create_new_vcf_line(line, annot):
    """
    Write a new line, containing the added FILTER value.
    :return:
    """
    nline = line[:6]
    nline.append(annot)
    nline.extend(line[7:])
    return '\t'.join(nline) +'\n'


def main():

    args = supply_args()
    handle_out = open(args.output_vcf, 'w')

    with open(args.input_vcf, 'rU') as invcf:
        for line in invcf:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                filt = line[6]
                annot = replace_filter(filt, args.anno)
                if check_filt_str(annot, args.size):
                    handle_out.write(create_new_vcf_line(line, annot))
                else:
                    raise ValueError("The FILTER string you are trying to enter is either greater than " +
                                     str(args.size) +  ", or contains non alphanumeric characters.")
            else:
                handle_out.write(line)
    handle_out.close()


if __name__ == "__main__":
    main()
