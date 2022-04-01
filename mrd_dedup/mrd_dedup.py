#!/usr/bin/env python

# DESCRIPTION: Dedup forced variant calls to be sent to CGD.  This will take away duplicates within the VCF you are trying to send,
# and will remove entries that have been seen in other VCF files imported in to CGD.
# USAGE:
# CODED BY: John Letaw

import argparse

VERSION = '0.1.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('mrd_vcf', help="MRD VCF to be sent to CGD")
    parser.add_argument('output_vcf', help="Output corrected VCF")
    parser.add_argument('input_vcfs', type=argparse.FileType('r'), nargs="+", help="Input VCFs to be compared against to look for dups")
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def all_vars(invcfs):
    """
    Create data structure to hold existing variants in input files (input_vcfs).
    :return:
    """
    all_vars = {}
    for my_vcf in invcfs:
        with my_vcf as vcf:
            for line in vcf:
                if not line.startswith('#'):
                    nline = line.rstrip('\n').split('\t')
                    chrom = nline[0]
                    pos = nline[1]
                    ref = nline[3]
                    alt = nline[4]
                    samp = nline[9]
                    uniq_key = (chrom, pos, ref, alt)
                    if uniq_key not in all_vars:
                        all_vars[uniq_key] = line
    return all_vars


def mrd_parse(invcf):
    """
    Get information we need from the MRD VCF.
    :return:
    """
    header = []
    var_dict = {}
    with open(invcf, 'rU') as myvcf:
        for line in myvcf:
            if not line.startswith('#'):
                nline = line.rstrip('\n').split('\t')
                chrom = nline[0]
                pos = nline[1]
                ref = nline[3]
                alt = nline[4]
                samp = nline[9]
                uniq_key = (chrom, pos, ref, alt)
                if uniq_key not in var_dict and samp != '.':
                    var_dict[uniq_key] = line
            else:
                header.append(line)

    return header, var_dict


def find_mrd_all(mrd, all_vars):
    """
    Check through the list of all variants to see if there are matches with MRD variants.
    :return:
    """
    to_del = []
    for key in mrd:
        if key in all_vars:
            to_del.append(key)
    for entry in to_del:
        mrd.pop(entry)

    return mrd


def main():

    args = supply_args()
    handle_out = open(args.output_vcf, 'w')
    header, mrd_dict = mrd_parse(args.mrd_vcf)
    variants = all_vars(args.input_vcfs)
    mrd_dict = find_mrd_all(mrd_dict, variants)

    for line in header:
        handle_out.write(line)
    for entry in mrd_dict.items():
        handle_out.write(entry)

    handle_out.close()

if __name__ == "__main__":
    main()
