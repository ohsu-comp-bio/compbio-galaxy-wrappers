#!/usr/bin/env python

# USAGE:
# CODED BY: John Letaw

import argparse
from natsort import natsorted

VERSION = '0.1.1'

def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcfs', type=file, nargs="+", help="Input VCF's to be merged based on the AD field.")
    parser.add_argument('output_vcf', help="Output VCF.")
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def compare_entries(var1, var2, field='AD'):
    """
    Compare two VCF entries and return the one that we want to keep.
    """
    ad1 = var1[9].split(':')[1]
    ad2 = var2[9].split(':')[1]
    if add_ad(ad1) >= add_ad(ad2):
        return var1
    return var2


def add_ad(ad_in):
    """
    Add values in AD field to get a total.
    values are comma-sep.
    """
    total = 0
    for entry in ad_in.split(','):
        total += int(entry)

    return total


def main():

    variants = {}
    args = supply_args()
    handle_out = open(args.output_vcf, 'w')

    for my_vcf in args.input_vcfs:
        with my_vcf as vcf:
            for variant in vcf:
                if variant[0] == '#' and args.input_vcfs.index(my_vcf) == 0:
                    handle_out.write(variant)
                elif variant[0] != '#':
                    variant_s = variant.rstrip('\n').split('\t')
                    chrom = variant_s[0]
                    coord = variant_s[1]
                    ref = variant_s[3]
                    alt = variant_s[4]
                    uniq_id = (chrom, coord)

                    if uniq_id not in variants:
                        variants[uniq_id] = variant_s
                    else:
                        curr_var = variants[uniq_id]
                        this_var = variant_s
                        variants[uniq_id] = compare_entries(curr_var, this_var)  

    for variant in natsorted(variants.itervalues(), key=lambda x: (x[0], x[1])):
        handle_out.write('\t'.join(variant))
        handle_out.write('\n')

if __name__ == "__main__":
    main()


