#!/usr/bin/env python

# DESCRIPTION:
# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
from file_types.bed import BedReader
import argparse

VERSION = '0.1.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('bed', help='Input BED')
    parser.add_argument('vcf', help='Input VCF')
    parser.add_argument('outfile', help='Output VCF')
    parser.add_argument('anno', help='Annotation to be added to FILTER '
                                       'column.')
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


def main():

    args = supply_args()
    handle_out = open(args.outfile, 'w')
    my_bed = BedReader(args.bed).split_coords()

    with open(args.vcf, 'rU') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                nline = line.rstrip('\n').split('\t')
                chrom = nline[0]
                pos = int(nline[1])
                filter = nline[6]

                if pos in my_bed[chrom]:
                    nline[6] = replace_filter(filter, args.anno)

                handle_out.write('\t'.join(nline))
                handle_out.write('\n')

            else:
                handle_out.write(line)

    handle_out.close()

if __name__ == "__main__":
    main()
