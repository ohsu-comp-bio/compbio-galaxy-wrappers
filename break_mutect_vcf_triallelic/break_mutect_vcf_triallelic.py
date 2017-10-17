#!/usr/bin/env python

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
    parser.add_argument('input_vcf', type=file, help='Input VCF file to process.')
    parser.add_argument('output_vcf', help='Output VCF with triallelic sites only.')
    parser.add_argument('output_vcf_tri', help='Output VCF without triallelic sites.')
    parser.add_argument('output_intervals', help='Output interval_list file to be processed through MuTect2.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def write_interval(locus):
    """
    From VCF entry, create an interval_list styled representation of the coordinate.
    chrom:start-stop
    Both start and stop are 1-based in this representation.
    """
    locus = locus.rstrip('\n').split('\t')
    chrom = locus[0]
    coord = locus[1]
    interval = chrom + ':' + coord + '-' + coord
    return interval


def main():

    args = supply_args()
    
    outfile = open(args.output_vcf, 'w')
    outfile_tri = open(args.output_vcf_tri, 'w')
    outfile_intervals = open(args.output_intervals, 'w')

    with args.input_vcf as vcf:
        for locus in vcf:
            if locus[0] == '#':
                outfile.write(locus)
                outfile_tri.write(locus)
            else:
                filter_col = locus.rstrip('\n').split('\t')[6]
                if 'triallelic_site' in filter_col:
                    outfile_tri.write(locus)
                    interval = write_interval(locus)
                    outfile_intervals.write(interval)
                    outfile_intervals.write('\n')
                else:
                    outfile.write(locus)

    outfile.close()
    outfile_tri.close()
    outfile_intervals.close()


if __name__ == "__main__":
    main()


