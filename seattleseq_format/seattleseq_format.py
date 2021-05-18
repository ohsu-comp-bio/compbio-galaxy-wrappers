#!/usr/bin/env python

import argparse
import vcfpy

VERSION = '0.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('infile', help='Input VCF to pre-format.')
    parser.add_argument('outfile', help='Output VCF.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def main():
    
    args = supply_args()
    vcf = vcfpy.Reader.from_path(args.infile)
    writer = vcfpy.Writer.from_path(args.outfile, vcf.header)

    for vrnt in vcf:
        # Get rid of phased genotype characters.
        if vrnt.calls[0].gt_phase_char == '|':
            for entry in vrnt.calls:
                new_geno = entry.data['GT'].replace('|', '/')
                entry.data['GT'] = new_geno
        writer.write_record(vrnt)

    vcf.close()
    writer.close()


if __name__ == "__main__":
    main()
