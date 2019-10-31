#!/usr/bin/env python

import argparse
import pysam
import vcf

VERSION = '0.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--bam', help='BAM file to find sample name from read group in header.')
    parser.add_argument('--vcf', help='VCF that should have sample column replaced.')
    parser.add_argument('--vcf_out', help='VCF with the sample column replaced.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

def main():
    args = supply_args()
    vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    bam = pysam.AlignmentFile(args.bam, 'r')
    new_sample_id = bam.header['RG'][0]['SM']
    # Will hard code this as 'sample' for the time being, so that this tool is only used for the scalpel vcfs.
    # After I have a chance to write this a little differently, will expand use cases.
    if vcf_reader.samples[0] == 'sample':
        vcf_reader.samples = [new_sample_id]

    vcf_writer = vcf.Writer(open(args.vcf_out, 'w'), vcf_reader)
    for record in vcf_reader:
        vcf_writer.write_record(record)

    vcf_writer.close()

if __name__ == "__main__":
    main()
