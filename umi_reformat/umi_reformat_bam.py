#!/usr/bin/env python

# USAGE: python umi_reformat_bam.py
# CODED BY: John Letaw
# TODO: Add @PG header entry for this tool.

import argparse
import pysam

VERSION = '0.3.0'

def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='Include the BX tag for UMIs in the proper position of the BAM file.')
    parser.add_argument('input_bam', help='Input BAM File')
    parser.add_argument('output_bam', help='Output BAM File')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def main():

    args = supply_args()

    mybam = pysam.AlignmentFile(args.input_bam, "rb")
    outbam = pysam.AlignmentFile(args.output_bam, "wb", template=mybam)

    for align in mybam:
        ident = align.query_name
        umi = ident.split(':')[-1]
        #        this_seq = align.query_sequence
#        umi = this_seq[:12]
#        align.query_sequence = this_seq[22:]
        align.set_tag('RX', umi)
        outbam.write(align)

        # Reorganize this.  This is for the situation where we are parsing bcl2fastq UMI trimmings.
        # ident = align.query_name
        # umi = ident.split(':')[-1]
        # align.query_name = ':'.join(ident.split(':')[:-1])

    outbam.close()
    mybam.close()


if __name__ == "__main__":
    main()


