#!/usr/bin/env python

import argparse
import json
import vcf

VERSION = '0.1.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--infile', help='Input VCF')
    parser.add_argument('--outfile', help='Output JSON')
    parser.add_argument('--gnomad_af', default=0.00001, type=float, help='GNOMAD VAF to filter list upon.  Values greater than this will not be utilized.')
    parser.add_argument('--m2_tlod', default=20.0, type=float, help='M2 TLOD score to filter list upon.  Values less than this will not be utilized.')
    parser.add_argument('--seq_space', default=0.61, type=float, help='Total genomic size of sequence targets, in Mb.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

class VcfRec(object):
    def __init__(self, rec):

        self.rec = rec
        self.chrom = str(rec.CHROM)
        self.coord = str(rec.POS)

        try:
            self.gnomad = float(rec.INFO['gnomad.AF'][0])
        except:
            self.gnomad = 0.0

        try:
            self.snpeff = rec.INFO['ANN']
        except:
            self.snpeff = None

        try:
            self.tlod = float(rec.INFO['TLOD'][0])
        except:
            self.tlod = None


def write_out(filename, val):
    """
    Prepare output json file.
    :return:
    """
    outfile = open(filename, 'w')
    out_metric = {'tmb': val}
    json.dump(out_metric, outfile)
    outfile.close()

def main():
    args = supply_args()
    vcf_reader = vcf.Reader(open(args.infile, 'r'))
    tmb_cnt = 0.0

    for record in vcf_reader:
        if record.is_snp:
            entry = VcfRec(record)
            if entry.gnomad < args.gnomad_af:
                if "missense_variant" in entry.snpeff[0]:
                    if entry.tlod > args.m2_tlod:
                        tmb_cnt += 1

    tmb_cnt = '{:.1f}'.format(tmb_cnt / args.seq_space)

    write_out(args.outfile, tmb_cnt)


if __name__ == "__main__":
    main()
