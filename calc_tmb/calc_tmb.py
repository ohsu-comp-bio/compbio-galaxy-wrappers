#!/usr/bin/env python

import argparse
import json
import vcf

VERSION = '0.1.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--infile', help='Input VCF')
    parser.add_argument('--outfile', help='Output JSON')
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
    tmb_cnt = 0

    for record in vcf_reader:
        if record.is_snp:
            entry = VcfRec(record)
            if entry.gnomad < 0.001:
                if "missense_variant" in entry.snpeff[0]:
                    if entry.tlod > 20.0:
                        tmb_cnt += 1
                        print(entry.rec)

    write_out(args.outfile, tmb_cnt)


if __name__ == "__main__":
    main()
