#!/usr/bin/env python

import argparse
import json
import vcf

VERSION = '0.0.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--infile',  help='Input VCF')
    parser.add_argument('--outfile', help='Output JSON')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

class VcfRec(object):
    def __init__(self, rec):

        self.rec = rec
        self.chrom = str(rec.CHROM)
        self.coord = str(rec.POS)
        self.ad = rec.samples[0]['AD']
        self.gt = rec.samples[0]['GT']
        if len(rec.alleles[1:]) > 1:
            self.is_biallelic = False
        else:
            self.is_biallelic = True
        if self.is_biallelic:
            self.ab = self._calc_ab()
        else:
            self.ab = None

    def _calc_ab(self):
        """
        Calculate for the variant
        :return:
        """
        af = self.ad[1] / (self.ad[0] + self.ad[1] + 0.0)
        ratio = None
        if self.gt == '0/1':
            if af > 0.35 and af < 0.65:
                ratio = abs(0.5 - af)
        elif self.gt == '0/0':
            ratio = abs(af)
        return ratio


class WholeVcf(object):
    """

    """
    def __init__(self, filename):
        self.vcf_reader = vcf.Reader(filename=filename)
        self.avg_ab = self._find_matches()
        print(self.avg_ab)

    def _find_matches(self):
        """
        Find the variants that match the criteria.
        :return:
        """
        total_ab = 0.0
        cnt = 0
        for record in self.vcf_reader:
            entry = VcfRec(record)
            if entry.ab:
                total_ab += entry.ab
                cnt += 1
        return total_ab / cnt


class Writer(object):
    """

    """
    def __init__(self, vcf_reader, outfile):
        self.outfile = outfile
        self.vcf_reader = vcf_reader

    def write_json_out(self, metric):
        """
        Prepare output json file.
        :return:
        """
        outfile = open(self.outfile, 'w')
        out_metric = {'allele_balance': metric}
        json.dump(out_metric, outfile)
        outfile.close()


def main():
    args = supply_args()
    my_vcf = WholeVcf(args.infile)
    writer = Writer(my_vcf.vcf_reader, args.outfile)
    writer.write_json_out(my_vcf.avg_ab)

if __name__ == "__main__":
    main()
