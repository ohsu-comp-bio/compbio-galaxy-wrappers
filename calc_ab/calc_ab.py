#!/usr/bin/env python

import argparse
import json
import numpy
import vcf

VERSION = '0.0.3'

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
        try:
            self.ad = rec.samples[0]['AD']
            self.af = self.ad[1] / (self.ad[0] + self.ad[1] + 0.0)
        except:
            self.ad = [0, 0]
            self.af = None
        self.total_ad = self.ad[0] + self.ad[1]
        self.gt = rec.samples[0]['GT']
        if len(rec.alleles[1:]) > 1:
            self.is_biallelic = False
        else:
            self.is_biallelic = True
        if self.is_biallelic and self.ad != [0, 0]:
            self.ab = self._calc_ab()
        else:
            self.ab = None

    def _calc_ab(self):
        """
        Calculate for the variant
        :return:
        """
        ratio = None
        # Need at least some number of reads to assess the site.
        if self.total_ad > 50:
            if self.gt == '0/1':
                # Hard to determine what this should be...
                if self.af > 0.1 and self.af < 0.9:
                    ratio = abs(0.5 - self.af)
        # if self.gt == '0/0':
        #     ratio = abs(af)
        return ratio


class WholeVcf(object):
    """

    """
    def __init__(self, filename):
        self.vcf_reader = vcf.Reader(filename=filename)
        # self.avg_ab = self._find_matches()
        self.stdev = self._find_stdev()

    def _find_stdev(self):
        """
        Find stdev of AF values.
        :return:
        """
        avg_afs = []
        for record in self.vcf_reader:
            entry = VcfRec(record)
            if entry.gt == '0/1' and entry.total_ad >= 50:
                af_diff = abs(0.5 - entry.af)
                if af_diff < 0.2:
                    avg_afs.append(af_diff)
            # if entry.gt == '0/0' and entry.total_ad >= 50:
            #     avg_afs.append(entry.af)

        num_arr = numpy.array(avg_afs)
        if num_arr != []:
            return numpy.std(num_arr, axis=0, ddof=1) * 100
        else:
            return 100.0

    def _find_matches(self):
        """
        Find the variants that match the criteria.
        :return:
        """
        total_ab = 0.0
        cnt = 0
        for record in self.vcf_reader:
            entry = VcfRec(record)
            if entry.ab is not None:
                total_ab += entry.ab
                cnt += 1

        if cnt != 0:
            val = (total_ab / cnt) * 100
        else:
            val = 100

        return val

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
    writer.write_json_out(my_vcf.stdev)

if __name__ == "__main__":
    main()
