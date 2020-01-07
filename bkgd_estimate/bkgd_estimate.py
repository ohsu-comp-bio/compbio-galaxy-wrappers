#!/usr/bin/env python

"""
According to binomial test, determine whether variant VAF is significantly different than estimated background.
If it is not, we will apply a tag to the VCF FILTER column.
"""

from argparse import ArgumentParser
from scipy import stats
import vcf

VERSION = '0.0.1'

def supply_args():
    parser = ArgumentParser(description='')
    ### Required
    parser.add_argument("--vcf", help="Input VCF to apply background estimate VCF FILTERs to.")
    parser.add_argument("--outfile", help="Output VCF.")
    parser.add_argument("--bkgd", default=0.02, type=float, help="Estimate of assay background VAF.")
    parser.add_argument("--pval", default=0.001, type=float, help="P-value at which VAF will be deemed indeterminate.")
    args = parser.parse_args()
    return args

class BkgdEst(object):

    def __init__(self, bkgd, pval):
        self.bkgd = bkgd
        self.pval = pval

    def find_binom_pval(self, success, tries):
        """

        :return:
        """
        return stats.binom_test(success, tries, self.bkgd, alternative='two-sided')

    def check_val(self, val):
        """
        FALSE is value is less that self.pval.
        TRUE otherwise.
        :return:
        """
        if float(val) <= self.pval:
            return False
        else:
            return True

def main():
    args = supply_args()
    myvcf = vcf.Reader(open(args.vcf, 'r'))
    bkgd_tester = BkgdEst(args.bkgd, args.pval)
    desc_str = 'Variant AF is within designated background range of %g%% (p<%g).' %(100*args.bkgd, args.pval)
    myvcf.filters['IN_BKGD'] = vcf.parser._Filter('IN_BKGD', desc_str)
    vcf_writer = vcf.Writer(open(args.outfile, 'w'), myvcf)
    for entry in myvcf:
        depth = entry.samples[0]['AD'][0] + entry.samples[0]['AD'][1]
        alt = entry.samples[0]['AD'][1]
        pval = bkgd_tester.find_binom_pval(alt, depth)
        if not entry.FILTER:
            if bkgd_tester.check_val(pval) or depth < 10:
                entry.FILTER = ['IN_BKGD']
        else:
            if bkgd_tester.check_val(pval) or depth < 10:
                entry.FILTER.append('IN_BKGD')
        vcf_writer.write_record(entry)
    vcf_writer.close()

if __name__ == "__main__":
    main()
