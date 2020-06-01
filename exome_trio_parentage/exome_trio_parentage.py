#!/usr/bin/env python3

# Ran the following command to prepare the gnomad VCF:
# gatk SelectVariants -R $HG19 -V ../../resources/gnomad/gnomad.exomes.r2.1.1.sites.vcf.bgz
# -L agilent_cre.intervals --select-type-to-include SNP --restrict-alleles-to BIALLELIC
# -O gnomad_agilent_snp_biallelic_10pct.vcf.gz --selectExpressions 'QUAL>1000000'
# --selectExpressions 'AF>0.05' --selectExpressions 'AF<0.95' --select-random-fraction 0.1

import argparse
import json
from scipy.stats import binom_test
import vcf

VERSION = '0.1.1'


def supply_args():
    parser = argparse.ArgumentParser(description='Determine parentage for exome trios, give a three-column VCF.')
    parser.add_argument('--input_vcf', help='Input three-column VCF to peform analysis on.')
    parser.add_argument('--output_vcf', help='Output three-column VCF containing discordant results.')
    parser.add_argument('--output_json', help='Output JSON with relevant parentage metrics.')
    parser.add_argument('--mother', help='VCF ID for mother.')
    parser.add_argument('--father', help='VCF ID for father.')
    parser.add_argument('--proband', help='VCF ID for proband.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION + '\n')
    args = parser.parse_args()
    return args


class TrioVcf:
    """
    Get information from the three-column vcf and organize it for further parentage calculations.
    Columns to be arranged in order mother, father, proband.
    """
    def __init__(self, filename, mother, father, proband, outfile):
        vcf_reader = vcf.Reader(open(filename, 'r'))
        self.mother = mother
        self.father = father
        self.proband = proband
        vcf_writer = vcf.Writer(open(outfile, 'w'), vcf_reader)

        self.total = 0
        self.disc = 0

        for record in vcf_reader:
            this_rec = VcfRec(record, self.mother, self.father, self.proband)
            if this_rec.valid:
                self.total += 1
                if not this_rec.geno_disc_valid:
                    self.disc += 1
                    vcf_writer.write_record(record)

        vcf_writer.close()

    def _binom_test(self, pval=0.001):
        """
        Calculate the binomial test.
        :return:
        """
        return binom_test(self.disc, self.total, pval)

    def write_json_out(self, filename):
        """
        Prepare output json file.
        :return:
        """
        outfile = open(filename, 'w')
        binom = self._binom_test()
        out_metric = {'parentage_sites': self.total,
                      'parentage_disc': self.disc,
                      'parentage_binom': binom,
                      'parentage_confirmed': self._parentage_confirm(binom)}
        json.dump(out_metric, outfile)
        outfile.close()

    def _parentage_confirm(self, binom, cutoff=0.0001, total_cutoff=200):
        """
        Decide whether the parentage is confirmed or not, and return a string that will be passed
        to sample metrics via json.
        0 means fail
        1 means pass
        999 means inconclusive due to total_cutoff
        :param binom:
        :return:
        """
        if self.total < total_cutoff:
            return "999"
        elif binom <= cutoff:
            return "0"
        elif binom > cutoff:
            return "1"
        else:
            raise Exception("parentage_confirm error, please check")


class VcfRec:
    """

    """
    def __init__(self, rec, mother, father, proband, het=0.05, hom=0.02):
        self.rec = rec
        self.het = het
        self.hom = hom
        self.mother = mother
        self.father = father
        self.proband = proband

        self.trio_valid = self._trio_valid()
        self.depth_valid = self._depth_check()
        self.ab_valid = self._allele_balance_check()
        self.geno_valid = self._geno_check()
        self.biallelic_valid = self._biallelic_check()
        self.base_valid = self._base_check()
        self.gq_valid = self._gq_check()
        self.ad_valid = self._ad_check()
        self.chrom_valid = self._chrom_check()
        self.checks = [self.trio_valid, self.depth_valid, self.ab_valid, self.geno_valid,
                       self.biallelic_valid, self.base_valid, self.gq_valid, self.ad_valid,
                       self.chrom_valid]

        if False in self.checks:
            self.valid = False
        else:
            self.valid = True

        self.geno_disc_valid = self._geno_disc_check()

    def _chrom_check(self):
        """
        Check to make sure there are not X or Y chroms.
        :return:
        """
        if self.rec.CHROM == 'X' or self.rec.CHROM == 'Y':
            return False
        return True

    def _ad_check(self):
        """
        Check the values in AD to see if they are consistent with what GATK is deciding in PL/GT.
        No reason to use these inconsistent records, we have plenty others.
        :return:
        """
        for samp in self.rec.samples:
            samp_ab = self._calc_ab(samp['AD'])
            if samp['GT'] == '0/0':
                if samp_ab > self.hom:
                    return False
            if samp['GT'] == '1/1':
                if samp_ab < (1 - self.hom):
                    return False
            if samp['GT'] == '0/1':
                if samp_ab < (.5 - self.het) or samp_ab > (.5 + self.het):
                    return False
        return True

    def _gq_check(self, gq=99):
        """
        If the value of GQ in the sample field of the VCF is not 99, do not use.
        :return:
        """
        for samp in self.rec.samples:
            if samp['GQ']:
                if int(samp['GQ']) < gq:
                    return False
            else:
                return False
        return True

    def _base_check(self):
        """
        Ensure the alleles are only made up of [ATCG].
        :return:
        """
        bases = ['A', 'T', 'C', 'G']
        if len(self.rec.ALT) == 1:
            if self.rec.REF in bases and self.rec.ALT[0] in bases:
                return True
        return False

    def _biallelic_check(self):
        """
        Ensure there is only one alternate allele, and there that the change is a SNP.
        :return:
        """
        if len(self.rec.REF) != 1:
            return False
        if len(self.rec.ALT) != 1:
            return False
        else:
            if len(self.rec.ALT[0]) != 1:
                return False
        return True

    def _geno_check(self):
        """
        Check to make sure only genotypes containing 1's or 0's exist.
        :return:
        """
        valid = ['0', '1']
        for samp in self.rec.samples:
            for allele in samp['GT'].split('/'):
                if allele not in valid:
                    return False
        return True

    @staticmethod
    def _split_gt(gt):
        """
        Split the GT field up by either / or |.
        :return:
        """
        if '/' in gt:
            return gt.split('/')
        elif '|' in gt:
            return gt.split('|')
        else:
            return gt

    def _geno_disc_check(self):
        """
        PyVCF provides samples ordered.
        :return:
        """
        mother = self._split_gt(self.rec.genotype(self.mother)['GT'])
        father = self._split_gt(self.rec.genotype(self.father)['GT'])
        proband = self._split_gt(self.rec.genotype(self.proband)['GT'])
        pro1 = proband[0]
        pro2 = proband[1]

        if pro1 in mother and pro2 in father:
            return True
        if pro2 in mother and pro1 in father:
            return True
        return False

    def _depth_check(self, depth=50):
        """
        If any of the samples contained in this entry are < $depth, result is False.
        :param depth:
        :return:
        """
        for samp in self.rec.samples:
            if samp['DP']:
                if int(samp['DP']) < depth:
                    return False
            else:
                return False
        return True

    @staticmethod
    def _calc_ab(ad):
        """
        Calculate the allele balance for this particular sample.
        :param ad:
        :return:
        """
        ref = float(ad[0]) + 0.0
        alt = float(ad[1]) + 0.0

        try:
            return alt / (alt + ref)
        except ZeroDivisionError:
            return 0

    def _allele_balance_check(self):
        """
        Everything should be less within +/- het of 0.5 and hom of 1 or 0.
        :return:
        """
        for samp in self.rec.samples:
            ab = self._calc_ab(samp['AD'])
            lower_het = 0.5 - self.het
            upper_het = 0.5 + self.het
            lower_hom = self.hom
            upper_hom = 1.0 - self.hom

            if lower_hom < ab < lower_het:
                return False
            if upper_hom > ab > upper_het:
                return False

        return True

    def _trio_valid(self, trio=3):
        """
        If the number of samples in this record are < $trio, result is False.
        :param trio:
        :return:
        """
        if len(self.rec.samples) != trio:
            return False
        return True


def main():

    args = supply_args()
    myvcf = TrioVcf(args.input_vcf, args.mother, args.father, args.proband, args.output_vcf)
    myvcf.write_json_out(args.output_json)


if __name__ == "__main__":
    main()
