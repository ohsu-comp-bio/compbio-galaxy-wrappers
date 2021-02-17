#!/usr/bin/env python

# DESCRIPTION: VCF files list multiple alternate alleles in the ALT columns,
#  separated by commas. This script will split entries with multiple alleles
# listed, and place each on a separate line.  This will allow us to import
# this data in to annotation software and properly receive prediction
# scores. Post-normalization (vt, bcftools) of this data is highly recommended.
# USAGE: split_mult_alleles.py -h

from copy import deepcopy
from file_types import vcfreader
from file_types import vcfwriter
import argparse
import sys

VERSION = '0.5.1'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(dest='input', help='')
    parser.add_argument(dest='output', help='')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class VcfRecDecomp(object):
    """
    0,1,2,... - No split
    A - alt1 = val1, alt2 = val2, etc.
    R - ref,alt1 = val1, ref,alt2 = val2, etc.
    G - 0/0,0/1,1/1;0/0,0/2,2/2;0/0,0/3,3/3
    . - no split

    Ordering of genotype fields, from spec
    F(j/k) = (k*(k+1)/2)+j
    AA,AB,BB,AC,BC,CC,AD,BD,CD,DD,AE,BE,CE,DE,EE

    idx1 - 0,1,2
    idx2 - 0,3,5
    idx3 - 0,6,9
    idx4 - 0,10,14

    i=0
    idx+idx*0
    """

    def __init__(self, vrnt, info_number, samples_number):
        self.info_number = info_number
        self.samples_number = samples_number
        self.vrnt_samples = self._decomp_samples_vrnt(vrnt)
        self.vrnt_info = self._decomp_info_vrnt(vrnt)
        self.decomp_vrnts = self._vrnt_update(vrnt)

    def _vrnt_update(self, vrnt):
        """
        Update VcfRecBase obj with split metrics.
        :return:
        """
        vrnt_updt = []
        print(self.vrnt_info)
        for allele, vals in self.vrnt_info.items():
            new_vrnt = deepcopy(vrnt)
            new_vrnt.ALT = [allele]
            for val in vals:
                new_vrnt.INFO[val] = vals[val]

            ind = 0
            for samp in self.vrnt_samples:
                for frmt in samp[allele]:
                    new_vrnt.samples[ind][frmt] = samp[allele][frmt]
                ind += 1
            vrnt_updt.append(new_vrnt)

        return vrnt_updt

    def _decomp_samples_vrnt(self, vrnt):
        """
        Decompose variant in to pieces, to be split.
        :return:
        """
        samp_changes = []
        for samp in vrnt.samples:
            split = {}
            for entry in samp:
                for alt in vrnt.ALT:
                    idx = vrnt.ALT.index(alt)
                    if alt not in split:
                        split[alt] = {}
                    if self.samples_number[entry] == 'A':
                        split[alt][entry] = self._decomp_alt(samp[entry], idx)
                    elif self.samples_number[entry] == 'R':
                        split[alt][entry] = self._decomp_ref(samp[entry], idx)
                    elif self.samples_number[entry] == 'G':
                        split[alt][entry] = self._decomp_geno(samp[entry], idx + 1)
                    elif entry == 'GT':
                        split[alt][entry] = self._decomp_gt(samp[entry])
            samp_changes.append(split)
        return samp_changes


    def _decomp_info_vrnt(self, vrnt):
        """
        Decompose variant in to pieces, to be split.
        :return:
        """
        split = {}
        for entry in vrnt.INFO:
            for alt in vrnt.ALT:
                idx = vrnt.ALT.index(alt)
                if alt not in split:
                    split[alt] = {}

                if self.info_number[entry] == 'A':
                    split[alt][entry] = self._decomp_alt(vrnt.INFO[entry], idx)
                elif self.info_number[entry] == 'R':
                    split[alt][entry] = self._decomp_ref(vrnt.INFO[entry], idx)
                elif self.info_number[entry] == 'G':
                    split[alt][entry] = self._decomp_geno(vrnt.INFO[entry], idx + 1)
                elif entry == 'GT':
                    split[alt][entry] = self._decomp_gt(vrnt.INFO[entry])

        return split


    def _decomp_gt(self, val):
        """
        Split the GT field apart.
        :param val:
        :return:
        """
        if len(val.split('/')) > 2:
            return "0/1"
        else:
            return "./."

    def _decomp_alt(self, val, idx):
        """
        Find fields in INFO that need to be split.
        Input is comma-deliminted value set, can be treated as a string.
        :return:
        """
        return val.split(',')[idx]

    def _decomp_ref(self, val, idx):
        """
        Find fields in SAMPLE columns that need to be split.
        :return:
        """
        new_val = ['0', val.split(',')[idx+1]]
        return ','.join(new_val)

    def _decomp_geno(self, val, idx):
        """

        :param val:
        :return:
        """
        geno_idxs = [self._calc_gl_pos(0, idx),
                     self._calc_gl_pos(idx, idx)]
        sval = val.split(',')
        new_val = [sval[0], sval[geno_idxs[0]], sval[geno_idxs[1]]]
        return ','.join(new_val)

    @staticmethod
    def _calc_gl_pos(ref, alt):
        """
        F(j/k) = (k*(k+1)/2)+j
        j = allele 1
        k = allele 2
        :return:
        """
        return int(((ref*(ref+1)/2)+alt))


def main():
    """
    # 1
    # 47402487
    # .
    # TTGTG
    # T, TTG
    # 1651.19
    # .
    # AC = 1,1;AF = 0.500,0.500;AN = 2;BaseQRankSum = 1.642;ClippingRankSum = 0.235;DP = 272;
    # FS = 0.000;MLEAC = 1, 1;MLEAF = 0.500, 0.500;MQ = 60.20;MQ0 = 0;MQRankSum = -1.718;QD = 6.07;ReadPosRankSum = 1.309
    # GT:AD:DP:GQ:PL
    # 1/2:6,66,39:111:99:2672,1142,1045,1425,0,1678

    # single
    # AC = 1;
    # AF = 0.500;
    # AN = 2;
    # BaseQRankSum = 0.510;
    # ClippingRankSum = 1.677;
    # DP = 215;
    # FS = 2.396;
    # MLEAC = 1;
    # MLEAF = 0.500;
    # MQ = 60.00;
    # MQ0 = 0;
    # MQRankSum = -1.866;
    # QD = 11.55;
    # ReadPosRankSum = -1.493

    # DP=1091;ECNT=5;POP_AF=5.000e-08,5.000e-08,5.000e-08;RPA=10,8,9,11;RU=T;STR;TLOD=3.01,188.84,16.83


    INFO fields:

    :return:
    """

    args = supply_args()
    myvcf = vcfreader.VcfReader(args.input)
    info_number = myvcf.info_number
    samples_number = myvcf.samples_number
    out_vcf_recs = []

    for vrnt in myvcf.myvcf.values():
        if len(vrnt.ALT) > 1:
            split_vrnt = VcfRecDecomp(vrnt, info_number, samples_number).decomp_vrnts
            for svrnt in split_vrnt:
                out_vcf_recs.append(svrnt.print_rec())
        else:
            out_vcf_recs.append(vrnt.print_rec())

    vcfout = vcfwriter.VcfWriter(args.output, out_vcf_recs, myvcf.raw_header).write_me()


if __name__ == "__main__":
    main()
