#!/usr/bin/env python

# DESCRIPTION: VCF files list multiple alternate alleles in the ALT columns,
#  separated by commas. This script will split entries with multiple alleles
# listed, and place each on a separate line.  This will allow us to import
# this data in to annotation software and properly receive prediction
# scores. Post-normalization (vt, bcftools) of this data is highly recommended.
# USAGE: split_mult_alleles.py -h

from copy import deepcopy
import vcfreader
import vcfwriter
import argparse

VERSION = '0.7.2'


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
        self.alts = self._choose_alts(vrnt)
        self.alt_count = len(vrnt.ALT)
        self.gl_map = self._create_gl_map()
        self.info_number = info_number
        self.samples_number = samples_number
        if len(vrnt.rec) >= 9:
            self.vrnt_samples = self._decomp_samples_vrnt(vrnt)
        else:
            self.vrnt_samples = None
        self.vrnt_info = self._decomp_info_vrnt(vrnt)
        self.decomp_vrnts = self._vrnt_update(vrnt)

    def _choose_alts(self, vrnt):
        """
        We don't want to do anything with the asterisks, other than ignore them.
        :return:
        """
        excl = ['*']
        alts = []
        for alt in vrnt.ALT:
            if alt not in excl:
                alts.append(alt)
        return alts

    def _assess_depths(self, samps, alle):
        """
        After a split, a variant may end up with zero depth.  For a single sample VCF, we don't want to
        write this entry if it is zero depth.  For a multi-sample VCF, we would not write the record if all samples
        were of depth zero.
        to file.

        [{'A': {'AD': '0,44', 'PL': '941,567,3942', 'GT': '0/1'}, 'ATT': {'AD': '0,66', 'PL': '941,3942,2650', 'GT': './.'}, 'ATTT': {'AD': '0,0', 'PL': '941,0,4630', 'GT': '0/1'}}, ...]

        :return:
        """
        for samp in samps:
            ad = samp[alle]['AD'].split(',')
            depth = int(ad[0]) + int(ad[1])
            if depth != 0:
                return False
        return True

    def _vrnt_update(self, vrnt, filt_phr='MAsite'):
        """
        Update VcfRecBase obj with split metrics.
        :return:
        """
        vrnt_updt = []
        has_zero_depth = False
        for allele, vals in self.vrnt_info.items():
            new_vrnt = deepcopy(vrnt)
            new_vrnt.ALT = [allele]
            new_vrnt.ID = None
            if vrnt.FILTER and vrnt.FILTER != ['.']:
                new_vrnt.FILTER.append(filt_phr)
            else:
                new_vrnt.FILTER = [filt_phr]
            for val in vals:
                new_vrnt.INFO[val] = vals[val]

            if self.vrnt_samples:
                ind = 0
                has_zero_depth = self._assess_depths(self.vrnt_samples, allele)
                for samp in self.vrnt_samples:
                    for frmt in samp[allele]:
                        new_vrnt.samples[ind][frmt] = samp[allele][frmt]
                    ind += 1

            if not has_zero_depth:
                vrnt_updt.append(new_vrnt)
            else:
                has_zero_depth = False

        return vrnt_updt

    def _decomp_samples_vrnt(self, vrnt):
        """
        Decompose variant in to pieces, to be split.
        split dict:
        {'C': {'GT': './.', 'AD': '0,2', 'PL': '437,258,504'},
        'CA': {'GT': './.', 'AD': '0,1', 'PL': '437,504,417'},
        'CAA': {'GT': './.', 'AD': '0,5', 'PL': '437,261,199'},
        'CAAA': {'GT': './.', 'AD': '0,8', 'PL': '437,378,109'},
        'CAAAA': {'GT': './.', 'AD': '0,0', 'PL': '437,417,569'},
        'CAAAAAA': {'GT': './.', 'AD': '0,0', 'PL': '437,126,569'}}
        :return:
        """
        samp_changes = []
        for samp in vrnt.samples:
            split = {}
            for entry in samp:
                for alt in self.alts:
                    idx = vrnt.ALT.index(alt)
                    if alt not in split:
                        split[alt] = {}
                    if self.samples_number[entry] == 'A':
                        split[alt][entry] = self._decomp_alt(samp[entry], idx)
                    elif self.samples_number[entry] == 'R':
                        split[alt][entry] = self._decomp_ref(samp[entry], idx + 1)
                    elif self.samples_number[entry] == 'G':
                        split[alt][entry] = self._decomp_geno(samp[entry], idx + 1)

            for alta, splits in split.items():
                try:
                    split[alta]['GT'] = self._decomp_gt(samp['GT'], splits['PL'])
                except KeyError:
                    split[alta]['GT'] = "0/1"
            samp_changes.append(split)

        return samp_changes

    def _decomp_info_vrnt(self, vrnt):
        """
        Decompose variant in to pieces, to be split.
        :return:
        """
        split = {}
        for entry in vrnt.INFO:
            for alt in self.alts:
                idx = vrnt.ALT.index(alt)
                if alt not in split:
                    split[alt] = {}

                if self.info_number[entry] == 'A':
                    split[alt][entry] = self._decomp_alt(vrnt.INFO[entry], idx)
                elif self.info_number[entry] == 'R':
                    split[alt][entry] = self._decomp_ref(vrnt.INFO[entry], idx + 1)
                elif self.info_number[entry] == 'G':
                    split[alt][entry] = self._decomp_geno(vrnt.INFO[entry], idx + 1)

        return split

    @staticmethod
    def _decomp_gt(val, pls):
        """
        Split the GT field apart.
        Take PLs:
        437,258,504
        :param val:
        :return:
        """
        if len(val.split('/')) > 2:
            return "0/1"
        else:
            pl = [int(x) for x in pls.split(',')]
            min_pl = min(pl)
            best_gt = pl.index(min_pl)
            if best_gt == 0:
                return "./."
            elif best_gt == 1:
                return "0/1"
            elif best_gt == 2:
                return "1/1"
            else:
                raise Exception("PL field should not contain more than three values.")

    @staticmethod
    def _decomp_alt(val, idx):
        """
        Find fields in INFO that need to be split.
        Input is comma-deliminted value set, can be treated as a string.
        :return:
        """
        return val.split(',')[idx]

    @staticmethod
    def _decomp_ref(val, idx):
        """
        Find fields in SAMPLE columns that need to be split.
        :return:
        """
        ref = val.split(',')[0]
        alt = val.split(',')[idx]
        new_val = [ref, alt]
        return ','.join(new_val)

    def _decomp_geno(self, val, idx):
        """
        Get the indexes that we will pull from the GL values, for each allele.
        :param val:
        :return:
        """
        geno_idxs = [self.gl_map.index((0, idx)),
                     self.gl_map.index((idx, idx))]
        sval = val.split(',')
        new_val = [sval[0], sval[geno_idxs[0]], sval[geno_idxs[1]]]
        return ','.join(new_val)

    def _create_gl_map(self):
        """
        Using the _calc_gl_pos equation, we need to get values based on all possible genotypes.
        :return:
        """
        gl_map = []
        top = 2
        while top <= (self.alt_count + 1):
            for j in range(top):
                for k in range(top):
                    if k >= j:
                        if (j, k) not in gl_map:
                            gl_map.append((j, k))
            top += 1
        return gl_map

    @staticmethod
    def _calc_gl_pos(ref, alt):
        """
        GL : genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible
        genotypes given the set of alleles defined in the REF and ALT fields. In presence of the GT field the same
        ploidy is expected and the canonical order is used; without GT field, diploidy is assumed. If A is the allele in
        REF and B,C,... are the alleles as ordered in ALT, the ordering of genotypes for the likelihoods is given by:
        F(j/k) = (k*(k+1)/2)+j. In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the
        ordering is: AA,AB,BB,AC,BC,CC, etc. For example: GT:GL 0/1:-323.03,-99.29,-802.53 (Floats)
        F(j/k) = (k*(k+1)/2)+j
        j = allele 1
        k = allele 2
        :return:
        """
        return int(((ref * (ref + 1) / 2) + alt))


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
    # FS = 0.000;MLEAC = 1, 1;MLEAF = 0.500, 0.500;MQ = 60.20;MQ0 = 0;MQRankSum = -1.718;QD = 6.07
        ;ReadPosRankSum = 1.309
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
    header_dict = {'id': 'MAsite',
                   'desc': 'Site was split from a multiallelic site.'}
    out_vcf_recs = []

    for vrnt in myvcf.myvcf.values():
        if len(vrnt.ALT) > 1:
            split_vrnt = VcfRecDecomp(vrnt, info_number, samples_number).decomp_vrnts
            for svrnt in split_vrnt:
                if svrnt.samples:
                    uniq_gts = set([x['GT'] for x in svrnt.samples])
                    if uniq_gts != {'./.'}:
                        out_vcf_recs.append(svrnt.print_rec())
                else:
                    out_vcf_recs.append(svrnt.print_rec())
        else:
            out_vcf_recs.append(vrnt.print_rec())

    new_header = vcfwriter.VcfHeader(myvcf.raw_header)
    new_header.add_header_line('FILTER', header_dict)
    vcfwriter.VcfWriter(args.output, out_vcf_recs, new_header.raw_header).write_me()


if __name__ == "__main__":
    main()
