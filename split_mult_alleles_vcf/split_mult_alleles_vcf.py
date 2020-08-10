#!/usr/bin/env python

# DESCRIPTION: VCF files list multiple alternate alleles in the ALT columns,
#  separated by commas. This script will split entries with multiple alleles
# listed, and place each on a separate line.  This will allow us to import
# this data in to annotation software and properly receive prediction
# scores. Post-normalization (vt, bcftools) of this data is highly recommended.
# USAGE: split_mult_alleles.py -h

from copy import deepcopy
import argparse
import vcfpy

VERSION = '0.6.0'


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

class VcfRec:
    """
    {'AC': [48, 25],
    'AF': [0.276, 0.144],
    'AN': 174,
    'BaseQRankSum': -0.74,
    'ClippingRankSum': 0.0,
    'DP': 25174,
    'ExcessHet': 160.0,
    'FS': 5.155,
    'InbreedingCoeff': -0.7756,
    'MLEAC': [50, 26],
    'MLEAF': [0.287, 0.149],
    'MQ': 58.95, 'MQRankSum': 0.0,
    'QD': 5.34,
    'ReadPosRankSum': 2.34,
    'SOR': 1.091}
    """
    def __init__(self, vcfpy_rec, idx, info_types, frmt_types):
        self.vcfpy_rec = vcfpy_rec
        self.idx = idx
        self.info_types = info_types
        self.frmt_types = frmt_types
        self.genos = self._genos_metrics()

    def create_new_rec(self):
        new_rec = deepcopy(self.vcfpy_rec)
        new_rec.ALT = self._alt()
        new_rec.INFO = self._info(self.vcfpy_rec.INFO)
        new_rec.calls = self._sample_upd(self.vcfpy_rec.calls)
        return new_rec

    def _sample_upd(self, calls):
        """

        :param calls:
        :return:
        """
        new_calls = []
        for call in deepcopy(calls):
            new_call = deepcopy(call)
            for field, metric in new_call.data.items():
                ntype = self.frmt_types[field]
                if field != 'GT':
                    new_call.data[field] = self._info_upd(ntype, metric)
            new_call.data['GT'] = self._assign_gt(new_call.data['PL'])
            new_calls.append(new_call)
        return deepcopy(new_calls)

    def _assign_gt(self, pl):
        """
        Reassign genotypes, where appropriate.  If there are too many
        :param pl:
        :return:
        """
        if pl.count(min(pl)) == 1:
            if pl.index(min(pl)) == 1:
                return '0/1'
            elif pl.index(min(pl)) == 2:
                return '1/1'
            else:
                return '0/0'
        return '0/0'

    def _info(self, info):
        """
        Copy and create a new INFO OrderedDict.
        :param info:
        :return:
        """
        new_info = deepcopy(info)
        for field in new_info:
            ntype = self.info_types[field]
            new_info[field] = self._info_upd(ntype, new_info[field])
        return new_info

    def _info_upd(self, ntype, metrics):
        """
        Based on the Number type defined in the header, return metrics for single variant.
        :param ntype:
        :param metrics:
        :return:
        """
        if ntype == 'A':
            return [metrics[self.idx]]
        elif ntype == 'R':
            return [metrics[0], metrics[self.idx+1]]
        elif ntype == 'G':
            return [metrics[x] for x in self.genos]
        else:
            return metrics

    def _genos_metrics(self):
        """
        Collect the actual indices for number=G metrics.
        :return:
        """
        geno_mets = []
        for entry in self._genos_allowed():
            geno_mets.append(self._calc_gl_pos(entry[0], entry[1]))
        return geno_mets

    def _genos_allowed(self):
        """
        Find tuples for all possible genotypes connected to this particular allele.
        :return:
        """
        geno_idx = self.idx+1
        return [(0, 0), (0, geno_idx), (geno_idx, geno_idx)]

    @staticmethod
    def _calc_gl_pos(ref, alt):
        """
        F(j/k) = (k*(k+1)/2)+j
        j = allele 1
        k = allele 2
        :return:
        """
        return int(((ref*(ref+1)/2)+alt))

    def _alt(self):
        """
        Replace the field in ALT with one of the alleles from a multiallelic entry.
        :return:
        """
        return [self.vcfpy_rec.ALT[self.idx]]

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
    vcfreader = vcfpy.Reader.from_path(args.input)

    info_types = {}
    for info in vcfreader.header.info_ids():
        info_types[info] = vcfreader.header.get_info_field_info(info).number

    frmt_types = {}
    for frmt in vcfreader.header.format_ids():
        frmt_types[frmt] = vcfreader.header.get_format_field_info(frmt).number

    vcfwriter = vcfpy.Writer.from_path(args.output, vcfreader.header)
    for entry in vcfreader:
        # if entry.POS == 48032740:
        # print(entry.ALT)
        # print(len(entry.ALT))
        if len(entry.ALT) > 1:
            # print(entry)
            # GL actually '.'
            # CA,GWAS_PUBMED,EXOME_CHIP,EA_AGE,AA_AGE,GRCh38_POSITION is unknown, check
            # per_all = ['EA_AC', 'AA_AC', 'TAC'] #looks good
            # # MAF is always 3
            # three = ['MAF']
            # dots = ['GL', 'FG', 'HGVS_CDNA_VAR', 'HGVS_PROTEIN_VAR', 'CDS_SIZES', 'GS', 'PH', 'GTS', 'EA_GTC', 'AA_GTC', 'GTC']
            # per_alt = []
            # per_geno = []
            # all = len(entry.ALT) + 1
            # a = len(entry.ALT)
            # c = ((a*(a+1)/2)+a)+1
            # for line in per_geno:
            #     b = len(entry.INFO[line])
            #     if int(c) != b:
            #         print(entry)
            #         print(line)
            #         print(entry.ALT)
            #         print(entry.INFO[line])
            #         print(c)
            #         print(b)
            #         exit(1)
            for alle in range(len(entry.ALT)):
                alt = entry.ALT[alle].value
                # We can ignore the asterisks, since they refer to upstream events, which we will be dealing with.
                if alt != '*':
                    new_entry = VcfRec(entry, alle, info_types, frmt_types).create_new_rec()
                    # print(new_entry)
                    vcfwriter.write_record(new_entry)
                    # print(new_entry)
                # new_entry = deepcopy(entry)
                # new_entry.ALT = [alle]
                # print(new_entry)
        # else:
        #     vcfwriter.write_record(entry)

    vcfwriter.close()

    # info_number = myvcf.info_number
    # samples_number = myvcf.samples_number
    # out_vcf_recs = []
    #
    # for vrnt in myvcf.myvcf.values():
    #     if len(vrnt.ALT) > 1:
    #         split_vrnt = VcfRecDecomp(vrnt, info_number, samples_number).decomp_vrnts
    #         for svrnt in split_vrnt:
    #             out_vcf_recs.append(svrnt.print_rec())
    #     else:
    #         out_vcf_recs.append(vrnt.print_rec())
    #
    # vcfout = vcfwriter.VcfWriter(args.output, out_vcf_recs, myvcf.raw_header).write_me()


if __name__ == "__main__":
    main()
