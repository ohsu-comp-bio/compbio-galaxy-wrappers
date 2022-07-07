#!/usr/bin/env python3

# DESCRIPTION: Given a panel of normals VCF, as created by
# panel_of_normals.py, annotate or remove entries in a sample VCF.

from __future__ import print_function
from collections import OrderedDict
from string import Template
import argparse

VERSION = '0.6.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('pon', help='Panel of normals VCF.')
    parser.add_argument('infile', help='Input VCF.')
    parser.add_argument('outfile', help='Output VCF.')
    parser.add_argument('outfile_bad', help='Output VCF with removed variants.')

    # Hard cutoff section
    parser.add_argument('--min_cnt', type=int, help='Variants seen at least this many times in the PON will be '
                                                    'removed, unless they are subject to background assessment.')
    parser.add_argument('--pon_flag_above', help='Flag to apply if variant is seen in PON, does not get removed '
                                                 'for other reasons, and is above min_cnt.')
    parser.add_argument('--pon_flag_below', help='Flag to apply if variant is seen in PON, does not get removed '
                                                 'for other reasons, and is below min_cnt.')
    parser.add_argument('--no_clinvar_rm', action='store_true', help='If the call is in ClinVar, and is not benign, '
                                                                     'do not remove the entry based on '
                                                                     'min_cnt thresholds.')

    # Background filtering option section
    parser.add_argument('--bkgd_avg', help='Average VAF for which we will assess whether variant rises above '
                                           'background. Site will be assessed if value falls below this. If this '
                                           'is specified, bkgd_std must also be specified.')
    parser.add_argument('--bkgd_std', help='VAF stdev for which we will assess whether variant rises above background. '
                                           'Site will be assessed if value falls below this.  '
                                           'If this is specified, bkgd_avg must also be specified.')
    parser.add_argument('--bkgd_pass_flag', help='Flag to apply if background assessment is passed.  No value '
                                                 'here means no flag will be applied.')
    parser.add_argument('--bkgd_min_cnt', default=0, type=int, help='PON variant must be seen at least this number '
                                                                    'of times to be assessed.')

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    if args.bkgd_avg and not args.bkgd_std:
        raise Exception('bkgd_avg and bkgd_std must be specified together.')
    elif not args.bkgd_avg and args.bkgd_std:
        raise Exception('bkgd_avg and bkgd_std must be specified together.')

    return args


class VcfWriter(object):
    """
    Write a VCF.  A header needs to be supplied, which is simply a list of all header lines.
    Also, the body needs to be supplied, which is additionally just a list of ordered lines.
    Creating these formatted lines is done by VcfRecBase.
    """
    def __init__(self, filename, vcf, header):
        self.vcf = vcf
        self.filename = filename
        self.header = header

    def write_me(self):
        """
        Perform writing of MyVcf structure.
        :return:
        """
        with open(self.filename, 'w') as outvcf:
            for entry in self.header:
                outvcf.write(entry)
                outvcf.write('\n')
            for entry in self.vcf:
                self._entry_len_check(entry)
                to_write = '\t'.join(entry)
                outvcf.write(to_write)
                outvcf.write('\n')

    @staticmethod
    def _entry_len_check(entry, length=10):
        """
        If the length of the list is smaller than what we expect, kill the program.
        :return:
        """
        if len(entry) < length:
            raise Exception("VCF record length does not have at least " + str(length) + " columns.")


class VcfRecBase(object):
    """
    Basic parsing of VCF records.
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	RDM-0C19G1-1
    1	16200729	.	A	AT	.	m2	DP=1180;ECNT=1;POP_AF=5e-08;RPA=9;RU=T;STR;TLOD=25.6;OLD_VARIANT
    =1:16200729:AT/ATT	GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:ORIGINAL_CONTIG_MISMATCH:SA_MAP_AF:
    SA_POST_PROB	0/1:906,38:0.031
    :1078:906,38:0,0:36,32:202,190:60,60:38,35:0:0.121,0.121,0.129:0.011,0.002304,0.986
    """
    def __init__(self, rec, args):

        self.rec = rec.rstrip('\n').split('\t')
        self.chrom = str(self.rec[0])
        self.args = args
        # 1-based
        self.coord = str(self.rec[1])
        self.ident = self.rec[2]
        self.ref = self.rec[3]
        self.alt = self.rec[4]
        self.qual = self.rec[5]
        self.filt = self.rec[6].split(';')
        self.info = self._create_info(self.rec[7])

        try:
            self.frmt = self.rec[8].split(':')
            self.samps = self._create_samps(self.rec[9:])
            self.vafs = self._get_samp_vafs()[0]
        except IndexError:
            self.frmt = None
            self.samps = None
            self.vafs = None

        self.uniq_key = (self.chrom, self.coord, self.ref, self.alt)

    def _get_samp_vafs(self):
        """
        Return a list of VAFs per sample in samps.
        :return:
        """
        vafs = []
        for samp, val in self.samps.items():
            try:
                vafs.append(float(val['AF']))
            except KeyError:
                vafs.append(self._get_from_ad(val['AD']))
        return vafs

    @staticmethod
    def _get_from_ad(depths):
        """
        Retrieve VAFs from the AD values, instead of and actual AF sample field.
        Depths are in [ref, alt] format.
        If there are additional alternate alleles, you would have additional entries in the depth list.
        :return:
        """
        ref_cnt = float(depths.split(',')[0])
        alt_cnt = float(depths.split(',')[1])
        vaf = alt_cnt / (alt_cnt + ref_cnt)
        return vaf

    @staticmethod
    def _retr_write(val, delim):
        """
        Return writable version of basic delimited strings.
        :return:
        """
        return delim.join(val)

    def _retr_info_write(self):
        """
        Return writable version of INFO string.
        :return:
        """
        info = []
        for key, val in self.info.items():
            if val:
                info.append('='.join([key, val]))
            else:
                info.append(key)
        return ';'.join(info)

    def _retr_samps_write(self):
        """
        Return SAMPLE columns in writable format.
        Looks like {<SAMPLE>: {<FRMT>: <SAMP_VAL>} }
        :return:
        """
        samps = []
        for data in self.samps.values():
            vals = []
            for val in data.values():
                vals.append(val)
            samps.append(':'.join(vals))
        return samps

    def retr_curr_values(self):
        """
        Return current state of vcfkdl rec.
        :return:
        """
        to_write = [self.chrom, self.coord, self.ident, self.ref, self.alt, self.qual,
                    self._retr_write(self.filt, delim=';'),
                    self._retr_info_write()]
        if self.frmt:
            to_write.append(self._retr_write(self.frmt, delim=':'))
            to_write.extend(self._retr_samps_write())
        return to_write

    def _create_samps(self, samps):
        """

        :return:
        """
        new_samps = OrderedDict()
        for samp in samps:
            if samp:
                new_samps[samp] = OrderedDict()
                samp_data = samp.split(':')
                for val in samp_data:
                    idx = samp_data.index(val)
                    new_samps[samp][self.frmt[idx]] = val
        return new_samps

    @staticmethod
    def _create_info(info):
        """
        Place INFO data in structure.
        :return:
        """
        new_info = OrderedDict()
        for entry in info.split(';'):
            key = entry.split('=')[0]
            try:
                val = entry.split('=')[1]
            except:
                val = None
            new_info[key] = val
        return new_info

    def _value_set(self, value):
        """
        Set the variable to something in the INFO field.
        :return:
        """
        try:
            return self.info[value]
        except:
            return None

    def var_req(self):
        """
        Print out the current value of each variable.
        :param self:
        :return:
        """
        print("PYVCF REC: {0}".format(self.rec))
        print("PYVCF INFO: {0}".format(self.info))


class MyVcf(object):
    """
    Hold VCF records, manage access to them based in chrom, pos, ref, alt.
    """
    def __init__(self, filename, vcfrec, args):
        self.filename = filename
        self.args = args
        self.vcfrec = vcfrec
        self.header = []
        self.myvcf = self._create_vcf()

    def _create_vcf(self):
        """
        VCF INFO field
        AVG_AF=0.585;STDEV_AF=0.357;MAX_AF=0.985;MIN_AF=0.024;COSMIC=F;SEEN=18;TLOD=1774.033;CLNSIG=F
        :return:
        """
        myvcf = OrderedDict()
        with open(self.filename, 'r') as vcf:
            for line in vcf:
                if not line.startswith('#'):
                    vrnt = self._sel_vcfrec(self.vcfrec, line, self.args)
                    if vrnt.uniq_key not in myvcf:
                        myvcf[vrnt.uniq_key] = vrnt
                    else:
                        raise Exception("No duplicates allowed.  Make sure your VCF doesn't contain them.")
                else:
                    self.header.append(line.rstrip('\n'))
        return myvcf

    @staticmethod
    def _sel_vcfrec(vcfrec, val, args):
        """
        Supply correct class call.
        :return:
        """
        return vcfrec(val, args)


class VcfRecPON(VcfRecBase):
    """
    Specific to PON INFO entries.
    For PONs, INFO fields currently look like this:
    AVG_AF=0.021;STDEV_AF=0.003;MAX_AF=0.027;MIN_AF=0.016;COSMIC=F;SEEN=35;TLOD=10.257;CLNSIG=F
    """
    def __init__(self, rec, args):
        super(VcfRecPON, self).__init__(rec, args)
        self.min_cnt = args.min_cnt
        self.bkgd_avg = args.bkgd_avg
        self.bkgd_std = args.bkgd_std
        self.bkgd_min_cnt = args.bkgd_min_cnt
        # Set variables that are included in INFO line.
        self.avg_ab = self._value_set('AVG_AF')
        self.stdev_ab = self._value_set('STDEV_AF')
        self.max_ab = self._value_set('MAX_AF')
        self.min_ab = self._value_set('MIN_AF')
        self.cosmic = self._value_set('COSMIC')
        self.seen = self._value_set('SEEN')
        self.tlod_avg = self._value_set('TLOD')
        self.clnsig = self._value_set('CLNSIG')

        # Assess variables based on additional criteria, if they have been passed.
        self.assess_bkgd = self._bkgd_assess()
        self.remove_me = self._remove_assess()
        self.pon_label_me = self._label_assess()
        self.is_cosmic = self._assess_bool(self.cosmic)
        self.is_clinvar = self._assess_clinvar(self.clnsig)

    @staticmethod
    def _assess_clinvar(val):
        """
        For ClinVar, the entries that are not benign are interesting.
        :return:
        """
        # If ClinVar labels are in this list, they can be treated just like CLNSIG=F PON variants.
        bad_clinvar = ['Likely_benign', 'Benign', 'Benign/Likely_benign']
        if val == 'F':
            return False
        elif val in bad_clinvar:
            return False
        else:
            return True

    @staticmethod
    def _assess_bool(val):
        """
        If value is set as T, return True.  Otherwise, return False.
        :return:
        """
        if val == 'T':
            return True
        else:
            return False

    def _label_assess(self):
        """
        Determine whether seeing this PON entry will cause it to be labeled with something in the FILTER column.
        This will happen if remove_me is False, and if assess_bkgd is False.  If assess_bkgd is True, there should only
        be two options, either removal, or inclusion with flag.
        :return:
        """
        if not self.assess_bkgd and not self.remove_me:
            return True
        return False

    def _remove_assess(self):
        """
        Determine whether the variant will be an absolute removal.
        :return:
        """
        if int(self.seen) >= int(self.min_cnt):
            return True
        else:
            return False

    def _bkgd_assess(self):
        """
        when avg_af < 0.2 and stdev_af < 0.06
          check whether vaf exceeds (3 * stdev_af + avg_af)
        if yes:
          add flag (PON_OV, user defined)
        if no:
          don't pass through
        Propose flags called (--bkgd_avg, --bkgd_std, bkgd_fail_flag, bkgd_pass_flag, bkgd_min_cnt)
        Where fail implies the call's VAF does not exceed threshold.
        If these flag args are not passed, don't add a flag
        Background min count will govern when this process can actually be applied, based on COUNT of variant in PON.
        Above this value will mean the call can be assessed.  Below, and it will not.
        :return:
        """
        if float(self.avg_ab) < float(self.bkgd_avg):
            if float(self.stdev_ab) < float(self.bkgd_std):
                if int(self.seen) >= int(self.bkgd_min_cnt):
                    return True
        return False

    def var_req(self):
        """
        Print out the current value of each variable.
        :param self:
        :return:
        """
        super(VcfRecPON, self).var_req()
        print("Assess Background?: {0}".format(self.assess_bkgd))
        print("Hard Remove?: {0}".format(self.remove_me))


class MyVcfEdited(MyVcf):
    """

    """
    def __init__(self, filename, vcfrec, args, low_seen):
        super(MyVcfEdited, self).__init__(filename, vcfrec, args)
        self.param_dict = {'bkgd_min_cnt': args.bkgd_min_cnt,
                           'bkgd_pass_flag': args.bkgd_pass_flag,
                           'pon_flag_below': args.pon_flag_below,
                           'pon_flag_above': args.pon_flag_above,
                           'min_cnt': args.min_cnt,
                           'low_seen': low_seen}
        self.new_headers = self._create_new_headers()
        self.new_header = self._insert_new_headers()

    def _insert_new_headers(self):
        """
        Locate the place these new header lines belong.  Currently, just place them before the first
        column heading line.
        :return:
        """
        return self.header[:-1] + self.new_headers + self.header[-1:]

    def _create_new_headers(self):
        """

        :return:
        """
        new_headers = []
        for param, val in self.param_dict.items():
            if val:
                try:
                    new_headers.append(self._filter_create(param))
                except:
                    pass
        return new_headers

    def _filter_create(self, tmpl):
        """
        Create FILTER lines to be added for each new label.
        :return:
        """
        filt_tmpls = {'bkgd_pass_flag': Template('##FILTER=<ID=${bkgd_pass_flag},Description=\"Variant found in '
                                                 'panel of normals in ${bkgd_min_cnt} or greater samples, '
                                                 'and exceeds background at this site.\">'),
                      'pon_flag_below': Template('##FILTER=<ID=${pon_flag_below},Description=\"Variant found in '
                                                 'panel of normals in ${low_seen} or more and fewer '
                                                 'than ${min_cnt} samples.\">'),
                      'pon_flag_above': Template('##FILTER=<ID=${pon_flag_above},Description=\"Variant found in '
                                                 'panel of normals in ${min_cnt} or more samples.\">')}
        return filt_tmpls[tmpl].substitute(self.param_dict)


class PanelOfNormals(MyVcf):
    """
    COSMIC variant handling
    non-COSMIC variant handling
    low VAF variant handling
    applying FILTER tags
    """
    def __init__(self, filename, vcfrec, args):
        super(PanelOfNormals, self).__init__(filename, vcfrec, args)
        self.min_cnt = args.min_cnt
        self.bkgd_min_cnt = args.bkgd_min_cnt
        self.bkgd_pass_flag = args.bkgd_pass_flag
        self.pon_flag_below = args.pon_flag_below
        self.pon_flag_above = args.pon_flag_above
        self.low_seen = self._find_low_seen()

    def _find_low_seen(self):
        """
        VCF INFO field
        AVG_AF=0.585;STDEV_AF=0.357;MAX_AF=0.985;MIN_AF=0.024;COSMIC=F;SEEN=18;TLOD=1774.033;CLNSIG=F
        Find the lowest SEEN value.
        :return:
        """
        low_seen = 1000000
        for entry in self.myvcf.values():
            if int(entry.info['SEEN']) < low_seen:
                low_seen = int(entry.info['SEEN'])
        return low_seen

    @staticmethod
    def _delim_to_dict(info, delim=';'):
        """
        Create a dict out of a typical VCF INFO field string.
        :return:
        """
        new_dict = {}
        for entry in info.split(delim):
            field = entry.split('=')[0]
            val = entry.split('=')[1]
            new_dict[field] = val
        return new_dict


class VcfRecComp(object):
    """
    Decide on whether the entry should be written.
    The rules:
        if the entry is in the PON, at above min_cnt, we can either remove it, or just label it
            min_cnt implies either removal or label of variants seen above min_cnt times
        if the entry is in the PON, at below min_cnt, it can be labeled, though probably not removed
        (otherwise no need to specify min_cnt)

    """
    def __init__(self, pon_rec, vcf_rec, bkgd_pass_flag, pon_flag_below, no_clinvar_rm):
        self.pon_rec = pon_rec
        self.vcf_rec = vcf_rec
        self.bkgd_pass_flag = bkgd_pass_flag
        self.pon_flag_below = pon_flag_below
        self.no_clinvar_rm = no_clinvar_rm

        self.add_flags = []
        self.write_me = True

        # If the record is set as a remove candidate, based on min_cnt, set write_me to False.
        if self.pon_rec.remove_me:
            # If arg has been sent to not remove ClinVar variants...
            if self.no_clinvar_rm:
                if self.pon_rec.is_clinvar:
                    self.write_me = True
                else:
                    self.write_me = False
            else:
                self.write_me = False
        else:
            self.write_me = True

        # If the entry is classified as needing assessment, do it here.
        if self.pon_rec.assess_bkgd:
            self.over_bkgd = self._assess_bkgd()
            if self.over_bkgd:
                self._add_flag(self.bkgd_pass_flag)
                self.write_me = True
            else:
                self.write_me = False

        # Check if the record is set for labeling.  This will occur if the number of times SEEN is less than
        # min_cnt and bkgd_min_cnt.
        if self.pon_rec.pon_label_me:
            self._add_flag(self.pon_flag_below)

        for entry in self.add_flags:
            self.vcf_rec.filt = self._replace_filter(self.vcf_rec.filt, entry)

    def _add_flag(self, arg):
        """
        Check if we are adding a flag.
        :return:
        """
        if arg:
            self.add_flags.append(arg)

    def _assess_bkgd(self, min_stdev=0.01):
        """
        Perform the background assessment.
        :return:
        """
        thresh = (3.0 * (max(float(self.pon_rec.stdev_ab), min_stdev)) + float(self.pon_rec.avg_ab))
        if float(self.vcf_rec.vafs) >= thresh:
            return True
        else:
            return False

    @staticmethod
    def _replace_filter(filt, anno):
        """
        Remove the dot if it's there, otherwise just append something with a
        semicolon.
        :return:
        """
        pass_flags = ['.', 'PASS']
        for pss in pass_flags:
            if pss in filt:
                return [anno]
        if anno not in filt:
            filt.append(anno)
        return filt


class CompVcf(object):
    """
    Deal with comparison portion.
    """
    def __init__(self, myvcf, pon, args):
        self.myvcf = myvcf
        self.pon = pon
        self.bkgd_pass_flag = args.bkgd_pass_flag
        self.pon_flag_below = args.pon_flag_below
        self.no_clinvar_rm = args.no_clinvar_rm
        self.good_vcf, self.bad_vcf = self._start_comp()

    def _do_comp(self, vcf, pon):
        """
        Do the actual comparisons, according to specified logic.
        :return:
        """
        comp = VcfRecComp(pon, vcf, self.bkgd_pass_flag, self.pon_flag_below, self.no_clinvar_rm)
        if comp.write_me:
            return True
        else:
            return False

    def _start_comp(self):
        """
        Kick off the comparison process.
        :return:
        """
        good_vcf = []
        bad_vcf = []
        for coord in self.myvcf.myvcf:
            if coord in self.pon.myvcf:
                vcfv = self.myvcf.myvcf[coord]
                ponv = self.pon.myvcf[coord]
                if self._do_comp(vcfv, ponv):
                    good_vcf.append(self.myvcf.myvcf[coord].retr_curr_values())
                else:
                    bad_vcf.append(self.myvcf.myvcf[coord].retr_curr_values())
            else:
                good_vcf.append(self.myvcf.myvcf[coord].retr_curr_values())
        return good_vcf, bad_vcf


def main():

    args = supply_args()
    my_pon = PanelOfNormals(args.pon, VcfRecPON, args)
    my_vcf = MyVcfEdited(args.infile, VcfRecBase, args, my_pon.low_seen)
    good_vcf = CompVcf(my_vcf, my_pon, args).good_vcf
    bad_vcf = CompVcf(my_vcf, my_pon, args).bad_vcf
    VcfWriter(args.outfile, good_vcf, my_vcf.new_header).write_me()
    VcfWriter(args.outfile_bad, bad_vcf, my_vcf.header).write_me()


if __name__ == "__main__":
    main()
