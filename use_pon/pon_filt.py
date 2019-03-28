#!/usr/bin/env python

# DESCRIPTION: Given a panel of normals VCF, as created by
# panel_of_normals.py, annotate or remove entries in a sample VCF.
# USAGE: use_pon.py <pon> <infile> <outfile>
# CODED BY: John Letaw

from __future__ import print_function
from string import Template
import argparse
import numpy
import vcf

VERSION = '0.5.0'

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
    parser.add_argument('--min_cnt', type=int, help='Variants seen at least this many times in the PON will be removed, unless they are subject to background assessment.')
    parser.add_argument('--pon_flag', help='Flag to apply if variant is seen in PON, but does not get removed for other reasons.')

    # Background filtering option section
    parser.add_argument('--bkgd_avg', help='Average VAF for which we will assess whether variant rises above background.  '
                                           'Site will be assessed if value falls below this.  '
                                           'If this is specified, bkgd_std must also be specified.')
    parser.add_argument('--bkgd_std', help='VAF stdev for which we will assess whether variant rises above background. '
                                           'Site will be assessed if value falls below this.  '
                                           'If this is specified, bkgd_avg must also be specified.')
    parser.add_argument('--bkgd_pass_flag', help='Flag to apply if background assessment is passed.  No value here means no flag will be applied.')
    parser.add_argument('--bkgd_min_cnt', default=0, type=int, help='PON variant must be seen at least this number of times to be assessed.')

    parser.add_argument('--hotspots', required=False, help='Hotspots VCF.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    if args.bkgd_avg and not args.bkgd_std:
        raise Exception('bkgd_avg and bkgd_std must be specified together.')
    elif not args.bkgd_avg and args.bkgd_std:
        raise Exception('bkgd_avg and bkgd_std must be specified together.')

    return args


class PanelOfNormals(object):
    """
    COSMIC variant handling
    non-COSMIC variant handling
    low VAF variant handling
    applying FILTER tags
    """
    def __init__(self, args):
        self.filename = args.pon
        self.args = args
        self.pon = self._create_pon()

    def _create_pon(self):
        """
        VCF INFO field
        AVG_AF=0.585;STDEV_AF=0.357;MAX_AF=0.985;MIN_AF=0.024;COSMIC=F;SEEN=18;TLOD=1774.033;CLNSIG=F
        :return:
        """
        pon = {}
        for variant in vcf.Reader(open(self.filename, 'r')):
            uniq_key = (str(variant.CHROM), str(variant.POS), str(variant.REF), str(variant.ALT[0]))
            if uniq_key not in pon:
                pon[uniq_key] = VcfRecPON(variant, self.args.min_cnt, self.args.bkgd_min_cnt, self.args.bkgd_avg, self.args.bkgd_std)
            else:
                raise Exception("No duplicates allowed.  Make sure PON doesn't contain them.")
        return pon

    def _delim_to_dict(self, info, delim=';'):
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


class VcfRec(object):
    """
    """
    def __init__(self, rec):

        self.rec = rec
        self.chrom = str(rec.CHROM)
        self.coord = str(rec.POS)

    def _value_set(self, value):
        """
        Set the variable to something in the INFO field.
        :return:
        """
        try:
            return self.rec.INFO[value]
        except:
            return None

    def _calc_ab(self):
        """
        Get the average allele balance for a specific genotype.
        :return:
        """
        all_vafs = []
        for entry in self.rec.samples:
            if entry['GT'] == '0/1':
                ref_cnt = entry['AD'][0]
                alt_cnt = entry['AD'][1]
                vaf = alt_cnt / (alt_cnt + ref_cnt + 0.0)
                all_vafs.append(vaf)

        try:
            return numpy.mean(all_vafs), numpy.std(all_vafs)
        except:
            return 0.0, 0.0

    def var_req(self):
        """
        Print out the current value of each variable.
        :param self:
        :return:
        """
        print("PYVCF REC: {0}".format(self.rec))
        print("PYVCF INFO: {0}".format(self.rec.INFO))
        # print("Chromosome: {0}".format(self.chrom))
        # print("Position: {0}".format(self.coord))

class VcfRecPON(VcfRec):
    """
    Specific to PON INFO entries.
    For PONs, INFO fields currently look like this:
    AVG_AF=0.021;STDEV_AF=0.003;MAX_AF=0.027;MIN_AF=0.016;COSMIC=F;SEEN=35;TLOD=10.257;CLNSIG=F
    """
    def __init__(self, rec, min_cnt, bkgd_min_cnt=0, bkgd_avg=None, bkgd_std=None):
        super(VcfRecPON, self).__init__(rec)
        self.min_cnt = min_cnt
        self.bkgd_avg = bkgd_avg
        self.bkgd_std = bkgd_std
        self.bkgd_min_cnt = bkgd_min_cnt
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
        try:
            if int(self.seen) >= int(self.min_cnt):
                if not self.assess_bkgd:
                    return True
        except:
            return False
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
        try:
            if float(self.avg_ab) < float(self.bkgd_avg):
                if float(self.stdev_ab) < float(self.bkgd_std):
                    if int(self.seen) >= int(self.bkgd_min_cnt):
                        return True
        except:
            return False
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


class PanelOfNormalsWriter():
    """
    Handle writing of new PON.
    """
    def __init__(self, infile, outfile, outfile_bad, pon, bkgd_min_cnt=None, bkgd_pass_flag=None, pon_flag=None, min_cnt=None):
        self.infile = infile
        self.outfile = outfile
        self.outfile_bad = outfile_bad
        self.pon = pon
        self.min_cnt = min_cnt
        self.bkgd_min_cnt = bkgd_min_cnt
        self.bkgd_pass_flag = bkgd_pass_flag
        self.pon_flag = pon_flag
        self.param_dict = {'bkgd_min_cnt': self.bkgd_min_cnt,
                           'bkgd_pass_flag': self.bkgd_pass_flag,
                           'pon_flag': self.pon_flag,
                           'min_cnt': self.min_cnt}


    def _replace_filter(self, filt, anno):
        """
        Remove the dot if it's there, otherwise just append something with a
        semicolon.
        :return:
        """
        if anno:
            if filt == '.' or filt == 'PASS':
                return anno
            else:
                return ';'.join([filt, anno])
        return filt

    def _assess_bkgd(self, rec, vaf):
        """
        Perform the background assessment.
        :return:
        """
        if vaf >= ((3 * rec.stdev_ab) + rec.avg_ab):
            return True
        else:
            return False

    def _filter_create(self, tmpl):
        """
        Create FILTER lines to be added for each new label.
        :return:
        """
        filt_tmpls = {'bkgd_pass_flag': Template('##FILTER=<ID=${bkgd_pass_flag},Description=\"Variant found in panel of normals in ${bkgd_min_cnt} or greater samples, and exceeds background at this site.\">\n'),
                      'pon_flag': Template('##FILTER=<ID=${pon_flag},Description=\"Variant found in panel of normals in less than ${min_cnt} samples.\">\n')}
        return filt_tmpls[tmpl].substitute(self.param_dict)

    def write_new_vcf(self):
        """
        Write a new VCF with the panel_of_normals tag in FILTER section.
        :return:
        """
        outfile = open(self.outfile, 'w')
        outfile_bad = open(self.outfile_bad, 'w')
        check = True
        with open(self.infile, 'rU') as infile:
            for line in infile:
                if not line.startswith('#'):
                    sline = line.rstrip('\n').split('\t')
                    labels = sline[8].split(':')
                    info = sline[9].split(':')
                    vaf = float(info[labels.index('AF')])
                    uniq_key = (sline[0], sline[1], sline[3], sline[4])
                    write_me = True
                    if uniq_key in self.pon:
                        this_rec = self.pon[uniq_key]
                        # Check if PON variant is labeled as needing assessment.
                        if this_rec.assess_bkgd:
                            # Check is VAF falls outside of PON background range.
                            if self._assess_bkgd(this_rec, vaf):
                                sline[6] = self._replace_filter(sline[6], self.bkgd_pass_flag)
                                write_me = True
                            else:
                                write_me = False
                        # Check to see if the record should be removed, based on min_cnt.
                        if this_rec.remove_me:
                            write_me = False
                        # Check to see if we should label the record.
                        if this_rec.pon_label_me:
                            sline[6] = self._replace_filter(sline[6], self.pon_flag)
                            write_me = True
                    if write_me:
                        outfile.write('\t'.join(sline))
                        outfile.write('\n')
                    else:
                        outfile_bad.write('\t'.join(sline))
                        outfile_bad.write('\n')

                else:
                    outfile.write(line)
                    outfile_bad.write(line)
                    # Move to own function
                    if line.startswith("##FILTER") and check:
                        if self.bkgd_pass_flag:
                            outfile.write(self._filter_create('bkgd_pass_flag'))
                        if self.pon_flag:
                            outfile.write(self._filter_create('pon_flag'))
                        check = False

        outfile.close()
        outfile_bad.close()


class Hotspots():
    """
    Import HOTSPOTS file, if given, to idential structure as for PON.
    """
    def __init__(self, filename):
        self.filename = filename
        self.hotspots = self._hotspots_grab()

    def _hotspots_grab(self):
        """

        :param filename:
        :return:
        """
        hotspots = []
        with open(self.filename, 'rU') as hots:
            for line in hots:
                if not line.startswith('#'):
                    line = line.rstrip('\n').split('\t')
                    chrom = line[0]
                    pos = line[1]
                    ref = line[3]
                    alt = line[4]
                    uniq_key = (chrom, pos, ref, alt)
                    if uniq_key not in hotspots:
                        hotspots.append(uniq_key)
        return hotspots


def main():

    # Possible filters:

    # Calls not assessed as above will need to go through normal filtering process. Parameters should be set for
    # number of times seen (SEEN), and flag names should be specified.  When looking at number of times seen, user can
    # choose to include or exclude variants in COSMIC.  Similar criteria can be set (eg SEEN) for these COSMIC or ClinVar
    # instances.


    args = supply_args()
    # if args.hotspots:
    #     hotspots = Hotspots(args.hotspots).hotspots
    # else:
    #     hotspots = None
    my_pon = PanelOfNormals(args)

    # test_var = ('1', '115251259', 'A', 'G')
    # my_var = my_pon.pon[test_var]
    # my_var.var_req()
    # print(my_var.seen)
    # print(my_var.min_cnt)
    # print(my_var.assess_bkgd)

    PanelOfNormalsWriter(args.infile, args.outfile, args.outfile_bad, my_pon.pon,
                         args.bkgd_min_cnt, args.bkgd_pass_flag, args.pon_flag,
                         args.min_cnt).write_new_vcf()


if __name__ == "__main__":
    main()
