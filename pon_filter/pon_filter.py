#!/usr/bin/env python3

"""
Given a panel of normals VCF, as created by panel_of_normals.py, annotate entries in a sample VCF.

Usage:
python pon_filter.py tests/test_data/input.vcf tests/test_data/pon.vcf tests/test_data/output.vcf \
    --min_cnt "40" \
    --pon_flag "PON" \
    --pon_flag_filtered "PON_FILT" \
    --bkgd_avg "0.2" \
    --bkgd_std "0.06" \
    --bkgd_min_cnt "4" \
    --pon_flag_above_bkgd "PON_OV"

Details:
Adapted from pon_filter/pon_filt.py to use vcfpy library and simplify

"""


import argparse
import vcfpy
import sys

VERSION = '1.0.0'


class VcfRecord:
    def __init__(self, record):
        self.record = record
        self.id = '{}:{}{}>{}'.format(self.record.CHROM, self.record.POS, self.record.REF, self.record.ALT)

    def merge(self, new_record):
        for filt in new_record.FILTER:
            if filt not in self.record.FILTER:
                self.record.FILTER.append(filt)
        for info in new_record.INFO:
            self.record.INFO[info] = new_record.INFO[info]
        return self.record

    def merge2(self, new_record):
        field_line_types = ['FILTER', 'INFO']
        for ftype in field_line_types:
            for fcol in new_record[ftype]:
                if fcol not in self.record[ftype]:
                    self.record[ftype][fcol] = new_record[ftype][fcol]
        return self.record


class PonRecord(VcfRecord):
    def __init__(self, record, min_cnt, bkgd_avg, bkgd_std, bkgd_min_cnt):
        super().__init__(record)
        self.avg_af = record.INFO['AVG_AF']
        self.stdev_af = record.INFO['STDEV_AF']
        self.max_af = record.INFO['MAX_AF']
        self.min_af = record.INFO['MIN_AF']
        self.cosmic = record.INFO['COSMIC']
        self.seen = record.INFO['SEEN']
        self.tlod = record.INFO['AVG_TLOD']
        self.clinsig = record.INFO['CLNSIG']

        self.min_cnt = min_cnt
        self.bkgd_avg = bkgd_avg
        self.bkgd_std = bkgd_std
        self.bkgd_min_cnt = bkgd_min_cnt

        self.above_min_cnt = self._above_min_cnt()
        self.below_bkgd_thresholds = self._below_bkgd_thresholds()
        self.get_bkgd_threshold = self._get_bkgd_threshold()

    def _above_min_cnt(self):
        return self.seen >= self.min_cnt

    def _below_bkgd_thresholds(self):
        if float(self.avg_af) < float(self.bkgd_avg):
            if float(self.stdev_af) < float(self.bkgd_std):
                if int(self.seen) >= int(self.bkgd_min_cnt):
                    return True
        return False

    def _get_bkgd_threshold(self, min_stdev=0.01):
        return 3.0 * (max(float(self.stdev_af), min_stdev)) + float(self.avg_af)

    def above_bkgd(self, sample):
        return float(self.record.call_for_sample[sample].data['AF']) >= self.get_bkgd_threshold


class VcfRecords:
    def __init__(self, input_vcf):
        self.reader = vcfpy.Reader.from_path(input_vcf)
        self.header = self.reader.header
        self.sample = self._get_sample()
        self.records = self._as_dict()

    def _as_dict(self):
        records = {}
        for record in self.reader:
            record = VcfRecord(record)
            if record.id not in records:
                records[record.id] = record.record
        return records

    def _get_sample(self):
        try:
            return self.reader.header.samples.names[0]
        except IndexError:
            # e.g No sample VCF. Resource VCF?
            return None

    def merge_header(self, new_header):
        header_line_types = ['contig', 'FILTER', 'INFO']
        for header_line_type in header_line_types:
            for hline in new_header.get_lines(header_line_type):
                if not self.header.has_header_line(header_line_type, hline.id):
                    self.header.add_line(hline)
        return self.header

    @staticmethod
    def update_header(header, hid, hdesc):
        header.add_filter_line(vcfpy.OrderedDict([('ID', hid),
                                                  ('Description', hdesc)]))
        return header


def get_args(args_input):
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', help='Input VCF.')
    parser.add_argument('pon_vcf', help='Panel of normals VCF.')
    parser.add_argument('output_vcf', help='Output VCF.')
    parser.add_argument('--min_cnt', type=int, help='A variant seen at least this many times in the PON will be removed, '
                                                    'unless it is bkgd assessed and found to be above background.')
    parser.add_argument('--pon_flag', type=str, help='A variant will have this flag if seen in the PON'
                                                     'but is below min_cnt and is not bkgd accessed.')
    parser.add_argument('--pon_flag_filtered', type=str, help='A variant will have this flag if seen in the PON, and should be removed'
                                                              'because it is (1)either above min_cnt and not bkgd assessed or '
                                                              '(2)bkgd assessed and below bkgd.')

    # Background filtering option section
    parser.add_argument('--bkgd_avg', help='AVG_AF value at and under which we will assess whether variant rises above background. '
                                           'If specified, bkgd_std must also be specified.')
    parser.add_argument('--bkgd_std', help='STDEV_AF value at and under which we will assess whether variant rises above background. '
                                           'If specified, bkgd_std must also be specified.')
    parser.add_argument('--bkgd_min_cnt', default=0, type=int, help='SEEN value at and under which we will assess whether variant rises above background. '
                                                                    'If specified, bkgd_std must also be specified.')
    parser.add_argument('--pon_flag_above_bkgd', help='A variant will have this flag if below filtering thresholds specified by'
                                                      'bkgd_avg, bkgd_std, bkgd_min_cnt and VAF is above the calculated bkgd threshold.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    if args.bkgd_avg and not args.bkgd_std:
        raise Exception('bkgd_avg and bkgd_std must be specified together.')
    elif not args.bkgd_avg and args.bkgd_std:
        raise Exception('bkgd_avg and bkgd_std must be specified together.')

    return args


def main(args_input=None):
    if args_input is None:
        args_input = sys.argv[1:]
    args = get_args(args_input)

    input_vcf = VcfRecords(args.input_vcf)
    pon_vcf = VcfRecords(args.pon_vcf)

    header = input_vcf.merge_header(pon_vcf.header)
    header = input_vcf.update_header(header, args.pon_flag, "A variant will have this flag if seen in the PON but is below min_cnt and is not bkgd accessed.")
    header = input_vcf.update_header(header, args.pon_flag_filtered, "A variant will have this flag if seen in the PON, and should be removed because it is (1)either above min_cnt and not bkgd assessed or (2)bkgd assessed and below bkgd.")
    header = input_vcf.update_header(header, args.pon_flag_above_bkgd, "A variant will have this flag if below filtering thresholds specified by bkgd_avg, bkgd_std, bkgd_min_cnt and VAF is above the calculated bkgd threshold.")

    writer = vcfpy.Writer.from_path(args.output_vcf, header=header)

    for record in input_vcf.records:
        record = VcfRecord(input_vcf.records[record])
        if record.id in pon_vcf.records:
            record = record.merge(pon_vcf.records[record.id])
            record = PonRecord(record, args.min_cnt, args.bkgd_avg, args.bkgd_std, args.bkgd_min_cnt)
            # these are flagged to be filtered out with pon_flag_filtered label
            # NOTE: might not be filtered out after bkgd assesment
            # if below bkgd assessment filtering thresholds and found to be above bkgd so labelled as pon_flag_above_bkgd
            if record.above_min_cnt:
                # TODO: add clinvar or cosmic exclusion or use vcf_filter to save variants?
                record.record.add_filter(args.pon_flag_filtered)
            else:
                # these are below above_min_cnt so flagged as just pon_flag
                record.record.add_filter(args.pon_flag)

            # TODO: is this needed with BKGD filtering?
            # Bkgd assessment
            if record.below_bkgd_thresholds:
                if record.above_bkgd(input_vcf.sample):
                    # these are above_min_cnt so already labelled as pon_flag_filtered, remove pon_flag_filtered label
                    if record.above_min_cnt:
                        record.record.FILTER.remove(args.pon_flag_filtered)
                    record.record.add_filter(args.pon_flag_above_bkgd)
                else:
                    record.record.add_filter(args.pon_flag_filtered)
        writer.write_record(record.record)


if __name__ == "__main__":
    main()
