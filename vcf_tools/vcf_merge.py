#!/usr/bin/env python

"""
Merge VCFs produced by various variant callers.

Usage:
<<<<<<< HEAD
    vcf_merge.py  --input_vcfs input1.vcf input2.vcf --caller_labels m2 fb --output_vcf output.vcf
=======
    python vcf_merge.py
        --input_vcf input1.vcf --caller_label m2
        --input_vcf input2.vcf --caller_label fb
        --output_vcf output.vcf
>>>>>>> 6e72adcd42eea3d7770a11504877ed90fb589e10

Details:
Given multiple VCFs produced by various variant callers, this will merge all vcfs into one.
For variants called by multiple callers, it will output a single record and the filter, info column and sample column
will contain those of all callers that made the call.
"""

import argparse
import vcfpy
from operator import attrgetter
from natsort import natsorted

VERSION = '2.0.0'


class VcfMerger:
    """
    Merge headers and variants in VCFs produced by different variant callers
    """

    def __init__(self, input_vcfs, callers, vcf_format="VCFv4.2"):
        self.readers = [vcfpy.Reader.from_path(vcf) for vcf in input_vcfs]
        self.callers = callers
        self.vcf_format = vcf_format
        self.samples = list(set([name for reader in self.readers for name in reader.header.samples.names]))
        self.new_filter_ids = list(
            set(['{}_{}'.format(caller, flt) for reader, caller in zip(self.readers, self.callers) for flt in
                 reader.header.filter_ids()]))
        self.new_info_ids = list(
            set(['{}_{}'.format(caller, info) for reader, caller in zip(self.readers, self.callers) for info in
                 reader.header.info_ids()]))
        self.new_format_ids = list(
            set(['{}_{}'.format(caller, fmt) for reader, caller in zip(self.readers, self.callers) for fmt in
                 reader.header.format_ids()]))

    @staticmethod
    def add_header(vfield, caller, field, ids):
        if vfield in ['INFO', 'FORMAT']:
            field = '{}_{}'.format(caller, field)
        if field not in ['.', 'PASS']:
            if field not in ids:
                return True
    
    @staticmethod
    def update_header(caller, field):
        field_odict = vcfpy.OrderedDict()
        for k, v in field.mapping.items():
            if k == 'ID':
                v = '{}_{}'.format(caller, v)
            if k == 'Description':
                v = '{} {}'.format(caller, v)
            field_odict[k] = v
        # k, v = field_odict['ID'], '<{}>'.format(','.join(str(k)+'='+str(v) for k, v in field_odict.items()))
        return field_odict

    def merge_headers(self):
        new_header = vcfpy.Header(samples=[reader.header.samples for reader in self.readers][0])
        new_header.add_line(vcfpy.HeaderLine('fileformat', self.vcf_format))

        new_contig_headers = {}
        for reader, caller in zip(self.readers, self.callers):
            contig_headers = {contig.id: contig for contig in reader.header.get_lines('contig')}
            for contig in contig_headers:
                if self.add_header('contig', caller, contig, new_contig_headers):
                    new_contig_headers[contig_headers[contig].id] = contig_headers[contig]
                    new_header.add_contig_line(contig_headers[contig].mapping)
            new_header.add_filter_line(vcfpy.OrderedDict([('ID', caller), ('Description', 'Variant caller label.')]))
            filter_headers = {flt.id: flt for flt in reader.header.get_lines('FILTER')}
            for flt in filter_headers:
                if self.add_header('FILTER', caller, flt, new_header.filter_ids()):
                    new_header.add_filter_line(filter_headers[flt].mapping)
            info_headers = {info.id: info for info in reader.header.get_lines('INFO')}
            for info in info_headers:
                if self.add_header('INFO', caller, info, new_header.info_ids()):
                    new_header.add_info_line(self.update_header(caller, info_headers[info]))
            format_headers = {fmt.id: fmt for fmt in reader.header.get_lines('FORMAT')}
            for fmt in format_headers:
                if self.add_header('FORMAT', caller, fmt, new_header.format_ids()):
                    new_header.add_format_line(self.update_header(caller, format_headers[fmt]))
        return new_header

    @staticmethod
    def update_fields(caller, field_dict, rec_field_dict):
        for field in field_dict.keys():
            c, k = field.split('_', 1)
            if c == caller:
                field_dict[field] = rec_field_dict.get(k, None)
        return field_dict

    def update_record(self, caller, new_record, record, samples):
        new_record.add_filter(caller)
        new_record.INFO = self.update_fields(caller, new_record.INFO, record.INFO)
        for sample in samples:
            for fmt in new_record.FORMAT.keys():
                new_record.add_format(fmt)
            new_record.call_for_sample[sample].data = self.update_fields(caller, new_record.call_for_sample[sample].data, record.call_for_sample[sample].data)
        return new_record

    def merge_records(self):
        records = {}
        for reader, caller in zip(self.readers, self.callers):
            for record in reader:
                new_record = vcfpy.Record(CHROM=record.CHROM,
                                          POS=record.POS,
                                          ID=record.ID,
                                          REF=record.REF,
                                          ALT=record.ALT,
                                          QUAL=(None if record.QUAL == float('inf') else record.QUAL),  # happens with fc
                                          FILTER=([flt for flt in record.FILTER if flt not in ['.', 'PASS']]),
                                          INFO=dict.fromkeys(self.new_info_ids, None),
                                          FORMAT=dict.fromkeys(self.new_format_ids, None),
                                          calls=[vcfpy.Call(sample, dict.fromkeys(self.new_format_ids, None)) for sample in reader.header.samples.names])
                var_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
                if var_id not in records:
                    records[var_id] = self.update_record(caller, new_record, record, reader.header.samples.names)
                else:
                    records[var_id] = self.update_record(caller, records[var_id], record, reader.header.samples.names)
        return list(records.values())


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input_vcf', action='append', help="Input VCF")
    parser.add_argument('--caller_label', action='append', help="Label for given input vcf")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    merge = VcfMerger(args.input_vcf, args.caller_label)

    writer = vcfpy.Writer.from_path(args.output_vcf, header=merge.merge_headers())
    for record in natsorted(merge.merge_records(), key=attrgetter('CHROM', 'POS')):
        writer.write_record(record)


if __name__ == "__main__":
    main()
