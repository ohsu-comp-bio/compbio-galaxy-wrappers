#!/usr/bin/env python

"""
Select fields of interest from a merged VCF.

Example usage: vcf_select.py --input_vcf 'file.vcf' --info_fields AF DP TLOD --sample_fields AF AD DP PUMI
--caller_priority m2 fb --output_vcf 'output.vcf'

Details: Given a merged VCF (VCF with calls from various variant callers), this will create a new vcf
containing info fields and sample fields passed to --info_fields and --sample_fields, respectively.
A field not in the input VCF will raise an error.
"""

import argparse
import vcfpy
import warnings
from operator import attrgetter
from natsort import natsorted

VERSION = '0.0.1'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(
        description="vcf_select.py --input_vcf 'file.vcf' "
                    "--info_fields AF DP TLOD "
                    "--sample_fields AF AD DP PUMI"
                    "--caller_priority label1 label2 labelN "
                    "--output_vcf 'output.vcf'")
    parser.add_argument('--input_vcf',
                        help="vcf file from which to select info.")
    parser.add_argument('--info_fields', nargs="+",
                        help="select INFO fields of interest.")
    parser.add_argument('--sample_fields', nargs="+",
                        help="select sample fields of interest.")
    parser.add_argument('--caller_priority', nargs="+",
                        help="set the priority of variant caller.")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()
    return args


class VcfSelectFields(object):

    def __init__(self, input_vcf, info_fields, sample_fields, caller_priority, output_vcf):
        self.reader = vcfpy.Reader.from_path(input_vcf)
        self.info_fields = info_fields
        self.sample_fields = sample_fields
        self.caller_priority = caller_priority
        self.write_header = vcfpy.Header(samples=self.reader.header.samples)
        self.select_contig_header()
        self.select_filter_header()
        self.select_info_header()
        self.select_format_header()
        self.records = self.select_record_fields()
        self.output_vcf = output_vcf
        self.write_merged(self.records)

    def select_contig_header(self):
        for contig in self.reader.header.get_lines('contig'):
            self.write_header.add_contig_line(
                vcfpy.OrderedDict(
                    [('ID', contig.id),
                     ('length', contig.length)]
                )
            )

    def select_filter_header(self):
        for filter in self.reader.header.get_lines('FILTER'):
            self.write_header.add_filter_line(
                vcfpy.OrderedDict(
                    [('ID', filter.id),
                     ('Description', filter.description)]
                )
            )

    def select_info_header(self):
        for info_field in self.info_fields:
            for caller in self.caller_priority:
                id = '{}_{}'.format(caller, info_field)
                if id in self.reader.header.info_ids():
                    info = self.reader.header.get_info_field_info(id)
                    self.write_header.add_info_line(
                        vcfpy.OrderedDict(
                            [('ID', info.id.split('_', 1)[1]),
                             ('Number', info.number),
                             ('Type', info.type),
                             ('Description', info.description.split(' ', 1)[1])]
                        )
                    )
                    break
                else:
                    print('{} not found for {} in INFO column.'.format(info_field, caller))
        not_found = list(set(list(self.info_fields)) - set(self.write_header.info_ids()))
        if len(not_found) > 0:
            raise Exception(', '.join(not_found) + ' INFO field(s) not found in VCF')



    def select_format_header(self):
        for sample_field in self.sample_fields:
            for caller in self.caller_priority:
                id = '{}_{}'.format(caller, sample_field)
                if id in self.reader.header.info_ids():
                    format = self.reader.header.get_info_field_info('{}_{}'.format(caller, sample_field))
                    self.write_header.add_format_line(
                        vcfpy.OrderedDict(
                            [('ID', format.id.split('_', 1)[1]),
                             ('Number', format.number),
                             ('Type', format.type),
                             ('Description', format.description.split(' ', 1)[1])]
                        )
                    )
                    break
                else:
                    print('{} not found for {} in sample column.'.format(sample_field, caller))
        not_found = list(set(list(self.sample_fields)) - set(self.write_header.info_ids()))
        if len(not_found) > 0:
            raise Exception(', '.join(not_found) + ' sample field(s) not found in VCF')


    def select_record_fields(self):
        records = []
        for record in self.reader:

            for info in self.write_header.info_ids():
                for caller in self.caller_priority:
                    value = record.INFO.get('{}_{}'.format(caller, info), None)
                    if value is not None:
                        record.INFO[info] = value
                        break
            for info in list(record.INFO):
                if info not in self.write_header.info_ids():
                    record.INFO.pop(info, None)

            for sample in self.reader.header.samples.names:
                for format in self.write_header.format_ids():
                    record.add_format(key=format, value=None)
                    for caller in self.caller_priority:
                        caller_key = '{}_{}'.format(caller, format)
                        value = record.call_for_sample[sample].data.get(caller_key, None)
                        if value is not None:
                            record.call_for_sample[sample].data[format] = value
                            break
                for call in record.calls:
                    for key in list(call.data):
                        if key not in self.write_header.format_ids():
                            call.data.pop(key, None)
                            if key in record.FORMAT:
                                record.FORMAT.remove(key)

            records.append(record)
        return records

    def write_merged(self, records):
        writer = vcfpy.Writer.from_path(self.output_vcf, header=self.write_header)
        for record in natsorted(records, key=attrgetter('CHROM', 'POS')):
            writer.write_record(record)


def main():
    args = supply_args()

    VcfSelectFields(args.input_vcf, args.info_fields, args.sample_fields, args.caller_priority, args.output_vcf)


if __name__ == "__main__":
    main()
