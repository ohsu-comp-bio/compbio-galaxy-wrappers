#!/usr/bin/env python

"""
Merge VCFs produced by various variant callers.

Example usage: vcf_merge.py  --input_vcfs file1 file2 --caller_labels m2 fb --output_vcf output.vcf

Details: Given multiple vcfs produced by various variant callers, this will merge all vcfs into one.
For variants called by multiple callers, it will output a single record and the info column and sample column
will contain those of all callers that made the call.
"""

import argparse
import vcfpy
from operator import attrgetter
from natsort import natsorted

VERSION = '0.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(
        description="vcf_merge.py --input_vcfs 'file1.vcf' 'file2.vcf' 'fileN.vcf' --caller_labels label1 label2 labelN --output_vcf 'output.vcf'")
    parser.add_argument('--input_vcfs', nargs="+",
                        help="vcf files from multiple variant callers")
    parser.add_argument('--caller_labels', nargs="+",
                        help="Labels for the variant callers used to produce the input vcfs.")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()
    return args


class Merger(object):
    """
    :param object:
    :return:
    """

    def __init__(self, readers, callers, output_vcf):
        self.readers = readers
        self.callers = callers
        self.samples = list(set([name for reader in self.readers for name in reader.header.samples.names]))
        # TODO: multi-sample? using first vcf samples here
        self.merge_header = vcfpy.Header(samples=[reader.header.samples for reader in self.readers][0])
        self.add_file_format()
        self.merge_contig_header()
        self.add_caller_filter_header()
        self.merge_filter_header()
        self.merge_info_header()
        self.merge_format_header()
        self.records = self.merge_records()
        self.output_vcf = output_vcf
        self.write_merged(self.records)

    def add_file_format(self):
        self.merge_header.add_line(vcfpy.HeaderLine('fileformat', 'VCFv4.2'))

    # TODO: maybe use DuplicateHeaderLineWarning instead
    def merge_contig_header(self):
        for contig in self.readers[0].header.get_lines('contig'):
            self.merge_header.add_contig_line(
                vcfpy.OrderedDict(
                    [('ID', contig.id),
                     ('length', contig.length)]
                )
            )

    def add_caller_filter_header(self):
        """
        Add caller info as FILTER to writer header
        """
        for reader, caller in zip(self.readers, self.callers):
            reader.header.add_filter_line(
                vcfpy.OrderedDict(
                    [('ID', caller),
                     ('Description', 'Variant caller label.')])
            )

    def merge_filter_header(self):
        seen = []
        for reader, caller in zip(self.readers, self.callers):
            for filter in reader.header.get_lines('FILTER'):
                if filter.id not in seen:
                    seen.append(filter.id)
                    if filter.id not in ['.', 'PASS']:
                        self.merge_header.add_filter_line(
                            vcfpy.OrderedDict(
                                [('ID', filter.id),
                                 ('Description', filter.description)]
                            )
                        )


    def merge_info_header(self):
        for reader, caller in zip(self.readers, self.callers):
            for info in reader.header.get_lines('INFO'):
                self.merge_header.add_info_line(
                    vcfpy.OrderedDict(
                        [('ID', '{}_{}'.format(caller, info.id)),
                         ('Number', info.number),
                         ('Type', info.type),
                         ('Description', '{} {}'.format(caller, info.description))]
                    )
                )

    def merge_format_header(self):
        for reader, caller in zip(self.readers, self.callers):
            for format in reader.header.get_lines('FORMAT'):
                self.merge_header.add_format_line(
                    vcfpy.OrderedDict(
                        [('ID', '{}_{}'.format(caller, format.id)),
                         ('Number', format.number),
                         ('Type', format.type),
                         ('Description', '{} {}'.format(caller, format.description))]
                    )
                )

    def merge_records(self):
        records = {}
        for reader, caller in zip(self.readers, self.callers):
            for record in reader:
                record.add_filter(caller)
                for info in self.merge_header.info_ids():
                    record.INFO[info] = None
                    prefix, old_key = info.split('_', 1)
                    if prefix == caller:
                        record.INFO[info] = record.INFO.get(old_key, None)
                        record.INFO.pop(old_key, None)
                for sample in reader.header.samples.names:
                    for format in self.merge_header.format_ids():
                        record.add_format(key=format, value=None)
                        prefix, old_key = format.split('_', 1)
                        if prefix == caller:
                            record.call_for_sample[sample].data[format] = record.call_for_sample[sample].data.get(
                                old_key, None)
                            if old_key in record.FORMAT:
                                record.FORMAT.remove(old_key)
                                record.call_for_sample[sample].data.pop(old_key, None)

                # Uniquify
                var_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
                if var_id not in records:
                    records[var_id] = record
                else:
                    seen_record = records[var_id]
                    for filter in record.FILTER:
                        try:
                            seen_record.add_filter(filter)
                        except TypeError:
                            pass
                    for info in record.INFO.items():
                        prefix, old_key = info[0].split('_', 1)
                        if prefix == caller:
                            seen_record.INFO[info[0]] = info[1]
                    for sample in reader.header.samples.names:
                        for key in record.call_for_sample[sample].data.items():
                            prefix, old_key = key[0].split('_', 1)
                            if prefix == caller:
                                seen_record.call_for_sample[sample].data[key[0]] = key[1]
                    # update dictionary
                    records[var_id] = seen_record

            reader.close()
        return records

    def write_merged(self, records):
        writer = vcfpy.Writer.from_path(self.output_vcf, header=self.merge_header)
        for record in natsorted(records.values(), key=attrgetter('CHROM', 'POS')):
            writer.write_record(record)


def main():
    args = supply_args()

    readers = [vcfpy.Reader.from_path(vcf) for vcf in args.input_vcfs]

    Merger(readers, args.caller_labels, args.output_vcf)


if __name__ == "__main__":
    main()
