#!/usr/bin/env python

"""
Merge VCFs produced by various variant callers.

Usage:
python vcf_merge.py
    --input_vcf input1.vcf --caller_label m2
    --input_vcf input2.vcf --caller_label fb
    --output_vcf output.vcf
    --caller_priority sc fb m2 fc mrdf
    --caller_priority_indel sc m2 fb fc mrdf
    --filter_fields MAsite PON PON_OV PON_FILT StrandBias strand_bias BelowBKGD
    --info_fields AF DP TLOD SRF SRR SAF SAR ITD_LENGTH
    --format_fields GT AF DP AD SB

Details:
Given multiple VCFs produced by various variant callers, this will merge all vcfs into one.
For variants called by multiple callers, it will output a merged single record.
The FILTER column will contain the labels of the input VCFs.
The INFO and FORMAT/SAMPLE column will contain either all fields or the fields passed to --info_fields and --sample_fields, respectively.
A field not in the input VCF will raise an error.
"""

import argparse
import sys
from operator import attrgetter
from natsort import natsorted
import vcfpy


VERSION = '2.2.0'


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


class VcfSelecter:
    """
    In a merged VCF, select VCF fields of interest
    """
    def __init__(self, records, header, filter_fields, info_fields, format_fields, caller_priority, caller_priority_indel, vcf_format="VCFv4.2"):
        self.records = records
        self.header = header
        self.filter_fields = filter_fields if filter_fields else self.header.filter_ids()
        self.info_fields = info_fields
        self.format_fields = format_fields
        self.caller_priority = caller_priority
        self.caller_priority_indel = caller_priority_indel
        self.vcf_format = vcf_format

    @staticmethod
    def add_header(vfield, field, fields, ids):
        if vfield in ['INFO', 'FORMAT']:
            field = field.split('_', 1)[1]
        if field in fields:
            if field not in ids:
                return True

    @staticmethod
    def update_header(field):
        field_odict = vcfpy.OrderedDict()
        for k, v in field.mapping.items():
            if k == 'ID':
                v = v.split('_', 1)[1]
            if k == 'Description':
                v = v.split(' ', 1)[1]
            field_odict[k] = v
        # k, v = field_odict['ID'], '<{}>'.format(','.join(str(k)+'='+str(v) for k, v in field_odict.items()))
        return field_odict

    def select_headers(self):
        new_header = vcfpy.Header(samples=self.header.samples)
        new_header.add_line(vcfpy.HeaderLine('fileformat', self.vcf_format))

        contig_headers = {contig.id: contig for contig in self.header.get_lines('contig')}
        for contig in contig_headers:
            new_header.add_contig_line(contig_headers[contig].mapping)
        filter_headers = {flt.id: flt for flt in self.header.get_lines('FILTER')}
        for flt in filter_headers:
            if self.add_header('FILTER', flt, self.filter_fields+self.caller_priority, new_header.filter_ids()):
                new_header.add_filter_line(filter_headers[flt].mapping)
        info_headers = {info.id: info for info in self.header.get_lines('INFO')}
        for info in info_headers:
            if self.add_header('INFO', info, self.info_fields, new_header.info_ids()):
                new_header.add_info_line(self.update_header(info_headers[info]))
        format_headers = {fmt.id: fmt for fmt in self.header.get_lines('FORMAT')}
        for fmt in format_headers:
            if self.add_header('FORMAT', fmt, self.format_fields, new_header.format_ids()):
                new_header.add_format_line(self.update_header(format_headers[fmt]))
        return new_header

    @staticmethod
    def add_field(rec_field, rec_vfield_infos):
        if rec_field in rec_vfield_infos:
            if isinstance(rec_vfield_infos, dict):
                if rec_vfield_infos[rec_field] is not None:
                    return True

    def update_fields(self, vfield, fields, caller, vfield_infos, record_vfield_infos, new_record):
        for field in fields:
            rec_field = '{}_{}'.format(caller, field)
            if field not in vfield_infos:
                if vfield == 'INFO':
                    if self.add_field(rec_field, record_vfield_infos):
                        vfield_infos[field] = record_vfield_infos[rec_field]
                else:
                    if field not in new_record.FORMAT:
                        new_record.add_format(field)
                        vfield_infos[field] = record_vfield_infos.get(rec_field, '.')

                        # TODO: AF in fb is Number = 1 but 'A' in m2. This avoids vcfpy writing m2 AF fields as []
                        # might consider changing m2 header number first
                        if rec_field == 'm2_AF':
                            vfield_infos[field] = record_vfield_infos[rec_field][0]

        return vfield_infos

    def update_record(self, new_record, record, samples):
        first_caller = new_record.FILTER[0]
        new_record.INFO = self.update_fields('INFO', self.info_fields, first_caller, new_record.INFO, record.INFO, new_record)
        for sample in samples:
            new_record.call_for_sample[sample] = vcfpy.Call(sample, vcfpy.OrderedDict())
            new_record.call_for_sample[sample].data = self.update_fields('FORMAT', self.format_fields, first_caller, new_record.call_for_sample[sample].data, record.call_for_sample[sample].data, new_record)
            new_record.calls.append(new_record.call_for_sample[sample])

    def select_records(self):
        records = {}
        for record in self.records:
            caller_priority = self.caller_priority
            if record.ALT[0].type != 'SNV':
                caller_priority = self.caller_priority_indel
            new_record = vcfpy.Record(CHROM=record.CHROM,
                                      POS=record.POS,
                                      ID=record.ID,
                                      REF=record.REF,
                                      ALT=record.ALT,
                                      QUAL=record.QUAL,
                                      FILTER=sorted([flt for flt in record.FILTER if flt in caller_priority+self.filter_fields], key=(caller_priority + record.FILTER).index),
                                      INFO={},
                                      FORMAT=[],
                                      calls=[])
            self.update_record(new_record, record, self.header.samples.names)
            var_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
            records[var_id] = new_record
        return list(records.values())


def get_args(args_input):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input_vcf', action='append', help="Input VCF")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('--caller_label', action='append', help="Label for given input vcf")

    parser.add_argument('--caller_priority', nargs="+", help="Priority of caller labels")
    parser.add_argument('--caller_priority_indel', nargs="+", help="Priority of caller labels for INDELs")
    parser.add_argument('--filter_fields', nargs="+", help="FILTER fields of interest")
    parser.add_argument('--info_fields', nargs="+", help="INFO fields of interest")
    parser.add_argument('--format_fields', nargs="+", help="FORMAT fields of interest")

    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args(args_input)
    return args


def main(args_input=sys.argv[1:]):
    args = get_args(args_input)

    merge = VcfMerger(args.input_vcf, args.caller_label)

    select = VcfSelecter(merge.merge_records(), merge.merge_headers(),
                         args.filter_fields, args.info_fields, args.format_fields,
                         args.caller_priority, args.caller_priority_indel)

    writer = vcfpy.Writer.from_path(args.output_vcf, header=select.select_headers())
    for record in natsorted(select.select_records(), key=attrgetter('CHROM', 'POS')):
        writer.write_record(record)


if __name__ == "__main__":
    main()
