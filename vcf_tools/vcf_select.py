#!/usr/bin/env python

"""
Select fields of interest from a merged VCF.

Usage:
    vcf_select.py 'input.vcf' 'output.vcf'
        --filter_fields MAsite PON PON_OV PON_FILT StrandBias strand_bias BelowBKGD
        --info_fields AF DP TLOD SRF SRR SAF SAR ITD_LENGTH
        --format_fields GT AF DP AD SB
        --caller_priority sc fb m2 fc mrdf
        --caller_priority_indel sc m2 fb fc mrdf

Details:
Given a merged VCF (VCF with calls from various variant callers), this will create a new vcf
containing the filter info and sample fields passed to --info_fields and --sample_fields, respectively.
A field not in the input VCF will raise an error.
"""

import argparse
import vcfpy
from operator import attrgetter
from natsort import natsorted

VERSION = '2.0.0'


class VcfSelecter:
    """
    In a merged VCF, select VCF fields of interest
    """
    def __init__(self, input_vcf, filter_fields, info_fields, format_fields, caller_priority, caller_priority_indel, vcf_format="VCFv4.2"):
        self.reader = vcfpy.Reader.from_path(input_vcf)
        self.filter_fields = filter_fields
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
        new_header = vcfpy.Header(samples=self.reader.header.samples)
        new_header.add_line(vcfpy.HeaderLine('fileformat', self.vcf_format))

        contig_headers = {contig.id: contig for contig in self.reader.header.get_lines('contig')}
        for contig in contig_headers:
            new_header.add_contig_line(contig_headers[contig].mapping)
        filter_headers = {flt.id: flt for flt in self.reader.header.get_lines('FILTER')}
        for flt in filter_headers:
            if self.add_header('FILTER', flt, self.filter_fields+self.caller_priority, new_header.filter_ids()):
                new_header.add_filter_line(filter_headers[flt].mapping)
        info_headers = {info.id: info for info in self.reader.header.get_lines('INFO')}
        for info in info_headers:
            if self.add_header('INFO', info, self.info_fields, new_header.info_ids()):
                new_header.add_info_line(self.update_header(info_headers[info]))
        format_headers = {fmt.id: fmt for fmt in self.reader.header.get_lines('FORMAT')}
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
                if self.add_field(rec_field, record_vfield_infos):
                    if vfield == 'INFO':
                        vfield_infos[field] = record_vfield_infos[rec_field]
                    else:
                        if field not in new_record.FORMAT:
                            new_record.add_format(field)
                            vfield_infos[field] = record_vfield_infos[rec_field]

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
        for record in self.reader:
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
            self.update_record(new_record, record, self.reader.header.samples.names)
            var_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
            records[var_id] = new_record
        return list(records.values())


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('input_vcf', help="Input VCF")
    parser.add_argument('output_vcf', help="Output VCF.")
    parser.add_argument('--filter_fields', nargs="+", help="Filter fields of interest")
    parser.add_argument('--info_fields', nargs="+", help="Info fields of interest")
    parser.add_argument('--format_fields', nargs="+", help="Format fields of interest")
    parser.add_argument('--caller_priority', nargs="+", help="Priority of caller labels")
    parser.add_argument('--caller_priority_indel', nargs="+", help="Priority of caller labels for INDELs")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    select = VcfSelecter(args.input_vcf, args.filter_fields, args.info_fields, args.format_fields, args.caller_priority, args.caller_priority_indel)

    writer = vcfpy.Writer.from_path(args.output_vcf, header=select.select_headers())
    for record in natsorted(select.select_records(), key=attrgetter('CHROM', 'POS')):
        writer.write_record(record)


if __name__ == "__main__":
    main()
