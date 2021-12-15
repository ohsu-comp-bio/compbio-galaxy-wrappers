
import vcfpy
from operator import attrgetter
from natsort import natsorted
import json

VERSION = '0.0.1'


class Records:
    def __init__(self, input_vcf):
        self.reader = vcfpy.Reader.from_path(input_vcf)

    @staticmethod
    def _get_dict_records(reader):
        records = {}
        for record in reader:
            rec_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
            if rec_id not in records:
                records[rec_id] = record
        return records

    def get_filter_records(self, filt, exclusive):
        records = []
        for record in self.reader:
            if exclusive:
                if set(filt) == set(record.FILTER):
                    records.append(record)
            else:
                for f in filt:
                    if f in record.FILTER:
                        records.append(record)
            return records

    def label_records(self, resource, label, description):
        """
        Label records based on resource VCF
        """
        self.reader.header.add_filter_line(vcfpy.OrderedDict([('ID', label), ('Description', description)]))

        records = {}
        for record in self.reader:
            rec_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
            if rec_id in self._get_dict_records(vcfpy.Reader.from_path(resource)):
                record.add_filter(label)
                records[rec_id] = record
            else:
                records[rec_id] = record
        return records


class RecordsWriter:
    def __init__(self, records):
        self.records = records

    def as_vcf(self, output_vcf, header):
        writer = vcfpy.Writer.from_path(output_vcf, header=header)
        for record in natsorted(self.records, key=attrgetter('CHROM', 'POS')):
            writer.write_record(record)

    def as_json(self, output_json):
        with open(output_json, 'w') as oj:
            to_dump = {k: v for k, v in self.records.items()}
            oj.write(json.dumps(to_dump))


# TODO: multi-sample instead of first vcf samples
class VcfMerger:
    def __init__(self, input_vcfs, callers, vcf_format="VCFv4.2"):
        self.readers = [vcfpy.Reader.from_path(vcf) for vcf in input_vcfs]
        self.callers = callers
        self.samples = list(set([name for reader in self.readers for name in reader.header.samples.names]))
        self.header = vcfpy.Header(samples=[reader.header.samples for reader in self.readers][0])
        self.header.add_line(vcfpy.HeaderLine('fileformat', vcf_format))

    def merge_headers(self):
        # using contig header from first sample only
        for contig in self.readers[0].header.get_lines('contig'):
            self.header.add_contig_line(vcfpy.OrderedDict([('ID', contig.id), ('length', contig.length)]))
        # add each caller to header as FILTER
        for reader, caller in zip(self.readers, self.callers):
            self.header.add_filter_line(vcfpy.OrderedDict([('ID', caller), ('Description', 'Variant caller label.')]))
        # add FILTER field info from each vcf to new header
        seen = []
        for reader, caller in zip(self.readers, self.callers):
            for filt in reader.header.get_lines('FILTER'):
                if filt.id not in seen:
                    seen.append(filt.id)
                    if filt.id not in ['.', 'PASS']:
                        self.header.add_filter_line(vcfpy.OrderedDict([('ID', filt.id),
                                                                       ('Description', filt.description)]))
        # add INFO field info from each vcf to new header
        for reader, caller in zip(self.readers, self.callers):
            for info in reader.header.get_lines('INFO'):
                self.header.add_info_line(vcfpy.OrderedDict([('ID', '{}_{}'.format(caller, info.id)),
                                                             ('Number', info.number), ('Type', info.type),
                                                             ('Description', '{} {}'.format(caller, info.description))]))
        # add FORMAT field info from each vcf to new header
        for reader, caller in zip(self.readers, self.callers):
            for fmt in reader.header.get_lines('FORMAT'):
                self.header.add_format_line(vcfpy.OrderedDict([('ID', '{}_{}'.format(caller, fmt.id)),
                                                               ('Number', fmt.number), ('Type', fmt.type),
                                                               ('Description', '{} {}'.format(caller, fmt.description))]))
        return self.header

    def merge_records(self):
        header = self.merge_headers()
        records = {}
        for reader, caller in zip(self.readers, self.callers):
            for record in reader:
                rec_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
                record.add_filter(caller)
                for info in header.info_ids():
                    record.INFO[info] = None
                    prefix, old_key = info.split('_', 1)
                    if prefix == caller:
                        record.INFO[info] = record.INFO.get(old_key, None)
                        record.INFO.pop(old_key, None)
                # for INFO and FORMAT ids, add caller as prefix, e.g m2_AD, fb_AD
                for sample in reader.header.samples.names:
                    for fmt in header.format_ids():
                        record.add_format(key=fmt, value=None)
                        prefix, old_key = fmt.split('_', 1)
                        if prefix == caller:
                            record.call_for_sample[sample].data[fmt] = record.call_for_sample[sample].data.get(old_key, None)
                            if old_key in record.FORMAT:
                                record.FORMAT.remove(old_key)
                                record.call_for_sample[sample].data.pop(old_key, None)
                # Uniquify records and for records called by multiple callers, append all callers to FILTER column
                if rec_id not in records:
                    records[rec_id] = record
                else:
                    seen = records[rec_id]
                    # typeError might occur if
                    for filt in record.FILTER:
                        try:
                            seen.add_filter(filt)
                        except TypeError:
                            pass
                    # for INFO and FORMAT ids, add caller as prefix, e.g m2_AD, fb_AD
                    for info in record.INFO.items():
                        prefix, old_key = info[0].split('_', 1)
                        if prefix == caller:
                            seen.INFO[info[0]] = info[1]
                    for sample in reader.header.samples.names:
                        for key in record.call_for_sample[sample].data.items():
                            prefix, old_key = key[0].split('_', 1)
                            if prefix == caller:
                                seen.call_for_sample[sample].data[key[0]] = key[1]
                    records[rec_id] = seen
            reader.close()
        return records


class VcfSelecter:
    def __init__(self, input_vcf, info_fields, format_fields, caller_priority, vcf_format="VCFv4.2"):
        self.reader = vcfpy.Reader.from_path(input_vcf)
        self.info_fields = info_fields
        self.format_fields = format_fields
        self.caller_priority = caller_priority
        self.samples = self.reader.header.samples.names
        self.header = vcfpy.Header(samples=[self.samples])
        self.header.add_line(vcfpy.HeaderLine('fileformat', vcf_format))

    def select_headers(self):
        for contig in self.reader.header.get_lines('contig'):
            self.header.add_contig_line(vcfpy.OrderedDict([('ID', contig.id), ('length', contig.length)]))
        for filt in self.reader.header.get_lines('FILTER'):
            self.header.add_filter_line(vcfpy.OrderedDict([('ID', filt.id), ('Description', filt.description)]))
        for info_field in self.info_fields:
            for caller in self.caller_priority:
                info_id = '{}_{}'.format(caller, info_field)
                if info_id in self.reader.header.info_ids():
                    info = self.reader.header.get_info_field_info(info_id)
                    self.header.add_info_line(vcfpy.OrderedDict([('ID', info.id.split('_', 1)[1]),
                                                                 ('Number', info.number),
                                                                 ('Type', info.type),
                                                                 ('Description', info.description.split(' ', 1)[1])]))
                    break
                else:
                    print('{} not found for {} in INFO column.'.format(info_field, caller))
        not_found = list(set(list(self.info_fields)) - set(self.header.info_ids()))
        if len(not_found) > 0:
            raise Exception(', '.join(not_found) + ' INFO field(s) not found in VCF')
        for format_field in self.format_fields:
            for caller in self.caller_priority:
                fmt_id = '{}_{}'.format(caller, format_field)
                if fmt_id in self.reader.header.format_ids():
                    fmt = self.reader.header.get_format_field_info(fmt_id)
                    self.header.add_format_line(vcfpy.OrderedDict([('ID', fmt.id.split('_', 1)[1]),
                                                                   ('Number', fmt.number),
                                                                   ('Type', fmt.type),
                                                                   ('Description', fmt.description.split(' ', 1)[1])]))
                    break
                else:
                    print('{} not found for {} in sample column.'.format(format_field, caller))
        not_found = list(set(list(self.format_fields)) - set(self.header.format_ids()))
        if len(not_found) > 0:
            raise Exception(', '.join(not_found) + ' FORMAT field(s) not found in VCF')
        return self.header

    def select_records(self):
        # TODO: dup?
        header = self.select_headers()
        records = {}
        for record in self.reader:
            rec_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
            for info in self.info_fields:
                for caller in self.caller_priority:
                    value = record.INFO.get('{}_{}'.format(caller, info), None)
                    if value is not None:
                        record.INFO[info] = value
                        break
            # TODO:
            # will come up if an INFO field on interest was not in any of the callers
            for info in list(record.INFO):
                if info not in header.info_ids():
                    record.INFO.pop(info, None)
            for sample in self.samples:
                for fmt in self.format_fields():
                    record.add_format(key=fmt, value=None)
                    for caller in self.caller_priority:
                        caller_key = '{}_{}'.format(caller, fmt)
                        value = record.call_for_sample[sample].data.get(caller_key, None)
                        if value is not None:
                            record.call_for_sample[sample].data[fmt] = value
                            break
                # TODO:
                # will come up if an INFO field on interest was not in any of the callers
                for call in record.calls:
                    for key in list(call.data):
                        if key not in header.format_ids():
                            call.data.pop(key, None)
                            if key in record.FORMAT:
                                record.FORMAT.remove(key)
            records[rec_id] = record
        return records
