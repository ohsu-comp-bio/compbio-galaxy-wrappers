#!/usr/bin/env python

"""
Filter a merged VCF based on FILTER column.

Usage:
vcf_filter.py 'input.vcf' 'output.vcf' 'output2_vcf'
    --callers "sc fb m2 fc mrdf"
    --exclude "PON_FILT" --include "mrdf fc clinvar.CLNSIG=Pathogenic clinvar.CLNSIG=Likely_pathogenic clinvar.CLNSIG=Pathogenic/Likely_pathogenic CLNSIGCONF=Pathogenic CLNSIGCONF=Likely_pathogenic" --snp_threshold 1.0 --indel_threshold 1.0
    --exclude "BelowBKGD" --include "mrdf fc clinvar.CLNSIG=Pathogenic clinvar.CLNSIG=Likely_pathogenic clinvar.CLNSIG=Pathogenic/Likely_pathogenic CLNSIGCONF=Pathogenic CLNSIGCONF=Likely_pathogenic" --snp_threshold 1.0 --indel_threshold 1.0
    --exclude "StrandBias" --include "mrdf fc clinvar.CLNSIG=Pathogenic clinvar.CLNSIG=Likely_pathogenic clinvar.CLNSIG=Pathogenic/Likely_pathogenic CLNSIGCONF=Pathogenic CLNSIGCONF=Likely_pathogenic" --snp_threshold 1.0 --indel_threshold 1.0
    --exclude "fb" --include "mrdf fc clinvar.CLNSIG=Pathogenic clinvar.CLNSIG=Likely_pathogenic clinvar.CLNSIG=Pathogenic/Likely_pathogenic CLNSIGCONF=Pathogenic CLNSIGCONF=Likely_pathogenic" --snp_threshold 0.1 --indel_threshold 0.02
    --exclude "m2"  --include "mrdf fc clinvar.CLNSIG=Pathogenic clinvar.CLNSIG=Likely_pathogenic clinvar.CLNSIG=Pathogenic/Likely_pathogenic clinvar.CLNSIGCONF=Pathogenic clinvar.CLNSIGCONF=Likely_pathogenic" --snp_threshold 0.1 --indel_threshold 0.1

Given a merged VCF, this will filter out variants that contain the exclude argument in the FILTER column.
"""

import argparse
import vcfpy

VERSION = '2.2.0'


def get_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', help='Input VCF')
    parser.add_argument('output_vcf', help='Output VCF')
    parser.add_argument('output2_vcf', help='Output VCF of filtered out variants')
    parser.add_argument('--callers', help='FILTER annotations of variant callers used in the VCF')
    parser.add_argument('--exclude',  action='append', help='Filters to use to exclude VCF records')
    parser.add_argument('--include', action='append', help='Filters to use to include VCF records')
    parser.add_argument('--inc_multicalled', action='append', help='Include multicalled variants')
    parser.add_argument('--snp_threshold', action='append', type=float, help='For SNPs, only perform filtering if VAF is below this threshold')
    parser.add_argument('--indel_threshold', action='append', type=float, help='For INDELs, only perform filtering if VAF is below this threshold')
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)

    args = parser.parse_args()

    if len(args.exclude) != len(args.snp_threshold) != len(args.snp_threshold):
        raise Exception('Must specify a snp_threshold and indel_threshold arguments for each exclude argument')
    return args


class VcfFilters:
    def __init__(self, filters, callers):
        self.filters = filters
        self.callers = callers
        self.callers = [filt for filt in filters if filt in callers]
        self.flags = [filt for filt in filters if filt not in callers]
        self.is_multicalled = True if len(self.callers) > 1 else False
        self.info = self._get_info()

    def _get_info(self):
        info_dict = {}
        for filt in self.filters:
            if '=' in filt:
                k, v = filt.split('=')
                if k not in info_dict:
                    info_dict[k] = [v]
                else:
                    info_dict[k].append(v)
        return info_dict

    def caller_prioritize(self):
        return self.callers + self.flags


class RecordKeeper:
    def __init__(self, record, callers):
        self.record = record
        self.callers = callers
        self.record_filters = VcfFilters(record.FILTER, callers)

    def check_record_infos(self, info_terms_dict):
        for d in [self.record.INFO, self.record.calls[0].data]:
            for k in info_terms_dict:
                if k in d:
                    for item in info_terms_dict[k]:
                        if any(item == item2 for item2 in d[k]) or item == '*':
                            return True

    def exclude_record(self, exclude, threshold):
        exclude = VcfFilters(exclude, self.callers)
        if self.record.calls[0].data['AF'] <= float(threshold):
            # if only given a caller check that it is the only caller present
            if exclude.callers and not (exclude.flags or exclude.info):
                if exclude.callers == self.record_filters.callers:
                    return True
            # if given only flag(s) (or info field(s))
            elif not exclude.callers and (exclude.flags or exclude.info):
                if any(filt in self.record_filters.filters for filt in exclude.flags) or self.check_record_infos(exclude.info):
                    return True
            # if given callers and flags or info fields
            elif exclude.callers and (exclude.flags or exclude.info):
                if exclude.callers == self.record_filters.callers and (any(filt in self.record_filters.filters for filt in exclude.flags) or self.check_record_infos(exclude.info)):
                    return True
            else:
                raise Exception('Filters must be callers, flags or info (given as ID=value)')
        return False

    def include_record(self, include, inc_multicalled):
        include = VcfFilters(include, self.callers)
        if any(filt in self.record_filters.filters for filt in include.filters) or self.check_record_infos(include.info) or (inc_multicalled == "True" and self.record_filters.is_multicalled):
            return True
        return False

    def keep_record(self, exclude, include, inc_multicalled, threshold):
        if exclude and not include:
            return not self.exclude_record(exclude, threshold)
        elif include and not exclude:
            return self.include_record(include, inc_multicalled)
        elif exclude and include:
            if self.include_record(include, inc_multicalled):
                return True
            else:
                return not self.exclude_record(exclude, threshold)
        else:
            raise Exception('Filters must be callers, flags or info (given as ID=value)')


def main():
    args = get_args()

    reader = vcfpy.Reader.from_path(args.input_vcf)
    header = reader.header
    writer = vcfpy.Writer.from_path(args.output_vcf, header=header)
    writer2 = vcfpy.Writer.from_path(args.output2_vcf, header=header)

    callers = args.callers.split(' ')

    for record in reader:
        keep = True
        if args.exclude:
            for i in range(len(args.exclude)):
                exclude = args.exclude[i].split(' ')
                try:
                    include = args.include[i].split(' ')
                except IndexError:
                    include = []
                try:
                    inc_multicalled = args.inc_multicalled[i]
                except IndexError:
                    inc_multicalled = False

                if keep:  # if a record has not been excluded by a preceding exclude argument
                    recordkeeper = RecordKeeper(record, callers)
                    if record.ALT[0].type == 'SNV':
                        keep = recordkeeper.keep_record(exclude, include, inc_multicalled, float(args.snp_threshold[i]))
                    else:
                        keep = recordkeeper.keep_record(exclude, include, inc_multicalled, float(args.indel_threshold[i]))
                else:
                    break

        record.FILTER = VcfFilters(record.FILTER, callers).caller_prioritize()
        if keep:
            writer.write_record(record)
        else:
            writer2.write_record(record)


if __name__ == "__main__":
    main()
