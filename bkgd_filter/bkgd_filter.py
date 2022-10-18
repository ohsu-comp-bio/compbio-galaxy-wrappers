#!/usr/bin/env python

"""
Label variants found to be below the background threshold.
"""

import argparse
import vcfpy

VERSION = '0.0.1'


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('input_vcf', help="Input VCF")
    parser.add_argument('output_vcf', help="Output VCF.")
    parser.add_argument('filter_id', help="FILTER ID to add to VCF record")
    parser.add_argument('filter_desc', help="FILTER Description to add to VCF header")
    parser.add_argument('--callers', help="Callers in a merged VCF.")
    parser.add_argument('--fc_label', help="Forced calls label.")
    parser.add_argument('--filter_on', help="VCF FORMAT field ID to filter on")
    parser.add_argument('--filter_from', help="VCF INFO field ID to filter on. If not given, threshold argument must be given")
    parser.add_argument('--threshold', help="Threshold at which filter_on must be under to perform filtering."
                                            "If filter_from is not given, threshold to filter filter_on on.")
    parser.add_argument('--filter_on_alt_allele', action='store_true', help="Select filter_on field of ALT allele e.g BE_A or BE_G")
    parser.add_argument('--compare_on_alt_allele', action='store_true', help="Compare filter_on and filter_from of ALT allele e.g BE_A or BE_G")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)

    args = parser.parse_args()
    return args


def get_filtering_val(record, field_type, field, use_alt_allele):
    if field_type == 'FORMAT' or field_type == 'SAMPLE':
        field_dict = record.calls[0].data
    elif field_type == 'INFO':
        field_dict = record.INFO
    else:
        raise Exception('field_type must either be INFO, FORMAT, SAMPLE')

    if use_alt_allele:
        field = field+'_'+record.ALT[0].value

    try:
        field_val = field_dict[field]
    except KeyError:
        print('{}:{} {} field not annotated.'.format(record.CHROM, record.POS, field))
        field_val = 0
    return field_val


def main():
    args = get_args()

    reader = vcfpy.Reader.from_path(args.input_vcf)
    header = reader.header
    header.add_filter_line(vcfpy.OrderedDict([('ID', args.filter_id), ('Description', args.filter_desc)]))
    writer = vcfpy.Writer.from_path(args.output_vcf, header=header)

    for record in reader:
        filt_on_val = get_filtering_val(record, 'FORMAT', args.filter_on, args.filter_on_alt_allele)
        blw_threshold = record.calls[0].data[args.filter_on] < float(args.threshold)
        if args.filter_from and args.threshold != 0:
            # filter filter_on based on filter_from
            filt_from_val = get_filtering_val(record, 'INFO', args.filter_from, args.compare_on_alt_allele)
        elif not args.filter_from and args.threshold != 0:
            # filter filter_on based on given threshold
            filt_from_val = args.threshold
        else:
            raise Exception('Missing required arguments')  # TODO

        if record.ALT[0].type == 'SNV':
            if blw_threshold and filt_on_val < filt_from_val:
                record.add_filter(args.filter_id)
        else:
            # TODO: for now, filtering fc only indels
            # if indel is fc only
            if args.fc_label in record.FILTER:
                if not any(caller in record.FILTER for caller in [caller for caller in args.callers.split(' ') if caller != args.fc_label]):
                    record.add_filter(args.filter_id)
            else:
                # if indel has no depth
                if record.calls[0].data['AD'][1] == 0:
                    record.add_filter(args.filter_id)
        writer.write_record(record)


if __name__ == "__main__":
    main()
