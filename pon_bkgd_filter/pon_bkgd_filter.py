#!/usr/bin/env python

"""
Filter variants in a VCF based on PON background allele frequencies

Example:
python pon_bkgd_filter/pon_bkgd_filter.py "input.vcf" "output.vcf"
    "BKGD" "Below Background"
    "fc_bkgd_metrics.txt" "fc"
    --per_base_bkgd_metrics "per_base_bkgd_metrics.txt" --vaf_threshold 0.1

1.1.0 - make both forced call and per bases metrics inputs optional
"""

import argparse
import vcfpy


VERSION = '1.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', help='Input VCF')
    parser.add_argument('output_vcf', help='Output VCF')
    parser.add_argument('filter_id', help='FILTER ID to add to VCF record')

    parser.add_argument('--fc_bkgd_metrics', help='Text file containing freq stats for forced call')
    parser.add_argument('--fc_label', help='Forced call label')

    parser.add_argument('--per_base_bkgd_metrics', help='Text file containing freq stats of locus in panel')

    parser.add_argument('--vaf_threshold', help='VAF threshold under to perform filtering')
    parser.add_argument('--alt_specific', action='store_true', help='Use ALT allele threshold to filter')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def add_metrics_header(field, prefix, metric, field_number, field_type, field_desc, header):
    if field == 'INFO':
        header.add_info_line(vcfpy.OrderedDict([('ID', '{}{}'.format(prefix, metric)),
                                                ('Number', field_number), ('Type', field_type),
                                                ('Description', field_desc)]))
    else:
        header.add_format_line(vcfpy.OrderedDict([('ID', '{}{}'.format(prefix, metric)),
                                                 ('Number', field_number), ('Type', field_type),
                                                 ('Description', field_desc)]))
    return header


def get_bkgd_dict(bkgd_metrics_file):
    bkgd_metrics_dict = {}
    with open(bkgd_metrics_file) as infile:
        columns = next(infile)
        metric_types = [*columns.strip('\n').split('\t')[4:]]
        for line in infile:
            roi_id = '{}:{}{}>{}'.format(*line.strip('\n').split('\t')[:4])
            metrics = line.split('\t')[4:]
            bkgd_metrics_dict[roi_id] = dict(zip(metric_types, metrics))
    return bkgd_metrics_dict


def main():
    args = supply_args()

    reader = vcfpy.Reader.from_path(args.input_vcf)
    header = reader.header

    if args.per_base_bkgd_metrics:
        per_base_mets = get_bkgd_dict(args.per_base_bkgd_metrics)
    else:
        per_base_mets = None

    if args.fc_bkgd_metrics:
        fc_mets = get_bkgd_dict(args.fc_bkgd_metrics)
    else:
        fc_mets = None

    # add header info
    prefix = 'BE'
    metrics = ['mean', 'std', 'median', 'threshold', 'zscore', 'pval']

    for metric in metrics:
        add_metrics_header('INFO', prefix, metric, 1, 'Float',
                           'Descriptive stats of locus allele frequencies from PON samples', header)
    add_metrics_header('INFO', '', args.filter_id, 1, 'Float',
                       'Below background', header)
    if fc_mets:
        add_metrics_header('INFO', 'fc', args.filter_id, 1, 'Float',
                           'Below background for forced calls', header)
    writer = vcfpy.Writer.from_path(args.output_vcf, header=header)

    # add info to VCF records
    for record in reader:
        roi_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT[0].value)
        if args.fc_label in record.FILTER:
            if fc_mets:
                for metric in metrics:
                    record.INFO['{}{}'.format(prefix, metric)] = fc_mets[roi_id][metric]
                if record.calls[0].data['AF'] < float(fc_mets[roi_id]['threshold']):
                    record.add_filter('{}{}'.format(args.fc_label, args.filter_id))
        else:
            if per_base_mets:
                if record.ALT[0].type == 'SNV' and record.calls[0].data['AF'] < float(args.vaf_threshold):
                    if not args.alt_specific:
                        roi_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, 'nontarget', 'nontarget')
                    else:
                        roi_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.ALT[0].value, record.ALT[0].value)
                    if roi_id in per_base_mets:
                        for metric in metrics:
                            record.INFO['{}{}'.format(prefix, metric)] = per_base_mets[roi_id]['threshold']
                        if record.calls[0].data['AF'] < float(per_base_mets[roi_id]['threshold']):
                            record.add_filter(args.filter_id)
                    else:
                        print('{}:{} not found in PON background estimates file'.format(record.CHROM, record.POS))
        writer.write_record(record)


if __name__ == "__main__":
    main()
