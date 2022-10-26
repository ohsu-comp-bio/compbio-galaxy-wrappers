#!/usr/bin/env python

"""
Given a sample VCF and sample Depth of Coverage file, collect the depth at each locus and base
and annotate on the sample VCF INFO field, the estimated bkgd threshold from a cohort of samples (given as a text file)
and on the FORMAT field, the calculated bkgd (base depth/total depth)

Example:
python variant_metrics_umi/bkgd_metrics.py "input_vcf.vcf" "input_doc.tsv" "bkgd_est_file.txt" "output_bkgd_metrics.vcf"
"""

import argparse
import vcfpy

from vcf_tools import VarWriter
from doc_tools import DepthOfCoverageReader, calc_bkgd_est


VERSION = '0.0.1'

BASES = ['A', 'C', 'G', 'T']


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', help='Input VCF')
    parser.add_argument('input_cov', type=DepthOfCoverageReader, help='Input Depth of Coverage file')
    parser.add_argument('bkgd_est_file', help='Text file containing depth stats of loci in panel.')
    parser.add_argument('output_vcf', help='Output VCF')
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def get_est_dict(bkgd_est_file):
    bkgd_threshold_dict = {}
    with open(bkgd_est_file) as f:
        header_line = next(f)
        for line in f:
            lsplit = line.strip().split('\t')
            bkgd_threshold_dict.setdefault(lsplit[0], {})
            for i, base in enumerate(BASES):
                bkgd_threshold_dict[lsplit[0]].setdefault(base, {})
                bkgd_threshold_dict[lsplit[0]][base] = {'upper': float(lsplit[i+1])}
            bkgd_threshold_dict[lsplit[0]]['general'] = {'upper': float(lsplit[-1])}
    return bkgd_threshold_dict


def main():
    args = supply_args()

    reader = vcfpy.Reader.from_path(args.input_vcf)
    header = reader.header
    sample = reader.header.samples.names[0]

    doc_reader = args.input_cov.doc

    res_est = get_est_dict(args.bkgd_est_file)

    info_id = 'BET'
    fmt_id = 'BE'

    # add header info
    header.add_info_line(vcfpy.OrderedDict([('ID', info_id),
                                            ('Number', 1), ('Type', 'Float'),
                                            ('Description', 'Background estimate threshold from cohort of normals (general)')]))
    header.add_format_line(vcfpy.OrderedDict([('ID', fmt_id),
                                              ('Number', 1), ('Type', 'Float'),
                                              ('Description', 'Background estimate')]))
    for base in BASES:
        header.add_info_line(vcfpy.OrderedDict([('ID', '{}_{}'.format(info_id, base)),
                                                ('Number', 1), ('Type', 'Float'),
                                                ('Description', 'Background estimate threshold from cohort of normals for base {}'.format(base))]))
        header.add_format_line(vcfpy.OrderedDict([('ID', '{}_{}'.format(fmt_id, base)),
                                                  ('Number', 1), ('Type', 'Float'),
                                                  ('Description', 'Background estimate for base {}'.format(base))]))

    records = {}
    # add info to VCF records
    for record in reader:
        pos = str(record.CHROM) + ':' + str(record.POS)
        if pos in res_est:
            record.INFO[info_id] = res_est[pos]['general']['upper']
            record.add_format(key=fmt_id, value=None)
            bkgd_est = calc_bkgd_est(doc_reader[pos]['off_target'], doc_reader[pos]['total_depth'])
            record.call_for_sample[sample].data[fmt_id] = float('{:0.3f}'.format(bkgd_est))
            for base in BASES:

                # add cohort bkgd estimates info VCF INFO field
                base_hid = '{}_{}'.format(info_id, base)
                record.INFO[base_hid] = res_est[pos][base]['upper']

                # add calculated AF/background estimates of variant to VCF SAMPLE field
                base_hid = '{}_{}'.format(fmt_id, base)
                record.add_format(key=base_hid, value=None)
                bkgd_est = calc_bkgd_est(doc_reader[pos][base], doc_reader[pos]['total_depth'])
                record.call_for_sample[sample].data[base_hid] = float('{:0.3f}'.format(bkgd_est))
        else:
            print('{} not found in PON background estimates file'.format(pos))

        var_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
        records[var_id] = record

    VarWriter(list(records.values())).as_vcf(args.output_vcf, header)


if __name__ == "__main__":
    main()
