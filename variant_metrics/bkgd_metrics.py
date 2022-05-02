#!/usr/bin/env python

"""
Given a sample VCF and sample Depth of Coverage file, collect the depth at each locus and base
and annotate on the sample VCF INFO field, the estimated bkgd threshold from a cohort of samples (given as a text file)
and on the FORMAT field, the calculated bkgd (base depth/total depth)

Example:
python variant_metrics/bkgd_metrics.py "/Users/onwuzu/Downloads/Galaxy115-[Sort_VCF_on_data_114__Sorted_VCF].vcf"
"/Users/onwuzu/Downloads/input_doc.tsv" --bkgd_est_file "/Users/onwuzu/Downloads/bkgd_est_output.txt" "/Users/onwuzu/Downloads/output_bkgd_metrics.vcf"
"""

import argparse
import vcfpy

from pon_bkgd_est import calc_bkgd_est
from vcf_tools import VarWriter
from doc_tools import DepthOfCoverageReader


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
    parser.add_argument('--bkgd_est_file', help='Text file containing depth stats of loci in panel.')
    parser.add_argument('output_vcf', help='Output VCF')
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def get_est_dict(bkgd_est_file):
    bkgd_threshold_dict = {}
    with open(bkgd_est_file) as f:
        i = 0
        for line in f:
            print(line)
            i += 1
            if i > 1:
                lsplit = line.strip().split('\t')
                bkgd_threshold_dict.setdefault(lsplit[0], {})
                for i, base in enumerate(BASES):
                    print(i)
                    print(lsplit)
                    bkgd_threshold_dict[lsplit[0]].setdefault(base, {})
                    bkgd_threshold_dict[lsplit[0]][base] = {
                        'lower': float(lsplit[i+1].split(',')[0]),
                        'upper': float(lsplit[i+1].split(',')[1])
                    }
    return bkgd_threshold_dict


def main():
    args = supply_args()

    hid = 'BE'
    htype = 'Float'
    hdesc = 'Background estimate lower and upper threshold'

    reader = vcfpy.Reader.from_path(args.input_vcf)

    header = reader.header

    doc_reader = args.input_cov.doc

    res_est = get_est_dict(args.bkgd_est_file)

    # add header info
    for base in BASES:
        base_hid = '{}_{}'.format(hid, base)
        header.add_info_line(vcfpy.OrderedDict([('ID', base_hid),
                                                ('Number', 2), ('Type', htype),
                                                ('Description', '{} for base {}'.format(hdesc, base))]))
        header.add_format_line(vcfpy.OrderedDict([('ID', base_hid),
                                                  ('Number', 1), ('Type', htype),
                                                  ('Description', '{} for base {}'.format(hdesc, base))]))

    records = {}
    for record in reader:

        var_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)

        pos = str(record.CHROM) + ':' + str(record.POS)
        if pos in res_est:
            for base in BASES:
                pos_base_est = res_est[pos][base]
                base_hid = '{}_{}'.format(hid, base)

                record.INFO[base_hid] = [pos_base_est['lower'], pos_base_est['upper']]

                record.add_format(key=base_hid, value=None)
                record.call_for_sample[reader.header.samples.names[0]].data[base_hid] = calc_bkgd_est(doc_reader[pos][base],
                                                                                                      doc_reader[pos]['depth'])
        else:
            raise KeyError('{} not found in PON background estimates file'.format(pos))

        records[var_id] = record

    VarWriter(list(records.values())).as_vcf(args.output_vcf, header)


if __name__ == "__main__":
    main()
