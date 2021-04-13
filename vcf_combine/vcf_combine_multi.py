#!/usr/bin/env python

# DESCRIPTION: vcf_combine_multi.py --input_vcfs "file1.vcf" "file2.vcf" "fileN.vcf" --output_vcf "output.vcf"
# BY: O.K

"""
Combine VCFs produced by various tools and keep only unique variants.
"""

import argparse
import pandas as pd
from functools import reduce

VERSION = '0.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description="vcf_combine_multi.py --input_vcfs 'file1.vcf' 'file2.vcf' 'fileN.vcf' --output_vcf 'output.vcf'")
    parser.add_argument('-i', '--input_vcfs', nargs="+", help="vcf files from several tools")
    parser.add_argument('-o', '--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()
    return args


def get_header(vcffile):
    i = 0
    with open(vcffile, 'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                i += 1
                if line.startswith('##source'):
                    source = line.split()[0].split('=')[1]
    return i, source

def main():

    args = supply_args()

    vcfs = args.input_vcfs

    vars = []
    for vcf in vcfs:
        header_info = get_header(vcf)
        var = pd.read_csv(vcf, sep='\t', header=header_info[0])
        #TODO: keep duplicates in a single vcf?
        var = var.drop_duplicates()
        var = var.astype({'QUAL': object})
        var[['QUAL', 'ID']] = '.'
        var.loc[var['FILTER'] != '.', 'FILTER'] = header_info[1] + '=' + var['FILTER']
        var.loc[var['FILTER'] == '.', 'FILTER'] = header_info[1]
        vars.append(var)

    match_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL']
    merged = reduce(lambda left, right: pd.merge(left, right, on=match_cols, how='outer'), vars)

    for col in var.columns[6:]:
        join_cols = [c for c in merged if c.startswith(col)]
        merged[col] = merged[join_cols].apply(lambda x: '|'.join(x.dropna()), axis=1)

    merged = merged[var.columns]
    merged = merged.sort_values(['#CHROM', 'POS'], ascending=(True, True))

    with open(args.output_vcf, 'w') as outfile, open(args.input_vcfs[0], 'r') as infile:
        for line in infile:
            if line.startswith('##'):
                outfile.write(line)

    merged.to_csv(args.output_vcf, sep='\t', mode='a', index=False)


if __name__ == '__main__':
    main()
