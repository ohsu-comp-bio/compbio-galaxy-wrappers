"""
Filter records in VCF based on filter of interest.

Example usage: var_filter.py 'file.vcf' --filter_from BET --filter_on BE --alt_allele 'output.vcf'

Details: Given a VCF, filter records using given INFO field and on FORMAT field.
Assumes that VCF records with multiple alternate alleles have been split.
"""

import argparse
import vcf_tools

VERSION = vcf_tools.VERSION + '.0'


def main():
    parser = argparse.ArgumentParser(
        description="var_filter.py --input_vcf 'file.vcf' "
                    "--filter_from BET"
                    "--filter_on BE "
                    "--alt_allele"
                    "output_vcf 'output.vcf'")
    parser.add_argument('input_vcf', help="Input VCF")
    parser.add_argument('--filter_from', help="Info fields to filter from")
    parser.add_argument('--filter_on', help="Format fields to filter on")
    parser.add_argument('--alt_allele', action='store_true', help="Select fields of ALT allele e.g BE_A")
    parser.add_argument('output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()

    to_filter = vcf_tools.VarFilter(args.input_vcf)
    header = to_filter.header
    above, below = to_filter.filter_records(args.filter_from, args.filter_on, args.alt_allele)

    vcf_tools.VarWriter(above).as_vcf(args.output_vcf, header)
    vcf_tools.VarWriter(below).as_vcf("filtered_out.vcf", header)


if __name__ == "__main__":
    main()
