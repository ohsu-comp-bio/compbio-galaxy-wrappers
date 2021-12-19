"""
Select variants in a VCF.

Example usage: var_select.py 'input.vcf' 'fc' 'output.vcf'
var_select.py '/Users/onwuzu/Downloads/test_output_var_label.vcf' 'fc' --exclusive '/Users/onwuzu/Downloads/test_output_var_select.vcf'

Details: Select variants in input VCF.
"""

import argparse
import vcf_tools

VERSION = vcf_tools.VERSION+'.0'


def main():
    parser = argparse.ArgumentParser(description="var_select.py 'file.vcf' 'fc' --exclusive 'output.vcf'")
    parser.add_argument('input_vcf',  help="Input VCF")
    parser.add_argument('filter', nargs="+", help="Filter(s) of interest")
    parser.add_argument('--exclusive', action='store_true', help="Select variants with ONLY filter(s) of interest")
    parser.add_argument('output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()

    var = vcf_tools.VarLabel(args.input_vcf, args.filter)
    records = var.get_filter_records(args.exclusive)
    vcf_tools.VarWriter(records).as_vcf(args.output_vcf, var.reader.header)


if __name__ == "__main__":
    main()
