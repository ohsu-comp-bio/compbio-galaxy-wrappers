"""
Label variants in a VCF.

Example usage: var_label.py "input.vcf" "resource.vcf" "ID" "description" "output.vcf"

Details: Label variants in input VCF found in a resource VCF.
"""

import argparse
import vcf_tools

VERSION = vcf_tools.VERSION+'.0'


def main():
    parser = argparse.ArgumentParser(description='var_label.py "input.vcf" "resource.vcf" "ID" "description" "output.vcf"')
    parser.add_argument('input_vcf', help="Input VCF")
    parser.add_argument('resource', help="Resource VCF")
    parser.add_argument('label', help="ID of label")
    parser.add_argument('description', help="Description of label")
    parser.add_argument('output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()

    var = vcf_tools.VarLabel(args.input_vcf, args.label)
    records = var.label_records(args.resource)
    var.add_label_info(args.description)
    vcf_tools.VarWriter(records).as_vcf(args.output_vcf, var.reader.header)


if __name__ == "__main__":
    main()
