"""
Label records/variants in VCF(s).

Example usage: var_label.py --input_vcf 'input.vcf' '--resource resource.vcf --label ID --description Description
--output_vcf 'output.vcf'

Details: Label variants in input VCF found in a resource VCF.
"""

import argparse
import vcf

VERSION = vcf.VERSION+'.0'


def main():
    parser = argparse.ArgumentParser(
        description="var_label.py --input_vcf 'input.vcf' "
                    "--resource resource.vcf "
                    "--label ID "
                    "--description Description"
                    "--output_vcf 'output.vcf'")
    parser.add_argument('--input_vcf',
                        help="vcf file(s) to be labelled")
    parser.add_argument('--resource', help="Resource VCF.")
    parser.add_argument('--label', help="ID of filter.")
    parser.add_argument('--description', help="Description of filter.")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()

    records = vcf.Records(args.input_vcf).label_records(args.resource, args.label, args.description)

    vcf.RecordsWriter(records).as_vcf(args.output_vcf, records.reader.header)


if __name__ == "__main__":
    main()
