"""
Select records/variants in VCF.

Example usage: record_var_select.py --input_vcf 'input.vcf' --filter fc --output_vcf 'output.vcf'

Details: Select variants in input VCF.
"""

import argparse
import vcf

VERSION = vcf.VERSION+'.0'


def main():
    parser = argparse.ArgumentParser(
        description="var_select.py --input_vcf 'file.vcf' "
                    "--filter 'fc' "
                    "--exclusive True "
                    "--output_vcf 'output.vcf'")
    parser.add_argument('--input_vcf',
                        help="VCF file from which to select info.")
    parser.add_argument('--filter', nargs="+",
                        help="select FILTER(s) of interest.")
    parser.add_argument('--exclusive',
                        help="select variants with only FILTER(s) of interest.")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()

    records = vcf.Records(args.input_vcf).get_filter_records(args.filter, args.exclusive)

    vcf.RecordsWriter(records).as_vcf(args.output_vcf, records.reader.header)


if __name__ == "__main__":
    main()
