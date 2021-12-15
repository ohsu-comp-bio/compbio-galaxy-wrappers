"""
Select fields of interest from a merged VCF.

Example usage: vcf_select.py --input_vcf 'file.vcf' --info_fields AF DP TLOD --sample_fields AF AD DP PUMI
--caller_priority m2 fb --output_vcf 'output.vcf'

Details: Given a merged VCF (VCF with calls from various variant callers), this will create a new vcf
containing info fields and sample fields passed to --info_fields and --sample_fields, respectively.
A field not in the input VCF will raise an error.
"""

import argparse
import vcf

VERSION = vcf.VERSION+'.0'


def main():
    parser = argparse.ArgumentParser(
        description="vcf_select.py --input_vcf 'file.vcf' "
                    "--info_fields AF DP TLOD "
                    "--sample_fields AF AD DP PUMI"
                    "--caller_priority label1 label2 labelN "
                    "--output_vcf 'output.vcf'")
    parser.add_argument('--input_vcf',
                        help="vcf file from which to select info.")
    parser.add_argument('--info_fields', nargs="+",
                        help="select INFO fields of interest.")
    parser.add_argument('--sample_fields', nargs="+",
                        help="select sample fields of interest.")
    parser.add_argument('--caller_priority', nargs="+",
                        help="set the priority of variant caller.")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()

    header = vcf.VcfSelecter(args.input_vcf, args.info_fields, args.sample_fields, args.caller_priority).select_headers()

    records = vcf.VcfSelecter(args.input_vcf, args.info_fields, args.sample_fields, args.caller_priority).select_records()

    vcf.RecordsWriter(records).as_vcf(args.output_vcf, header)


if __name__ == "__main__":
    main()
