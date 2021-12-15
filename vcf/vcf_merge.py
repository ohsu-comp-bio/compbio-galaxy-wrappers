"""
Merge VCFs produced by various variant callers.

Example usage: vcf_merge.py  --input_vcfs file1 file2 --caller_labels m2 fb --output_vcf output.vcf

Details: Given multiple vcfs produced by various variant callers, this will merge all vcfs into one.
For variants called by multiple callers, it will output a single record and the info column and sample column
will contain those of all callers that made the call.
"""

import argparse
import vcf

VERSION = vcf.VERSION+'.0'


def main():
    parser = argparse.ArgumentParser(
        description="vcf_merge.py --input_vcfs 'file1.vcf' 'file2.vcf' 'fileN.vcf' "
                    "--caller_labels label1 label2 labelN "
                    "--output_vcf 'output.vcf'")
    parser.add_argument('--input_vcfs', nargs="+",
                        help="vcf files from multiple variant callers")
    parser.add_argument('--caller_labels', nargs="+",
                        help="Labels for each input vcf.")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()

    header = vcf.VcfMerger(args.input_vcfs, args.caller_labels).merge_headers()

    records = vcf.VcfMerger(args.input_vcfs, args.caller_labels).merge_records()

    vcf.RecordsWriter(records).as_vcf(args.output_vcf, header)


if __name__ == "__main__":
    main()
