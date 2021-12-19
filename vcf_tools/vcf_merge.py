"""
Merge VCFs produced by various variant callers.

Example usage: vcf_merge.py  --input_vcfs file1 file2 --caller_labels m2 fb --output_vcf output.vcf

Details: Given multiple vcfs produced by various variant callers, this will merge all vcfs into one.
For variants called by multiple callers, it will output a single record and the info column and sample column
will contain those of all callers that made the call.
"""

import argparse
import vcf_tools

VERSION = vcf_tools.VERSION+'.0'


def main():
    parser = argparse.ArgumentParser(
        description="vcf_merge.py --input_vcfs 'file1.vcf' 'file2.vcf' 'fileN.vcf' "
                    "--caller_labels label1 label2 labelN "
                    "--output_vcf 'output.vcf'")
    parser.add_argument('--input_vcfs', nargs="+", help="Input VCFs")
    parser.add_argument('--caller_labels', nargs="+", help="Labels for each input vcf")
    parser.add_argument('--output_vcf', help="Output VCF.")
    parser.add_argument('-v', '--version', action='version', version="%(prog)s " + VERSION)
    args = parser.parse_args()

    merge = vcf_tools.VcfMerger(args.input_vcfs, args.caller_labels)
    header = merge.merge_headers()
    records = merge.merge_records()
    vcf_tools.VarWriter(records).as_vcf(args.output_vcf, header)


if __name__ == "__main__":
    main()
