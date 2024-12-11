"""
Reformat a VCF by replacing occurrences of the GT value '1' with the value '1/1'.
"""

import vcfpy
import argparse

VERSION = '0.0.1'

def supply_args():
    """
    Populate arguments
    """
    parser = argparse.ArgumentParser(description='Reformat GT field in DRAGEN VCFs')
    parser.add_argument('input', help='Input DRAGEN VCF')
    parser.add_argument('output', help='Output Reformatted DRAGEN VCF')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

def main():
    # Collect args, open input, open output
    args = supply_args()
    in_handle = open(args.input)
    out_handle = open(args.output, 'w')

    # read input vcf and prepare to write output vcf
    with in_handle, out_handle:
        reader = vcfpy.Reader(in_handle)
        writer = vcfpy.Writer.from_path(args.output, reader.header)

        # reformat any GT calls that are '1' to '1/1' and write to output vcf
        for vrnt in reader:
            for call in vrnt.calls:
                if call.data['GT'] == '1':
                    call.data['GT'] = '1/1'
            writer.write_record(vrnt)

if __name__ == "__main__":
    main()



