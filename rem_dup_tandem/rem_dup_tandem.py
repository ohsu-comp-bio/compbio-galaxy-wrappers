#!/usr/bin/env python

### Remove <DUP:TANDEM>, <DUP>, <INV>, and <RPL> tags that mucks things up in Oncotator and SeattleSeq.
### Usage: python rem_dup_tandem.py <input_vcf> <output_vcf>

import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("input_vcf", help="VCF with DUP:TANDEM tag in ALT column.")
    parser.add_argument("output_vcf", help="VCF output.")
    args = parser.parse_args()

    handle_out = open(args.output_vcf, 'w') 

    with open(args.input_vcf, 'rU') as vcf:
        for line in vcf:
            if line[0] != "#":
                alt = line.split('\t')[4]
                if alt[0] != "<":
                    handle_out.write(line)
            else:
                handle_out.write(line)

    handle_out.close()

if __name__ == "__main__":
    main()

