#!/usr/bin/env python

### Remove entries where SAMPLE DP is 0.
### USAGE: python remove_zero_dp.py <in_vcf> <out_vcf>
# VERSION: 0.1.1

import vcf
import sys

def main():
    
    vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))
    vcf_writer = vcf.Writer(open(sys.argv[2], 'w'), vcf_reader)

    for record in vcf_reader:
        try:
            if record.samples[0].data.DP:
                dp_normal = record.samples[0].data.DP
                dp_tumor = None
                if len(record.samples) > 1:
                    dp_tumor = record.samples[1].data.DP
                if dp_normal == 0 or dp_tumor == 0:
                    print("SAMPLE DP = 0.  Do not write this record to output VCF.")
                    print(record)
                    print(record.samples)
                else:
                    vcf_writer.write_record(record)

        except AttributeError:
            print("No DP value within the SAMPLE field of this VCF.  Record follows.")
            print(record)

if __name__ == "__main__":
    main()


