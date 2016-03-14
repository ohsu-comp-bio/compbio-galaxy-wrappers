#!/usr/bin/env python

### VCF files list multiple alternate alleles in the ALT columns, separated by commas.
### This script will split entries with multiple allese listed, and place each on a separate
### line.  This will allow us to import this data in to annotation software and properly
### receive prediction scores.
### John Letaw 06/09/15

import argparse

def writeVcfLine(handle, line, allele, alt):

    counts = line[9].split(':')[1].split(',')[0] + ',' + line[9].split(':')[1].split(',')[alt]
    geno = ':'.join(line[9].split(':')[2:])

    for i in range(len(line)):
        if i == len(line)-1:
#            handle.write(line[i] + '\n')
            handle.write("0/1:" + counts + ":" + geno + '\n')
        elif i == 4:
            handle.write(allele + '\t')
        elif i != 4:
             handle.write(line[i] + '\t')

def main():
    
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(dest='input', help='')
    parser.add_argument(dest='output', help='')
    parser.add_argument(dest='output_alt', help='VCF Output composed of alternate alleles.')
    args = parser.parse_args()

    handle_in_vcf = open(args.input, 'rU')
    handle_out_vcf = open(args.output, 'w')
    handle_out_vcf_alt = open(args.output_alt, 'w')

    with handle_in_vcf as vcf:
        for line in vcf:
            if line[0] != "#":                
                new_line = line.rstrip('\n').split('\t')
                if ',' in new_line[4]:
                    alt_allele = new_line[4].split(',')
                    for allele in alt_allele:
                        if alt_allele.index(allele) == 0:
                            writeVcfLine(handle_out_vcf, new_line, allele, 1)
                        else:
                            writeVcfLine(handle_out_vcf_alt, new_line, allele, 2)
                else:
                    handle_out_vcf.write(line)
            else:
                handle_out_vcf.write(line)
                handle_out_vcf_alt.write(line)

    handle_out_vcf.close()
    handle_out_vcf_alt.close()

if __name__ == "__main__":
    main()
