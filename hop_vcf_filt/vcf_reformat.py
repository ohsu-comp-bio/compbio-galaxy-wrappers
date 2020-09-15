#!/usr/bin/env python

import argparse

VERSION = '0.0.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('infile', help='Input VCF')
    parser.add_argument('outfile', help='Output VCF')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

def find_field(info, field):
    """
    Find a field from the INFO column.
    :return:
    """
    for entry in info.split(';'):
        key = entry.split('=')[0]
        val = entry.split('=')[1]
        if key == field:
            return val
    return ""

def main():
    args = supply_args()
    handle_out = open(args.outfile, 'w')
    header = ['SAMPLE ID', 'GENOTYPE', 'REF COUNT', 'ALT COUNT', 'LOCAL AF', 'CHROM', 'COORD', 'REF ALLELE', 'ALT ALLELE',
              'HGMD', 'GNOMAD AF', 'CLINVAR SIG', 'CLINVAR CONFLICTING', 'HGVS G', 'HGVS C', 'HGVS P', 'SNPEFF']
    handle_out.write('\t'.join(header))
    handle_out.write('\n')
    with open(args.infile, 'rU') as myfile:
        for line in myfile:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                chrom = line[0]
                coord = line[1]
                ref = line[3]
                alt = line[4]
                info = line[7]
                samples = line[9:]
                af = find_field(info, 'AF')
                hgmd = find_field(info, "hgmd.CLASS")
                clnsig = find_field(info, "clinvar.CLNSIG")
                clnsigconf = find_field(info, "clinvar.CLNSIGCONF")
                snpeff = find_field(info, "ANN")
                hgvs_g = find_field(info, "HGVS_G")
                hgvs_c = find_field(info, "HGVS_C")
                hgvs_p = find_field(info, "HGVS_P")
                gnomad = find_field(info, "gnomad.AF")
                for sample in samples:
                    geno = sample.split(':')[0]
                    ref_cnt = sample.split(':')[1].split(',')[0]
                    alt_cnt = sample.split(':')[1].split(',')[1]
                    if geno != "0/0" and geno != "./.":
                        to_write = [sample_names[samples.index(sample)]]
                        to_write.extend([geno, ref_cnt, alt_cnt, af])
                        to_write.extend([chrom, coord, ref, alt, hgmd, gnomad, clnsig, clnsigconf, hgvs_g, hgvs_c, hgvs_p, snpeff])
                        handle_out.write('\t'.join(to_write))
                        handle_out.write('\n')

            elif line.startswith('#CHROM'):
                sample_names = line.rstrip('\n').split('\t')[9:]



    handle_out.close()

if __name__ == "__main__":
    main()
