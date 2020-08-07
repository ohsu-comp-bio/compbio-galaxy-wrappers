#!/usr/bin/python

### Sort VCFs based on contig list below.
### ./sort_vcf.py <in_vcf> <out_sorted_vcf>

import sys

VERSION="0.1.0"

CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 NC_007605"

handle = open(sys.argv[1], 'rU')
handle_out = open(sys.argv[2], 'w')

vcf_dict = {}

with handle as vcf:
    for line in vcf:
        if line[0] == "#":
            handle_out.write(line)
        else:
            line = line.rstrip('\n').split('\t')
            chrom = line[0]
            coord = line[1]
            if chrom not in vcf_dict:
                vcf_dict[chrom] = [line]
            else:
                vcf_dict[chrom].append(line)

for chrom in CHROMS.split(' '):
    if chrom in vcf_dict:
        for entry in sorted(vcf_dict[chrom], key=lambda x:int(x[1])):
            handle_out.write('\t'.join(entry))
            handle_out.write('\n')

handle_out.close()
