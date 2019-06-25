#!/usr/bin/env python

# Create VCF according to vcf_avg_af_hc.sh.  Then, loop through this on all output VCFs, writing to file.  Like:
# for line in $(ls output/*.vcf); do python avg_af.py $line; done > temp
# Also will create an output json blob for import to sample_metrics.

import argparse
import json
import numpy
import re

VERSION = '0.1.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--infile', help='Input VCF containing profile SNP loci.')
    parser.add_argument('--outfile', help='Output file with json string.')
    parser.add_argument('--runid', help='Supply runid, as opposed to parsing from header.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

def find_run(line):
    """
    Get run id from header.
    ##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller  --genotyping-mode GENOTYPE_GIVEN_ALLELES --alleles /home/gro
    ups/clinical/users/letaw/heme_snp_profiles/qiaseq_common.vcf.gz --output /home/groups/clinical/users/letaw/heme_snp_profiles/out
    put/18KD-022M0040-1.vcf --intervals /home/groups/clinical/users/letaw/heme_snp_profiles/qiaseq_common.interval_list --interval-p
    adding 50 --input /home/exacloud/clinical/PressLab/BAM/180130_NS500390_0196_AH7GVCBGX5/18KD-022M0040-1.bam
    """
    my_regex = "([0-9]{6}_NS500[0-9]{3}_[0-9]{4}_[A-Z0-9]{10})/"
    my_match = re.findall(my_regex, line)
    return my_match[0]


def main():

    args = supply_args()
    geno_cnts = {'0/0': [0.0, 0.0], '0/1': [0.0, 0.0], '1/1': [0.0, 0.0]}

    all_ab = {'0/0': [], '0/1': [], '1/1': []}

    with open(args.infile, 'rU') as myvcf:
        for line in myvcf:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                geno = line[9].split(':')[0]
                ref = float(line[9].split(':')[1].split(',')[0])
                alt = float(line[9].split(':')[1].split(',')[1])
                af = alt / (ref + alt)

                if geno in geno_cnts:
                    geno_cnts[geno][0] += ref
                    geno_cnts[geno][1] += alt

                if geno in all_ab:
                    all_ab[geno].append(af)

            else:
                if line.startswith('##GATKCommandLine') and not args.runid:
                    run_id = find_run(line)
                elif args.runid:
                    run_id = args.runid
                elif line.startswith('#CHROM'):
                    sample_id = line.rstrip('\n').split('\t')[9]

    hom_ref = geno_cnts['0/0'][1]/(geno_cnts['0/0'][0]+geno_cnts['0/0'][1])
    het = geno_cnts['0/1'][1]/(geno_cnts['0/1'][0]+geno_cnts['0/1'][1])
    hom_alt = geno_cnts['1/1'][0]/(geno_cnts['1/1'][0]+geno_cnts['1/1'][1])

    het_std = numpy.std(numpy.array(all_ab['0/1']))

    if hom_ref >= 0.02 or hom_alt >= 0.02:
        to_write = [sample_id, run_id, str(hom_ref), str(hom_alt), str(het)]
        print('\t'.join(to_write))

    print("0/0: " + str(hom_ref))
    print("0/1: " + str(het))
    print("1/1: " + str(hom_alt))

    print("0/0 std: " + str(numpy.std(numpy.array(all_ab['0/0']))))
    print("0/1 std: " + str(het_std))
    print("1/1 std: " + str(numpy.std(numpy.array(all_ab['1/1']))))

    handle_out = open(args.outfile, 'w')
    out_metric = {'allele_balance_hom_ref': hom_ref, 'allele_balance_het_std': het_std}
    json.dump(out_metric, handle_out)
    handle_out.close()

main()
