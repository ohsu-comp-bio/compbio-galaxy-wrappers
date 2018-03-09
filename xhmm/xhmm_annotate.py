#!/usr/bin/env python

# DESCRIPTION:
# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
import argparse
from collections import OrderedDict

VERSION = '0.0.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('in_vcf', help='Input VCF')
    parser.add_argument('probeqc', help='Input Probe QC')
    parser.add_argument('sample_id', help='Sample ID')
    parser.add_argument('outfile', help='Output VCF')
    parser.add_argument('outfile_table', help='Output Table')

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def find_sample_index(header, sample_id):
    """
    Get the index of the sample ID from the VCF column header.
    :return:
    """
    try:
        return header.index(sample_id)
    except:
        raise ValueError("Sample ID does not exist in input VCF.")


def split_line(line, idx):
    """
    Split apart the line based on index.
    :param line:
    :return:
    """
    nline = line[:9]
    nline.append(line[idx])
    return nline


def trans_call(value):
    """
    Take the raw 1 or 2 value from the VCF and translate it to DEL or DUP.
    :return:
    """
    if value == '1':
        return 'DEL'
    elif value == '2':
        return 'DUP'
    elif value == '.':
        return 'NO CALL'
    else:
        raise Exception('Call is not del or dup, please investigate.')


def parse_probeqc(filename):
    """
    Parse the ProbeQC.
    1	1634904	1635073	NM_024011.2;NM_033529.2;XM_005244780.1;XM_005244781.1;XM_005244782.1;XM_005244783.1;XM_005244784.1;XM_005244785.1;XM_005244786.1;XM_005244787.1;XM_005244788.1;XM_005244789.1;XM_005244790.1
    CDK11A	30.0	86.7	0.0	0.0	0.0	80.0	100.0
    :return:
    """
    pqc = OrderedDict()
    with open(filename, 'rU') as probeqc:
        for line in probeqc:
            line = line.rstrip('\n').split('\t')
            uniq_key = (line[0], line[1])
            hgnc = line[4]
            pqc[uniq_key] = hgnc
    return pqc


def find_genes(coord, numt, pqc):
    """
    Find all gene names that overlap the regions.
    :param coord:
    :param pqc:
    :return:
    """
    chrom = coord.split(':')[0]
    start = coord.split(':')[1].split('-')[0]
    to_find = (chrom, start)
    genes = []

    idx = pqc.keys().index(to_find)

    for i in range(idx, idx+numt):
        this_gene = pqc.values()[i]
        if this_gene not in genes:
            genes.append(this_gene)

    return genes


def manip_line(nline):
    """
    Take the information from a VCF line and adjust it to be more useful.
    1	1634904	1:1634904-1639033	<DIP>	<DEL>,<DUP>	.	.	AC=7,3;AF=0.12,0.05;AN=59;END=1639033;IMPRECISE;SVLEN=4130;SVTYPE=CNV;TPOS=16
    34904;TEND=1639033;NUMT=8;GQT=17;PREVTARGSTART=1268875;PREVTARGEND=1269851;POSTTARGSTART=1639608;POSTTARGEND=1639694	GT:NDQ:DQ:EQ:SQ:NQ:LQ:RQ:PL:RD:ORD:DS
    CVR	1:99:0:17,0:99,0:0,99:17,0:6,0:104,0,255:-3.38:30.08:Y
    :return:
    """
    region = nline[2]
    this_call = nline[9].split(':')[0]
    cnv_call = trans_call(this_call)
    pop_dels = nline[7].split(';')[0].split('=')[1].split(',')[0]
    pop_dups = nline[7].split(';')[0].split('=')[1].split(',')[1]
    pop_dels_af = nline[7].split(';')[1].split('=')[1].split(',')[0]
    pop_dups_af = nline[7].split(';')[1].split('=')[1].split(',')[1]
    svlen = nline[7].split(';')[5].split('=')[1]
    zscore = nline[9].split(':')[9]
    avg_depth = nline[9].split(':')[10]
    numt = nline[7].split(';')[9].split('=')[1]
    allele_count = nline[7].split(';')[2].split('=')[1]
    new_info = [region, cnv_call, svlen, zscore, avg_depth, numt, allele_count, pop_dels, pop_dels_af, pop_dups, pop_dups_af]
    return new_info


def tbl_header():
    """
    Provide the table header.
    :return:
    """
    header = ['REGION', 'DEL/DUP', 'CNV LENGTH', 'ZSCORE', 'MEAN DEPTH', 'NUMBER OF PROBES', 'TOTAL ALLELES',
              'POP DEL COUNT', 'POP DEL AF', 'POP DUP COUNT', 'POP DUP AF', 'GENES']
    return header


def main():

    args = supply_args()
    probeqc = parse_probeqc(args.probeqc)
    handle_out = open(args.outfile, 'w')
    handle_out_tbl = open(args.outfile_table, 'w')
    handle_out_tbl.write('\t'.join(tbl_header()))
    handle_out_tbl.write('\n')

    with open(args.in_vcf, 'rU') as myvcf:
        for line in myvcf:
            if line.startswith('#CHROM'):
                line = line.rstrip('\n').split('\t')
                samp_idx = find_sample_index(line, args.sample_id)
                nline = split_line(line, samp_idx)
                handle_out.write('\t'.join(nline))
                handle_out.write('\n')
            elif line.startswith('#'):
                handle_out.write(line)
            else:
                line = line.rstrip('\n').split('\t')
                nline = split_line(line, samp_idx)
                if nline[9].split(':')[-1] == 'Y':
                    handle_out.write('\t'.join(nline))
                    handle_out.write('\n')
                    tbl_line = manip_line(nline)
                    tbl_line.append(','.join(find_genes(tbl_line[0], int(tbl_line[5]), probeqc)))
                    handle_out_tbl.write('\t'.join(tbl_line))
                    handle_out_tbl.write('\n')

    handle_out.close()
    handle_out_tbl.close()


if __name__ == "__main__":
    main()
