#!/usr/bin/env python

# USAGE:
# CODED BY: Janice Patterson

from __future__ import print_function
import argparse
import os
import sys

VERSION = '0.1.0'


def supply_args():
    """
    what pirates say
    """
    parser = argparse.ArgumentParser(description='Converts Ensembl transcript identifiers to Refseq using'
                                                 ' Based on BioMart table grch37.ensembl.org (biomaRt_2.34.1)')
    parser.add_argument('input', help='tab delmited matrix with ensembl tx in 1st column')
    parser.add_argument('txIDtrans', help='ensembl2refseq.txt conversion file')
    parser.add_argument('output', help='tab delimited with Refseq ids and HGVS ids')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                                                               VERSION)
    args = parser.parse_args()
    return args


def main():
    args = supply_args()
    coords, refseq = read_txIDtrans(args.txIDtrans)
    converted = convert_enst2refseq(args.input, refseq)
    write_converted(args.output, converted)

def read_txIDtrans(filename):
    """
    :param filename: txIDtrans name of file
    :return: dictionary of ENST:refseq_mrna, ensembl_gene_id, hgnc_symbol, chromosome_name, start_position,end_position
    """
    coords={}
    refseq={}
    with open(filename, 'rU') as trans:
        for line in trans:
            if not line.startswith('ensembl'):
                sline = line.rstrip('\n').split('\t')
                coords[sline[0]] = [sline[4:], sline[3]]
                refseq[sline[0]] = [sline[1], sline[3]]
    return coords, refseq


def convert_enst2refseq(in_filename, refseq):
    """
    :param enst_filename:
    :return:
    """
    converted={}
    with open(in_filename, 'rU') as infile:
        for line in infile:
            if '_' not in line:
                sline = line.rstrip('\n').split('\t')
                refseqID = refseq[sline[0]][0]
                if not refseqID == "":
                    enst_line = sline[1:]
                    sline.insert(1,refseq[sline[0]][1])
                    converted[refseqID] = sline[1:]
    return converted

def write_converted(out_filename, converted):
    with open(out_filename, 'w') as outfile:
        for key, value in sorted(converted.iteritems()):
            to_write = '\t'.join([key, '\t'.join(value)])
            outfile.write(to_write)
            outfile.write('\n')



if __name__ == "__main__":
    main()


