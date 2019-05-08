#!/usr/bin/env python

# DESCRIPTION: Create sample level metrics to be passed to the CGD for RNA.  Metrics
#  are passed as a json dump.
# usage: sample_metrics.py -h
# by: Janice
# todo:
#   optional arguments when tpm is empty.

from __future__ import print_function

import argparse
import json
import numpy
from pybedtools import BedTool

VERSION = '0.22.0'


def main():
    args = supply_args()
    target_genes = targetbed_GTF(args.filtered_gtf)
    sumrz_cnts, hk_counts = summarize_starcounts(args.counts, target_genes, args.hk_file)

    if args.tpm:
        sumrz_tpm, hk_tpm = summarize_tpm(args.tpm, target_genes, args.hk_file)

    write_output(sumrz_cnts, sumrz_tpm, hk_counts, hk_tpm, args.out, args.outjson)


def supply_args():
    """

    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--counts', help='STAR_output_reads_per_gene')
    parser.add_argument('--tpm', help="tpmfile converted from counts, ENSG ids only")
    parser.add_argument('--filtered_gtf', help = 'ref_GRCh37.p13.truseq-rna-exome-targets.gtf')
    parser.add_argument('--hk_file', help='List of housekeeping genes to output')
    parser.add_argument('--out', help='Output file in human readable text format.')
    parser.add_argument('--outjson', help='json output for CGD')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args



def targetbed_GTF(filtered_gtf):
    """
    :return list of HGVS for genes only in the target file
    """

    gtf = BedTool(filtered_gtf)
    targetgenes = set()
    #gtf.filter(lambda x: x[2] == 'gene')
    for feature in gtf:
        attrib = str(feature[8].strip().split(';')[2])
        targetgenes.add(attrib.split('"')[1])


    return(targetgenes)


def get_hkgenes(filename_hk, target_genes, dict_of_abundance):
    """
    Housekeeping genes, filter for those genes are in the target file
    :param list of housekeeping genes
    :param target_genes (list) HGVS name
    :param dict_of_abundance (dict)  of gene names: counts/tpm
    :return: dict_of_hkgenes (dict)
    """
    hk_genes = {}

    with open(filename_hk, 'rU') as fh_hk:
        for lines in fh_hk:
            geneid = lines.rstrip('\n').split('\t')
            if geneid[0] in target_genes:
                hk_genes[geneid[0]] = [geneid[0], dict_of_abundance[geneid[0]]]
            else:
                print(geneid[0] + ": ", "Check gene id format (Hugo/ENSG) or gene is not in target intervals.")
    return(hk_genes)

def summarize_starcounts(infile, target_genes, filename_hk):
    """
    Summarize the reads.per.gene from STAR to TPM for table
    :param infile reads.per.gene
    :param list of target_genes
    :return sumrz_counts (dict)
    """
    rev_counts = {}
    other = {}
    sumrz_counts = {}
    with open(infile, 'rU') as star_counts:
        for lines in star_counts:
            line = lines.rstrip('\n').split('\t')
            if '_' not in line[0]:
                rev_counts[line[0]] = line[3]
            else:
                other[line[0]] = line[3]
    target_counts = [rev_counts[item] for item in rev_counts.keys() if item in target_genes]

    np_counts = numpy.asarray(target_counts, dtype=numpy.float64)
    sumrz_counts["COUNT_0"] = len(np_counts[np_counts >= 0])
    for i in range(0,6):
        idx = pow(10, i)
        sumrz_counts["COUNT_" + str(idx)] = len(np_counts[np_counts >= idx])

    hk_counts = get_hkgenes(filename_hk, target_genes, rev_counts)
    return(sumrz_counts, hk_counts)


def summarize_tpm(filename_tpm, target_genes, filename_hk):
    """
    Summarize the converted TPM (calc_norm.py) reads.per.gene from STAR to TPM for table
    :param filename_tpm (string) filename of tpm translated output
    :param sumrz_counts (dict)
    :param target_genes (list)  target HUGO gene ids
    :param filename_hk (string) filename of housekeeping gene list
    :return sumrz_tpm (dict)
    :return hk_genes (list of list)
    """
    tpm = {}
    sumrz_tpm = {}  #tuple, ordered
    with open(filename_tpm, 'rU') as fh_tpm:
        for lines in fh_tpm:
            line = lines.rstrip('\n').split('\t')
            if line[0] in target_genes:
                tpm[line[0]] = line[1]


    #Summarize TPM
    np_tpm = numpy.asarray(tpm.values(), dtype=numpy.float64)
    sumrz_tpm["TPM_0"] = len(np_tpm[np_tpm >= 0])
    for i in range(-2,4):
        idx = pow(10, i)
        sumrz_tpm["TPM_" + str(idx)] = len(np_tpm[np_tpm >= idx])



    hk_tpm = get_hkgenes(filename_hk, target_genes, tpm)

    return (sumrz_tpm, hk_tpm)

def write_output(sumrz_cnts, sumrz_tpm, hk_counts, hk_tpm, outfile, outjson):
    """
    Write metrics to a text file, mainly to be viewed in Galaxy.
    json output for CGD, works for now, <<shrug>>.
    :return:
    """

    with open(outfile, 'w') as out:
        out.write('COUNTS_greaterorequalto_Limit\tnumber_of_genes\n')
        for k in sorted(sumrz_cnts):
            out.write('%s\t%s\n' % (k, sumrz_cnts[k]))

        out.write('\n')
        out.write('TPM_greaterorequalto_Limit\tnumber_of_genes\n')

        for k in sorted(sumrz_tpm):
            out.write('%s\t%s\n' % (k, sumrz_tpm[k]))

        out.write('\n')
        out.write('HOUSEKEEPING GENES\n')
        out.write('hgvs_gene_id\tensembl_id\tcounts\ttpm\n')

        #sort by HUGO gene name, because life is hard
        for key, value in sorted(hk_counts.iteritems(), key=lambda (k, v): (v, k)):
            out.write('%s\t%s\t%s\t%s\n' % (value[0], key, str(value[1]), hk_tpm[key][1]))

        out.write('\n')


    with open(outjson, 'w') as jsonout:
        json_dict = {"COUNTS": sumrz_cnts, "TPM": sumrz_tpm}
        jsonout.write(json.dumps(json_dict, sort_keys=True, indent=4, separators=(',', ': ')))



if __name__ == "__main__":
    main()
