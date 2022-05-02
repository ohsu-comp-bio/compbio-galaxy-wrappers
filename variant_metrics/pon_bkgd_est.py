#!/usr/bin/env python

"""
Given a directory containing GATK Depth of Coverage files (with base counts) from a cohort,
calculate the bkgd (base depth/total depth) at each locus for each base
and for each locus and base, from each sample, find the lower and upper fence of the interquartile range.

Example:
python variant_metrics/pon_bkgd_est.py "/Users/onwuzu/Downloads/PON" "/Users/onwuzu/Downloads/bkgd_est_output.txt"
"""

import argparse
import os
import numpy as np
from doc_tools import DepthOfCoverageReader

VERSION = '0.0.1'

BASES = ['A', 'C', 'G', 'T']


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('doc_dir', help='Path to directory containing Depth of Coverage files')
    parser.add_argument('bkgd_est_output',
                        help='Output file of the estimated background stats.')
    parser.add_argument('--version', action='version',
                        version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def calc_bkgd_est(base_depth, total_depth):
    bkgd = int(base_depth) / int(total_depth)
    return float('{:0.3f}'.format(bkgd))


def calc_threshold(bkgd_list):
    q1 = np.quantile(bkgd_list, 0.25)
    q3 = np.quantile(bkgd_list, 0.75)
    iqr = q3 - q1
    lower = q1 - (1.5 * iqr)
    upper = q3 + (1.5 * iqr)
    return lower, upper


def main():

    args = supply_args()

    for dirpath, dirnames, filenames in os.walk(args.doc_dir):
        coll_bkgd = {}
        bkgd_threshold = {}
        for file in filenames:
            pon_file = os.path.join(dirpath, file)
            print(pon_file)
            if file.endswith('.tsv'):
                doc = DepthOfCoverageReader(pon_file).doc
                bkgd_dict = {}
                for locus in doc.keys():
                    for base in BASES:
                        bkgd = calc_bkgd_est(doc[locus][base], doc[locus]['depth'])
                        print(bkgd)
                        bkgd_dict.setdefault(locus, {})
                        bkgd_dict[locus].setdefault(base, bkgd)

                        coll_bkgd.setdefault(locus, {})
                        coll_bkgd[locus].setdefault(base, []).append(bkgd)

        for locus in coll_bkgd:
            for base in BASES:
                bkgd_threshold.setdefault(locus, {})
                bkgd_threshold[locus].setdefault(base, {})
                bkgd_threshold[locus][base]['lower'],  bkgd_threshold[locus][base]['upper'] = calc_threshold(coll_bkgd[locus][base])

    outfile = open(args.bkgd_est_output, 'w')
    outfile.write('Locus\tA\tC\tG\tT\n')
    for locus in bkgd_threshold:
        outfile.write(
            '{}\t{},{}\t{},{}\t{},{}\t{},{}\n'.format(
                locus,
                bkgd_threshold[locus]['A']['lower'],
                bkgd_threshold[locus]['A']['upper'],
                bkgd_threshold[locus]['C']['lower'],
                bkgd_threshold[locus]['C']['upper'],
                bkgd_threshold[locus]['G']['lower'],
                bkgd_threshold[locus]['G']['upper'],
                bkgd_threshold[locus]['T']['lower'],
                bkgd_threshold[locus]['G']['upper']))
    outfile.close()


if __name__ == "__main__":
    main()