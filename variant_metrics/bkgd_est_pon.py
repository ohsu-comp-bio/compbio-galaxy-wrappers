#!/usr/bin/env python

"""
Given a directory containing GATK Depth of Coverage files (with base counts) from a cohort,
calculate the bkgd (base depth/total depth) at each locus for each base
and for each locus and base, from each sample, find the lower and upper fence of the interquartile range.

Example:
python variant_metrics/bkgd_est_pon.py "/Users/onwuzu/Downloads/PON" "/Users/onwuzu/Downloads/bkgd_est_output.txt"
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
    return float(bkgd)


def calc_threshold(bkgd_list):
    q1 = np.quantile(bkgd_list, 0.25)
    q3 = np.quantile(bkgd_list, 0.75)
    iqr = q3 - q1
    lower = q1 - (1.5 * iqr)
    upper = q3 + (1.5 * iqr)
    return float('{:0.3f}'.format(upper))


def main():

    args = supply_args()

    coll_bkgd = {}
    bkgd_threshold = {}
    for dirpath, dirnames, filenames in os.walk(args.doc_dir):
        for file in filenames:
            pon_file = os.path.join(dirpath, file)
            if file.endswith('.tsv'):
                doc = DepthOfCoverageReader(pon_file).doc

                bkgd_dict = {}
                for locus in doc.keys():
                    base_depth_list = []
                    for base in BASES:
                        base_depth_list.append(doc[locus][base])

                        # calculate base bkgd estimates
                        bkgd = calc_bkgd_est(doc[locus][base], doc[locus]['total_depth'])
                        bkgd_dict.setdefault(locus, {})
                        bkgd_dict[locus].setdefault(base, bkgd)

                        # collect all bkgd at a locus
                        coll_bkgd.setdefault(locus, {})
                        coll_bkgd[locus].setdefault(base, []).append(bkgd)

                    # calculate general bkgd estimates
                    bkgd = calc_bkgd_est(doc[locus]['off_target'], doc[locus]['total_depth'])
                    coll_bkgd[locus].setdefault('general', []).append(bkgd)

        # calculate the upper limit threshold for outliers in collected bkgd
        for locus in coll_bkgd:
            for base in BASES:
                bkgd_threshold.setdefault(locus, {})
                bkgd_threshold[locus].setdefault(base, {})
                bkgd_threshold[locus][base]['upper'] = calc_threshold(coll_bkgd[locus][base])
            bkgd_threshold[locus].setdefault('general', {})
            bkgd_threshold[locus]['general']['upper'] = calc_threshold(coll_bkgd[locus]['general'])

    # write calculated upper limit threshold to file
    outfile = open(args.bkgd_est_output, 'w')
    outfile.write('Locus\tA\tC\tG\tT\tGeneral_BE\n')
    for locus in bkgd_threshold:
        outfile.write(
            '{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                locus,
                bkgd_threshold[locus]['A']['upper'],
                bkgd_threshold[locus]['C']['upper'],
                bkgd_threshold[locus]['G']['upper'],
                bkgd_threshold[locus]['T']['upper'],
                bkgd_threshold[locus]['general']['upper']))
    outfile.close()


if __name__ == "__main__":
    main()
