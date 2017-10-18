#!/usr/bin/env python

# DESCRIPTION: Create sample level metrics to be passed to the CGD.  Metrics
#  are passed as a json dump.
# USAGE: sample_metrics.py -h
# CODED BY: John Letaw

from __future__ import print_function

import argparse
import json
import sys
sys.path.append('/home/groups/clinical/users/letaw/jhl_tools')

from file_types.gatk_intervals import ProbeQcRead
from file_types.picard import AlignSummaryMetrics

VERSION = '0.1.4'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--probeqc_after', type=ProbeQcRead,
                        help='Probe coverage QC after UMI deduplication '
                             'metrics.')
    parser.add_argument('--probeqc_before', type=ProbeQcRead,
                        help='Probe coverage QC before UMI deduplication '
                             'metrics.')
    parser.add_argument('--picard_summary', help='Picard alignment summary '
                                                 'metrics file.')
    parser.add_argument('--outfile', help='Output file with json string.')
    parser.add_argument('--outfile_txt', help='Output file in human readable '
                                              'text format.')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                                                               VERSION)
    args = parser.parse_args()
    return args


def calc_total_bp(probeqc):
    """
    Get the total number of base pairs covered by your targeted region set.
    :return:
    """
    total_bp = 0
    for line in probeqc.probeqc.itervalues():
        try:
            curr = int(line['STOP']) - int(line['START']) + 1.0
            total_bp += curr
        except ValueError:
            pass

    return total_bp


def calc_cov(probeqc, metric):
    """
    Calculate total coverage across sample.
    :param probeqc:
    :param metric:
    :return:
    """
    total_cov = 0
    for line in probeqc.probeqc.itervalues():
        try:
            curr_bp = int(line['STOP']) - int(line['START']) + 1.0
            curr_cov = curr_bp * float(line[metric])
            total_cov += curr_cov
        except ValueError:
            pass

    return total_cov


def calc_metric(bp, cov):
    """
    Calculate the sample level AVGD from total bp anc coverage.
    :param bp:
    :param cov:
    :return:
    """
    return '{:0.1f}'.format(cov / bp)


def add_on_target(picard, total_cov):
    """
    Include the percent on target reads metric, mainly for amplicon assays.
    :return:
    """
    if 'PAIR' in picard:
        pf_bases_aligned = int(picard['PAIR']['PF_ALIGNED_BASES'])
    elif 'UNPAIRED' in picard:
        pf_bases_aligned = int(picard['UNPAIRED']['PF_ALIGNED_BASES'])
    else:
        pf_bases_aligned = None

    on_target = str("{:.4}".format((total_cov * 100.0) / pf_bases_aligned))

    return on_target


def write_to_text(sample_metrics, outfile_txt):
    """
    Write metrics to a text file, mainly to be viewed in Galaxy.
    :return:
    """
    for key, value in sorted(sample_metrics.iteritems()):
        outfile_txt.write(key)
        outfile_txt.write(': ')
        outfile_txt.write(value)
        if key != 'AVGD':
            outfile_txt.write('%')
        outfile_txt.write('\n')

    outfile_txt.close()


def map_fields(sample_metrics):
    """
    CGD recognizes fields as described here.
    Q30: qthirty
    D20: depthTwenty
    D100: depthOneHundred
    D500: depthFiveHundred
    D2000: depthTwoThousand
    AVGD: averageDepth
    pumi: percentUmi
    on_target: percentOnTarget
    :return:
    """
    mapped_metrics = {}
    mapping = {'Q30': 'qthirty',
                      'D10': 'depthTen',
                      'D20': 'depthTwenty',
                      'D50': 'depthFifty',
                      'D100': 'depthOneHundred',
                      'D250': 'depthTwoHundredFifty',
                      'D500': 'depthFiveHundred',
                      'D700': 'depthSevenHundred',
                      'D1250': 'depthTwelveHundredFifty',
                      'D2000': 'depthTwoThousand',
                      'AVGD': 'averageDepth',
                      'pumi': 'percentUmi',
                      'on_target': 'percentOnTarget'}

    for key, value in sample_metrics.iteritems():
        if key in mapping:
            mapped_metrics[mapping[key]] = value

    return mapped_metrics


def main():
    args = supply_args()
    sample_metrics = {}

    probeqc = args.probeqc_after
    total_cov = calc_cov(probeqc, 'AVGD')
    total_bp = calc_total_bp(probeqc)
    write_me = open(args.outfile, 'w')
    write_me_txt = open(args.outfile_txt, 'w')

    for label in probeqc.headers[5:]:
        this_cov = calc_cov(probeqc, label)
        sample_metrics[label] = calc_metric(total_bp, this_cov)

    # Calculate PUMI metric
    if args.probeqc_before:
        probeqc_before = args.probeqc_before
        total_cov_before = calc_cov(probeqc_before, 'AVGD')
        pumi = calc_metric(total_cov_before, (total_cov * 100))
        sample_metrics['pumi'] = pumi
    else:
        total_cov_before = None

    # Calculate on target reads, or amplicon efficiency
    this_picard = AlignSummaryMetrics(args.picard_summary)
    if total_cov_before:
        on_target = add_on_target(this_picard.metrics, total_cov_before)
    else:
        on_target = add_on_target(this_picard.metrics, total_cov)
    sample_metrics['on_target'] = on_target


    write_to_text(sample_metrics, write_me_txt)
    write_me.write(json.dumps(map_fields(sample_metrics)))
    write_me.close()


if __name__ == "__main__":
    main()
