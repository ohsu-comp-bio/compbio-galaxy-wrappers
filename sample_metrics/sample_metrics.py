#!/usr/bin/env python

# DESCRIPTION: Create sample level metrics to be passed to the CGD.  Metrics
#  are passed as a json dump.
# USAGE: sample_metrics.py -h
# CODED BY: John Letaw

from __future__ import print_function

import argparse
import copy
import json
import os
import sys

from file_types.gatk_intervals import ProbeQcRead
from file_types.picard import AlignSummaryMetrics

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

VERSION = '0.5.2'


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
    parser.add_argument('--primers_bed', help='BED file containing primer coordinates only.')
    parser.add_argument('--primers_bam', help='BAM file to calculate primer reads on target.')
    parser.add_argument('--msi', help='TSV file containing MSI results')
    parser.add_argument('--json_in', nargs='*', help='Arbitrary number of files to be included in sample metrics that are in json format.')
    parser.add_argument('--outfile', help='Output file with json string.')
    parser.add_argument('--outfile_new', help='Output file with new style json string.')
    parser.add_argument('--outfile_txt', help='Output file in human readable '
                                              'text format.')
    parser.add_argument('--workflow', help='Pass the Galaxy workflow name, if applicable.')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                                                               VERSION)
    args = parser.parse_args()
    return args


def run_cmd(cmd):
    """
    Run command.
    """
    print('Running the following command:')
    print('\t'.join(cmd))

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return stdout


def get_target_count_cmd(samfile, bed):
    """
    From a BAM file, use samtools view to count the number of reads with targets in primers_bed.
    samtools-1.3.1 view -L primers_only.bed test.bam -c
    Core dumps occur trying to run below pysam command.  Command work great for several versions of samtools on command
    line, so we're doing subprocess until better solution.
    :return:
    """
    #print(pysam.view(samfile, "-L", bed, "-c"))
    cmd = ['samtools', 'view', '-L', bed, '-c', samfile]
    return cmd


def calc_total_bp(probeqc):
    """
    Get the total number of base pairs covered by your targeted region set.
    :return:
    """
    total_bp = 0
    for line in probeqc.probeqc.values():
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
    for line in probeqc.probeqc.values():
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


def add_on_target(picard, total_cov, total_lbl='PF_ALIGNED_BASES'):
    """
    Include the percent on target reads metric, mainly for amplicon assays.
    Also include the on primer frag count percentage.  Rename variables...
    :return:
    """
    if 'PAIR' in picard:
        pf_bases_aligned = int(picard['PAIR'][total_lbl])
    elif 'UNPAIRED' in picard:
        pf_bases_aligned = int(picard['UNPAIRED'][total_lbl])
    else:
        pf_bases_aligned = None

    on_target = str("{:.4}".format((total_cov * 100.0) / pf_bases_aligned))

    return on_target

def write_to_text(sample_metrics, outfile_txt):
    """
    Write metrics to a text file, mainly to be viewed in Galaxy.
    :return:
    """
    for key, value in sorted(sample_metrics.items()):
        outfile_txt.write(key)
        outfile_txt.write(': ')
        outfile_txt.write(value)
        if key != 'AVGD' and key != 'on_primer_frag_count':
            outfile_txt.write('%')
        outfile_txt.write('\n')

    outfile_txt.close()

def req_metrics(wf, metrics):
    """
    Labs only want to see certain metrics...
    :return:
    """
    wf_metric_list = {'QIAseq_V3_RNA': ['Q30', 'AVGD', 'pumi', 'on_primer_frag_count_pct']}
    new_metrics = {}
    if wf in wf_metric_list:
        for metric, value in metrics.items():
            if metric in wf_metric_list[wf]:
                new_metrics[metric] = value
    else:
        return metrics

    return new_metrics

def map_fields(value):
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
    mapping = {'Q30': 'qthirty',
                      'D10': 'depthTen',
                      'D20': 'depthTwenty',
                      'D50': 'depthFifty',
                      'D100': 'depthOneHundred',
                      'D250': 'depthTwoHundredFifty',
                      'D500': 'depthFiveHundred',
                      'D700': 'depthSevenHundred',
                      'D1000': 'depthOneThousand',
                      'D1250': 'depthTwelveHundredFifty',
                      'D2000': 'depthTwoThousand',
                      'AVGD': 'averageDepth',
                      'pumi': 'percentUmi',
                      'on_target': 'percentOnTarget',
                      'on_primer_frag_count': 'total_on_target_transcripts',
                      'on_primer_frag_count_pct': 'total_on_target_transcripts_pct'}
    if value in mapping:
        return mapping[value]
    else:
        return value

def iter_samp_met(sample_metrics):
    """
    Iterate through sample metrics and return the mapping.
    :return:
    """
    mapped_metrics = {}
    for key, value in sample_metrics.items():
        mapped_metrics[map_fields(key)] = value
    return mapped_metrics

def new_iter_samp_met(sample_metrics):
    """
    Iterate through new sample metrics and return the mapping.
    :param sample_metrics:
    :return:
    """
    for met_type, values in sample_metrics.items():
        for metric in values:
            metric["metric"] = map_fields(metric["metric"])
    return sample_metrics

def create_sample_metrics(args):
    """

    :return:
    """
    sample_metrics = {}
    probeqc = args.probeqc_after
    total_cov = calc_cov(probeqc, 'AVGD')
    total_bp = calc_total_bp(probeqc)
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
    # if total_cov_before:
    #     on_target = add_on_target(this_picard.metrics, total_cov)
    # else:
    on_target = add_on_target(this_picard.metrics, total_cov)
    sample_metrics['on_target'] = on_target

    # # Only doing the below for new sample metrics until further notice.
    # if args.primers_bed and args.primers_bam:
    #     on_primer_frag_count = run_cmd(get_target_count_cmd(args.primers_bam, args.primers_bed))
    #     sample_metrics['on_primer_frag_count'] = on_primer_frag_count.rstrip('\n')
    return sample_metrics


def create_new_sample_metrics(args):
    """
    Same as above, but this time with the new format:
    {
    "sampleRunMetrics": [
        {
            "metric": "total_on_target_reads",
            "value": 1230411
        },
        {
            "metric": "percent_on_target_reads",
            "value": 0.91
        }
    ],
    "geneMetrics": [
        {
                "gene": "ASXL1",
            "metric": "total_on_target_reads",
            "value": 469012
        },
        {
                "gene": "BRCA1",
            "metric": "total_on_target_reads",
            "value": 362330
        }
    ]
}

    TODO: Simplify/refactor this using decorators.
    :param args:
    :return:
    """
    sample_metrics = {"sampleRunMetrics": [], "geneMetrics": []}
    # # Right now, they just want the new metrics to be passed to the new endpoint.  This is fine since I don't like
    # # the way I'm doing this right now.  Stuff below commented out until future notice.
    # probeqc = args.probeqc_after
    # total_cov = calc_cov(probeqc, 'AVGD')
    # total_bp = calc_total_bp(probeqc)
    # for label in probeqc.headers[5:]:
    #     this_cov = calc_cov(probeqc, label)
    #     sample_metrics["sampleRunMetrics"].append({"metric": label, "value": calc_metric(total_bp, this_cov)})
    # # Calculate PUMI metric
    # if args.probeqc_before:
    #     probeqc_before = args.probeqc_before
    #     total_cov_before = calc_cov(probeqc_before, 'AVGD')
    #     pumi = calc_metric(total_cov_before, (total_cov * 100))
    #     sample_metrics['sampleRunMetrics'].append({"metric": "pumi", "value": pumi})
    # else:
    #     total_cov_before = None
    #
    # # Calculate on target reads, or amplicon efficiency
    this_picard = AlignSummaryMetrics(args.picard_summary)
    # if total_cov_before:
    #     on_target = add_on_target(this_picard.metrics, total_cov)
    # else:
    #     on_target = add_on_target(this_picard.metrics, total_cov)
    # sample_metrics["sampleRunMetrics"].append({"metric": "on_target", "value": on_target})

    if args.msi:
        msi_res = msi_parse(args.msi)
        sample_metrics['sampleRunMetrics'].append({"metric": "msi_sites", "value": msi_res[0]})
        sample_metrics['sampleRunMetrics'].append({"metric": "msi_somatic_sites", "value": msi_res[1]})
        sample_metrics['sampleRunMetrics'].append({"metric": "msi_pct", "value": msi_res[2]})
    if args.primers_bed and args.primers_bam:
        on_primer_frag_count = run_cmd(get_target_count_cmd(args.primers_bam, args.primers_bed))
        sample_metrics["sampleRunMetrics"].append({"metric": "on_primer_frag_count", "value": on_primer_frag_count.rstrip('\n')})
        on_primer_frag_count_pct = add_on_target(this_picard.metrics, int(on_primer_frag_count), 'PF_HQ_ALIGNED_READS')
        sample_metrics["sampleRunMetrics"].append({"metric": "on_primer_frag_count_pct", "value": on_primer_frag_count_pct.rstrip('\n')})
    if args.json_in:
        for filename in args.json_in:
            with open(filename, 'rU') as myfile:
                for line in myfile:
                    for k, v in json.loads(line).items():
                        sample_metrics['sampleRunMetrics'].append({"metric": str(k), "value": v})

    return sample_metrics

def msi_parse(filename):
    """
    Get the MSI results from the TSV file.
    Total_Number_of_Sites	Number_of_Somatic_Sites	%
    230	60	26.09

    :return:
    """
    with open(filename, 'rU') as mymsi:
        for line in mymsi:
            if not line.startswith('Total'):
                line = line.rstrip('\n').split('\t')
                return (int(line[0]), int(line[1]), float(line[2]))


def main():
    args = supply_args()
    sample_metrics = req_metrics(args.workflow, create_sample_metrics(args))
    new_sample_metrics = create_new_sample_metrics(args)

    # Write the output.
    write_me = open(args.outfile, 'w')
    write_me_new = open(args.outfile_new, 'w')
    write_me_txt = open(args.outfile_txt, 'w')
    write_to_text(sample_metrics, write_me_txt)
    write_me.write(json.dumps(iter_samp_met(sample_metrics)))
    write_me_new.write(json.dumps(new_iter_samp_met(new_sample_metrics)))
    write_me.close()
    write_me_new.close()


if __name__ == "__main__":
    main()
