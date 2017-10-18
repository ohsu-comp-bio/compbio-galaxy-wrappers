#!/usr/bin/env python
"""
Wrapper for xhmm --mergeGATKdepths mode.

usage: xhmm_mergeGATKdepths.py [OPTIONS]...

XHMM version:   https://bitbucket.org/statgen/xhmm 72f7891
Author:         Adam Struck <strucka@ohsu.edu>
Last Modified:  29/09/2015
"""

import argparse
import logging
import subprocess


def collectArgs():
    descr = 'xhmm_mergeGATKdepths.py: Sort all target intervals, merge overlapping ones, \
    and print the resulting interval list'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--GATKdepths', help='GATK sample_interval_summary output file(s) to \
    be merged [must have *IDENTICAL* target lists]')
    parser.add_argument('--GATKdepthsList', help='A file containing a list of GATK \
    sample_interval_summary output files to be merged [must have *IDENTICAL* \
    target lists]')
    parser.add_argument('--sampleIDmap', help='File containing mappings of sample \
    names to new sample names (in columns designated by fromID, toID)')
    parser.add_argument('--fromID', type=int, default=1, help='Column number of OLD sample IDs \
    to map')
    parser.add_argument('--toID', type=int, default=2, help='Column number of NEW sample IDs \
    to map')
    parser.add_argument('--columnSuffix', default='_mean_cvg', help='Suffix of columns to be \
    used for merging [where columns are in the form: SAMPLE + columnSufix]')
    parser.add_argument('--rdPrecision', type=int, default=2, help='Decimal precision of \
    read depths output')
    parser.add_argument('--outputTargetsBySamples', help='Output targets x samples  \
    (instead of samples x targets)')
    args = parser.parse_args()
    return args


def buildXhmmCmd(args):
    baseCmd = 'xhmm --mergeGATKdepths'
    formattedCmd = baseCmd + ' --GATKdepths {}'.format(args.GATKdepths)
    return formattedCmd


def execute(cmd):
    logging.info('RUNNING: %s' % (cmd))
    print 'running', cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if len(stderr):
        print stderr
    return p.returncode


def main():
    args = collectArgs()
    print vars(args)
    if args.GATKdepths is None:
        return '--GATKdepths option required.'
    else:
        cmd = buildXhmmCmd(args)
        print cmd
        # execute(cmd)


if __name__ == '__main__':
    main()
