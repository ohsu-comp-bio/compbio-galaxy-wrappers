#!/usr/bin/env python
"""
Wrapper for xhmm --prepareTargets mode.

usage: xhmm_prepareTargets.py [OPTIONS]...

XHMM version:   https://bitbucket.org/statgen/xhmm 72f7891
Author:         Adam Struck <strucka@ohsu.edu>
Last Modified:  29/09/2015
"""

import argparse
import logging
import subprocess


def collectArgs():
    descr = "xhmm_prepareTargets.py: Sort all target intervals, merge overlapping ones, \
    and print the resulting interval list"
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--targets', help='Input targets lists')
    parser.add_argument('--mergedTargets', help="Output targets list")
    args = parser.parse_args()
    return args


def buildXhmmCmd(args):
    baseCmd = "xhmm --prepareTargets"
    formattedCmd = baseCmd + " --targets {}".format(args.targets)
    if args.mergedTargets is not None:
        formattedCmd = formattedCmd + " --mergedTargets {}".format(
            args.mergedTargets)
    return formattedCmd


def execute(cmd):
    logging.info("RUNNING: %s" % (cmd))
    print "running", cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if len(stderr):
        print stderr
    return p.returncode


def main():
    args = collectArgs()
    if args.targets is None:
        return "--targets option required."
    else:
        cmd = buildXhmmCmd(args)
        execute(cmd)


if __name__ == "__main__":
    main()
