#!/usr/bin/env python

# USAGE: Galaxy wrapper for Oncotator.
# CODED BY: Jaclyn Smith and John Letaw

from argparse import ArgumentParser
import logging
import argparse
import os
import sys

# Including this because is it quite useful...
# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

VERSION = '1.9.9.0.1'

def oncotator_argparse():
    parser = ArgumentParser(description='Run Oncotator')
    parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: 5]", default=5)
    parser.add_argument('-i', '--input_format', type=str, default="MAFLITE", choices=['VCF', 'SEG_FILE', 'MAFLITE'], help='Input format.  Note that MAFLITE will work for any tsv file with appropriate headers, so long as all of the required headers (or an alias -- see maflite.config) are present.  [default: %s]' % "MAFLITE")
    parser.add_argument('--db-dir', dest='db-dir', help='Main annotation database directory.')
    parser.add_argument('input_file', type=str, help='Input file to be annotated.  Type is specified through options.')
    parser.add_argument('output_file', type=str, help='Output file name of annotated file.')
    parser.add_argument('genome_build', metavar='genome_build', type=str, help="Genome build.  For example: hg19", choices=["hg19"])
    parser.add_argument('--output_format', type=str, default="TCGAMAF", choices=['TCGAMAF', 'VCF', 'SIMPLE_TSV', 'TCGAVCF', 'SIMPLE_BED', 'GENE_LIST'])
    parser.add_argument('--skip-no-alt', dest="skip-no-alt", action='store_true', help="If specified, any mutation with annotation alt_allele_seen of 'False' will not be annotated or rendered.  Do not use if output format is a VCF.  If alt_allele_seen annotation is missing, render the mutation.")
    parser.add_argument('--prepend', dest="prepend", action='store_true', help="If specified for TCGAMAF output, will put a 'i_' in front of fields that are not directly rendered in Oncotator TCGA MAFs")
    parser.add_argument('--infer-onps', dest="infer-onps", action='store_true', help="Will merge adjacent SNPs,DNPs,TNPs,etc if they are in the same sample.  This assumes that the input file is position sorted.  This may cause problems with VCF -> VCF conversion, and does not guarantee input order is maintained.")
    parser.add_argument('-c', '--canonical-tx-file', dest="canonical-tx-file", type=str, help="Simple text file with list of transcript IDs (one per line) to always select where possible for variants.  Transcript IDs must match the ones used by the transcript provider in your datasource (e.g. gencode ENST00000123456).  If more than one transcript can be selected for a variant, uses the method defined by --tx-mode to break ties.  Using this list means that a transcript will be selected from this list first, possibly superseding a best-effect.  Note that transcript version number is not considered, whether included in the list or not.")
    parser.add_argument('--workdir', default='/tmp')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    # Process arguments
    args = parser.parse_args()
    
    return args

def create_opts(arg_dict):
    args = []
    flag_opts = set(('canonical-tx-file', 'skip-no-alt','prepend','infer-onps', 'input_format', 'output_format', 'db-dir'))
    wrapper_arguments = set(('input_file', 'output_file', 'genome_build'))
    for option, value in arg_dict.items():
        if option in wrapper_arguments:
            continue
        if (option in flag_opts):
            if value == True:
                args.append("--" + option)
            elif type(value) != bool and value != None:
                args.append("--" + option + " " + str(value))
    return " ".join(args)


def build_cmd(options, input_file, output_file, build):
    return " ".join([
        "oncotator",
        '-v',
        '--log_name /dev/null',
        input_file,
        output_file,
        build,
        create_opts(options)])

def execute(cmd):

    logging.info("RUNNING: %s" % (cmd))
    print "running", cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    stdout, stderr = p.communicate()
    if len(stderr):
        print stderr
    return p.returncode

if __name__ == "__main__":
    args = oncotator_argparse()
    this_cmd = build_cmd(vars(args), args.input_file, args.output_file, args.genome_build)
    execute(this_cmd) 
