#!/usr/bin/env python

# USAGE: exomiser.py <vcf> <outfile> <java8_path> <exomiser_jar> [optional args]
# CODED BY: John Letaw

from __future__ import print_function

import argparse
import os
import sys

# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

VERSION = '7.2.3.3'

def supply_args():
    """                                                                                                                        
    Populate arguments.
    """
    parser = argparse.ArgumentParser(description='A Tool to Annotate and Prioritize Exome Variants')

    # Required
    parser.add_argument('vcf', help='')
    parser.add_argument('outfile_html', help='HTML')
    parser.add_argument('outfile_vcf', help='VCF')
    parser.add_argument('outfile_tab_gene', help='TAB-GENE')
    parser.add_argument('outfile_tab_variant', help='TAB-VARIANT')
    parser.add_argument('java8_path', help='Defined in $JAVA8_PATH.  Look in Galaxy environmental variable script.')
    parser.add_argument('exomiser_jar', help='Defined in $EXOMISER_JAR.  Look in Galaxy environmental variable script.')

    # Input options
    parser.add_argument('--prioritiser', choices=['hiphive', 'phive', 'phenix', 'exomewalker'], help='')
    parser.add_argument('-I', '--inheritance-mode', choices=['AR', 'AD', 'X'], help='')
    parser.add_argument('--hpo-ids', help='')
    parser.add_argument('--settings-file', help='')
    parser.add_argument('-R', '--restrict-interval', help='')
    parser.add_argument('-p', '--ped', help='')
    parser.add_argument('--full-analysis', action='store_true', help='')

    # Output options
    parser.add_argument('-f', '--out-format', default='HTML', help='')
    parser.add_argument('--num-genes', help='')
    parser.add_argument('--output-pass-variants-only', action='store_true', help='')

    # Filtering options
    parser.add_argument('-T', '--keep-off-target', action='store_true', help='')
    parser.add_argument('-P', '--keep-non-pathogenic', action='store_true', help='')
    parser.add_argument('--remove-known-variants', action='store_true', help='')
    parser.add_argument('-F', '--max-freq', metavar='max-freq', help='')
    parser.add_argument('-Q', '--min-qual', metavar='min-qual', help='')

    # Pedigree options
    parser.add_argument('--proband_id', help='')
    parser.add_argument('--proband_sex', choices=['male', 'female'], help='')
    parser.add_argument('--mother_id', help='')
    parser.add_argument('--father_id', help='')
    parser.add_argument('--proband_aff', action='store_true', help='')
    parser.add_argument('--mother_aff', action='store_true', help='')
    parser.add_argument('--father_aff', action='store_true', help='')

    # Advanced options
    parser.add_argument('--candidate-gene', help='')
    parser.add_argument('-D', '--disease-id', help='')
    parser.add_argument('--genes-to-keep', help='')
    parser.add_argument('-S', '--seed-genes', help='')
    parser.add_argument('-E', '--hiphive-params', help='')
    
    # Wrapper version information.
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION + '\n')
    args = parser.parse_args()
    return args


def build_cmd(args, pedigree):
    """
    Build the command that will be sent via subprocess.
    Create based on args definition.
    """
    dont_write = ('java8_path', 'exomiser_jar', 'outfile_html', 'outfile_vcf', 
                  'outfile_tab_gene', 'outfile_tab_variant', 'ped', 'proband_id', 
                  'mother_id', 'father_id', 'proband_sex', 'proband_aff', 
                  'father_aff', 'mother_aff')
    cmd = [args.java8_path, '-jar', args.exomiser_jar]

    for flag, arg, in vars(args).iteritems():
        form_flag = '--{0}'.format(flag.replace('_', '-'))
        if arg is not None and not isinstance(arg, bool) and flag not in dont_write:
            cmd.extend([form_flag, arg])
        elif arg == True and flag not in dont_write:
            cmd.extend([form_flag])

    if pedigree:
        cmd.extend(['--ped', 'pedigree'])
    elif args.ped:
        cmd.extend(['--ped', args.ped])
    else:
        pass

    return cmd


def run_cmd(cmd):
    """
    Run command.
    """
    print('Running the following command:')
    print('\t'.join(cmd))

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()

    print("From Exomiser Output: " + stdout)

    if p.wait() == 0:
        print("From Exomiser Error: " + stderr)
    else:
        eprint("From Exomiser Error: " + stderr)

    return p.wait()


def rename_output(args, outfile, res_dir, out_ext):
    """
    Galaxy will provide a .dat output file, while Exomiser creates files using
    user-defined prefixes.
    """
    just_file = args.vcf.split('/')[-1]
    exomiser_file = res_dir + '/' + just_file + '-exomiser-results.' + out_ext
    os.rename(exomiser_file, outfile)


def mkdir_p(path):
    """
    Check if PATH exists, mkdir if it doesn't.
    """
    if not os.path.isdir(path):
        os.mkdir(path)


def assign_aff(aff_bool):
    """
    Take a True or False value, return a 2 or 1.
    """
    if aff_bool == True:
        return '2'
    elif aff_bool == False:
        return '1'
    return None


def build_ped(args):
    """
    Build a pedigree file from inputs.
    """
    fam = 'FAM1'

    mother_aff = assign_aff(args.mother_aff)
    father_aff = assign_aff(args.father_aff)
    proband_aff = assign_aff(args.proband_aff)

    if args.proband_id and args.mother_id and args.father_id and args.proband_sex:
        handle_ped = open('pedigree', 'w')
        linem = [fam, args.mother_id, '0', '0', '2', mother_aff]
        linef = [fam, args.father_id, '0', '0', '1', father_aff]
        if args.proband_sex == 'male':
            linep = [fam, args.proband_id, args.father_id, args.mother_id, '1', proband_aff]
        else:
            linep = [fam, args.proband_id, args.father_id, args.mother_id, '2', proband_aff]

        handle_ped.write('\t'.join(linem))
        handle_ped.write('\n')
        handle_ped.write('\t'.join(linef))
        handle_ped.write('\n')
        handle_ped.write('\t'.join(linep))
        handle_ped.write('\n')

        handle_ped.close()
        return True
    
    return False


def eprint(*args, **kwargs):
    """
    http://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)
    

def main():

    # I don't know of a way to change the value of this in Exomiser, so it is hard-coded.
    res_dir = 'results'

    args = supply_args()
    pedigree = build_ped(args)
    cmd = build_cmd(args, pedigree)
    # For some reason this directory is not automatically created by Exomiser.
    mkdir_p(res_dir)
    returncode = run_cmd(cmd)
    rename_output(args, args.outfile_html, res_dir, 'html')
    rename_output(args, args.outfile_vcf, res_dir, 'vcf')
    rename_output(args, args.outfile_tab_gene, res_dir, 'genes.tsv')
    rename_output(args, args.outfile_tab_variant, res_dir, 'variants.tsv')
    os.rmdir(res_dir)

if __name__ == "__main__":
    main()


