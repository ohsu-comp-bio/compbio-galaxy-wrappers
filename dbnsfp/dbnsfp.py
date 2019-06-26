#!/usr/bin/env python

# USAGE: dbnsfp.py <input> <output_1> <output_2> <dbnsfp_db_path>
# CODED BY: John Letaw

import argparse
import os
import subprocess

PROGRAM_VERSION = '3.5c'
VERSION = '0.1.3'


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', help='')
    parser.add_argument('output_dbnsfp', help='')
    parser.add_argument('output_dbscsnv', help='')
    parser.add_argument('dbnsfp_path', help='')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + VERSION)
    parser.add_argument('-p', '--program_version', action='version', version='dbNSFP ' + PROGRAM_VERSION)
    args = parser.parse_args()
    return args


def create_cmd(args, input_vcf):
    """
    Create command to be run.
    """
    cmd = ['java', remove_class_name(args.dbnsfp_path), '-i', input_vcf, '-o', args.output_dbnsfp, '-v', 'hg19', '-s', args.output_dbscsnv, '-p']
    print(cmd)
    return cmd


def remove_class_name(path):
    """
    We are being sent PATH plus binary name.  Strip off the PATH.
    """
    return path.split('/')[-1]


def remove_prog_path(path):
    """
    We are being sent PATH plus binary name.  Strip off the binary name.
    """
    return '/'.join(path.split('/')[:-1])


def run_cmd(cmd):
    """
    Run the command built in create_cmd().
    """
    my_proc = subprocess.call(cmd)
    return my_proc


def main():

    args = supply_args()
    cwd = os.getcwd()
    symfile = '/'.join([cwd, 'input.vcf'])
    os.symlink(args.input_vcf, symfile)
    cmd = create_cmd(args, symfile)
    with cd(remove_prog_path(args.dbnsfp_path)):
        try:
            run_cmd(cmd)
            os.remove(symfile)
        except:
            os.remove(symfile)


if __name__ == "__main__":
    main()


