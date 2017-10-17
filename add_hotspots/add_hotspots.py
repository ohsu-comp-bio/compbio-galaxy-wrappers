#!/usr/bin/env python

# DESCRIPTION:
# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
import argparse
import os
import sys
sys.path.append('/home/groups/clinical/users/letaw/jhl_tools')

# Including this because is it quite useful...
# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess


VERSION = '0.1.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('vcf', help='')
    parser.add_argument('hotspots', help='')
    parser.add_argument('outfile', help='')

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
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
    print("From Output: " + stdout)
    eprint("From Error: " + stderr)

    return p.wait()


def eprint(*args, **kwargs):
    """
    http://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)


def replace_filter(filt, anno):
    """
    Remove the dot if it's there, otherwise just append something with a
    semicolon.
    :return:
    """
    if filt == '.' or filt == 'PASS':
        return anno
    else:
        return ';'.join([filt, anno])


def main():

    args = supply_args()
    handle_out = open(args.outfile, 'w')

    all_vcf = {}
    with open(args.vcf, 'rU') as vcf:
        for line in vcf:
            if not line.startswith('#'):
                nline = line.rstrip('\n').split('\t')
                uniq_key = (nline[0], nline[1], nline[3], nline[4])
                if uniq_key not in all_vcf:
                    all_vcf[uniq_key] = line
                handle_out.write(line)
            else:
                handle_out.write(line)

    with open(args.hotspots, 'rU') as hots:
        for line in hots:
            if not line.startswith('#'):
                nline = line.rstrip('\n').split('\t')
                uniq_key = (nline[0], nline[1], nline[3], nline[4])
                if uniq_key not in all_vcf:
                    all_vcf[uniq_key] = line
                    nline[6] = replace_filter(nline[6], 'hotspot')
                    if not 'DP=0;' in nline[7]:
                        handle_out.write('\t'.join(nline))
                        handle_out.write('\n')


    handle_out.close()

if __name__ == "__main__":
    main()
