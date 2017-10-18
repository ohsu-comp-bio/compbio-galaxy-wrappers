#!/usr/bin/env python

# USAGE: gather_vcfs.py [<input_vcf>...] <output_vcf>
# CODED BY: John Letaw

import argparse

VERSION = '0.1.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', nargs="+", help='Input files to output in one file.')
    parser.add_argument('output_file', help='Output file.')
    
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def main():

    args = supply_args()

    handle_out = open(args.output_file, 'w')
    first_file = True

    for infile in args.input_file:
        with open(infile, 'ru') as to_concat:
            for line in to_concat:
                if not line.startswith('#'):
                    first_file = False
                    handle_out.write(line)
                elif first_file == True:
                    handle_out.write(line)

    handle_out.close()

if __name__ == "__main__":
    main()


