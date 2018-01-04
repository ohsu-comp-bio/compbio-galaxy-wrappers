#!/usr/bin/env python

# DESCRIPTION: Remove duplicate columns in the SeattleSeq tsv output.
# USAGE: dedup_seattleseq_header.py ...
# CODED BY: John Letaw

import argparse

VERSION = '0.1.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('seattleseq', help="SeattleSeq input tsv file")
    parser.add_argument('outfile', help="Output corrected TSV")
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def main():

    args = supply_args()
    handle_out = open(args.outfile, 'w')
    dedup = False
    check = False
    with open(args.seattleseq, 'rU') as infile:
        for line in infile:
            sline = line.rstrip('\n').split('\t')
            print(sum(y is not -1 for y in [x.find('Gtype') for x in sline]))
            print(len([x.find('Qual') for x in sline]))
            if sum(y is not -1 for y in [x.find('Gtype') for x in sline]) == 2 and \
               sum(y is not -1 for y in [x.find('Qual') for x in sline]) == 2 and \
               not check:
                dedup = True
            check = True
            if dedup:
                nline = sline[:9]
                nline.append(sline[10])
                nline.append(sline[12])
                nline.extend(sline[14:])
                handle_out.write('\t'.join(nline))
                handle_out.write('\n')
            else:
                handle_out.write(line)

    handle_out.close()

if __name__ == "__main__":
    main()

