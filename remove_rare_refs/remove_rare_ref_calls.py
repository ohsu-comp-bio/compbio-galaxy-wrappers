#!/usr/bin/env python3

import argparse
from collections import OrderedDict

VERSION = '0.2.0'

def supply_args():
    """                                                                                                                        
    Populate arguments.
    """
    parser = argparse.ArgumentParser(description='SeattleSeq Formatter')
    parser.add_argument(dest='input', help='Input SeattleSeq TSV.')
    parser.add_argument(dest='output', help='Output SeattleSeq TSV.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class SeattleSeqEntry():
    def __init__(self, line):
        self.line = line
        self.chrom = self.line[0]
        self.pos = self.line[1]
        self.ref = self.line[3]
        self.alt = self.line[4]
        self.depth = self.line[8]
        self.uniq_key = (self.chrom, self.pos, self.ref, self.alt)


class SeattleSeqReader():
    def __init__(self, filename):
        self.filename = filename
        self.sseq = self._parse_sseq()

    def _parse_sseq(self):
        """
        Prepare the structure that will contain SeattleSeq WG data.
        """
        sseq = OrderedDict()
        with open(self.filename, 'r') as infile:
            for line in infile:
                line = line.rstrip('\n').split('\t')
                if line[0].startswith('#'):
                    self.header = line
                else:
                    rec = SeattleSeqEntry(line)
                    if rec.uniq_key not in sseq:
                        sseq[rec.uniq_key] = rec
        return sseq


class SeattleSeqWriter():
    def __init__(self, filename):
        self.outfile = open(filename, 'w')

    def write_line(self, line):
        self.outfile.write('\t'.join(line))
        self.outfile.write('\n')


def main():
    
    args = supply_args()
    rare_refs = [("1","169519049"),("7","6026775"),("13","32929387"),("14","75513883")]
    fake_depth = "100(100,0)"

    my_sseq = SeattleSeqReader(args.input)
    writer = SeattleSeqWriter(args.output)
    writer.write_line(my_sseq.header)
    for entry in my_sseq.sseq.values():
        if (entry.chrom, entry.pos) not in rare_refs:
            if entry.depth != fake_depth:
                if len(entry.alt) < 256:
                    if len(entry.ref) < 256:
                        writer.write_line(entry.line)
    writer.outfile.close()

if __name__ == "__main__":
    main()
