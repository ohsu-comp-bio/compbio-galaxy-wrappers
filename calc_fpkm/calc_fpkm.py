#!/usr/bin/env python

# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
import argparse
import numpy
import os
import sys

# Including this because is it quite useful...
# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess


VERSION = '0.2.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('infile', help='Input counts file from htseq-count.')    
    parser.add_argument('outfile_fpkm', help='Output FPKM normalized counts '
                                             'file.')
    parser.add_argument('outfile_fpkm_uq', help='Output FPKM-UQ normalized '
                                                'counts file.')
    parser.add_argument('outfile_tpm', help='Output TPM normalized counts '
                                            'file.')
    parser.add_argument('gtf', help='Input GTF file.')
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
    print("From Output: " + stdout)
    eprint("From Error: " + stderr)

    return p.wait()


def eprint(*args, **kwargs):
    """
    http://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)


class GTFReader(object):
    """
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1";  level 2; havana_gene "OTTHUMG00000000961.2";

    For this particular GTF, we want gene and Selenocysteine gene_id entries.  These are always the first entry of the last column
    """
    def __init__(self, filename):
        self.filename = filename
        self.gene_ids = self.find_gene_id()
        self.gene_coords = self.find_gene_coords()


    def find_gene_coords(self):
        """
        """
        gene_coords = {}
        with open(self.filename, 'rU') as my_gtf:
            for line in my_gtf:
                if not line.startswith('#'):
                    line = line.rstrip('\n').split('\t')
                    gene_type = line[2]
                    if gene_type == "exon" or gene_type == "start_codon" or gene_type == "stop_codon" or gene_type == "UTR":
                        gene_id = line[8].split(';')[0].split(' ')[1][1:-1]
                        if gene_id not in gene_coords:
                            gene_coords[gene_id] = []
                            
                        start = int(line[3])
                        stop = int(line[4])
                        for i in range(start, stop):
                            gene_coords[gene_id].append(i)
                        
        return gene_coords

    
    def find_gene_id(self):
        """
        """
        gene_ids = []
        with open(self.filename, 'rU') as my_gtf:
            for line in my_gtf:
                if not line.startswith('#'):
                    line = line.rstrip('\n').split('\t')
                    gene_type = line[2]
                    if gene_type == 'gene' or gene_type == 'Selenocysteine':
                        gene_id = line[8].split(';')[0].split(' ')[1][1:-1]
                        if gene_id not in gene_ids:
                            gene_ids.append(gene_id)
        return gene_ids
    

class HtseqReader(object):
    """
    """
    def __init__(self, filename, gtf):
        self.filename = filename
        self.gtf = gtf
        self.total, self.rpk_total = self.count_total()
        self.counts = self.make_counts(None)
        self.np_arr = numpy.array(self.counts.values())
        self.total_uq = numpy.percentile(self.np_arr, 75)
        self.counts_fpkm = self.make_counts('fpkm')
        self.counts_fpkm_uq = self.make_counts('fpkm_uq')
        self.counts_tpm = self.make_counts('tpm')

    def count_total(self):
        """
        Total number of counts in the HTSeq table of counts.
        Also count the RPK total, for TPM calculation.
        """
        total = 0
        rpk_total = 0.0
        with open(self.filename, 'rU') as my_htseq:
            for line in my_htseq:
                if line[0] != '_':
                    line = line.rstrip('\n').split('\t')
                    ensg_id = line[0]
                    gene_len = len(set(self.gtf.gene_coords[ensg_id])) / 1000.0
                    count = int(line[1])
                    total += count
                    rpk_total += float(count/gene_len)
        return total, rpk_total
    
        
    def make_counts(self, count_type):
        """
        """
        counts_dict = {}
        with open(self.filename, 'rU') as my_htseq:
            for line in my_htseq:
                if line[0] != '_':
                    line = line.rstrip('\n').split('\t')
                    ensg_id = line[0]
                    count = float(line[1])

                    if count_type == 'fpkm':
                        count = self.calc_fpkm(ensg_id, count)
                    elif count_type == 'fpkm_uq':
                        count = self.calc_fpkm_uq(ensg_id, count)
                    elif count_type == 'tpm':
                        count = self.calc_tpm(ensg_id, count)

                    counts_dict[ensg_id] = count

        return counts_dict

    
    def calc_fpkm(self, ensg_id, count):
        """
        """
        gene_len = len(set(self.gtf.gene_coords[ensg_id]))
        return float((count*1000000000.0)/(self.total*gene_len))

    
    def calc_fpkm_uq(self, ensg_id, count):
        """
        """
        gene_len = len(set(self.gtf.gene_coords[ensg_id]))
        return float((count*1000000000.0)/(self.total_uq*gene_len))


    def calc_tpm(self, ensg_id, count):
        """
        Calculate Transcripts per Million normalized counts.
        :return:
        """
        gene_len = len(set(self.gtf.gene_coords[ensg_id])) /1000.0
        return float((count/gene_len)/(self.rpk_total/1000000.0))


class HtseqWriter(object):
    """
    """
    def __init__(self, htseq, filename):
        self.htseq = htseq
        self.filename = open(filename, 'w')
        self.write_file()
        
    def write_file(self):
        """
        """
        for key, value in sorted(self.htseq.iteritems()):
            value = '{:.3f}'.format(value)
            to_write = '\t'.join([key, value])
            self.filename.write(to_write)
            self.filename.write('\n')

        self.filename.close()

        
def main():

    args = supply_args()
    my_gtf = GTFReader(args.gtf)

    my_htseq = HtseqReader(args.infile, my_gtf)

    HtseqWriter(my_htseq.counts_fpkm, args.outfile_fpkm)
    HtseqWriter(my_htseq.counts_fpkm_uq, args.outfile_fpkm_uq)
    HtseqWriter(my_htseq.counts_tpm, args.outfile_tpm)

if __name__ == "__main__":
    main()


