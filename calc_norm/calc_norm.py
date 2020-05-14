#!/usr/bin/env python3

from __future__ import print_function
import argparse
import numpy
import os
import sys

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

VERSION = '0.6.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='Takes an input file from htseq with ENST or ENSG IDs '
                                                 'and outputs 3 files: FPKM, '
                                                 'FPKM with upper quartile normalization,'
                                                 'and a pseudo-TPM file (using the gene/transcript length)')
    parser.add_argument('input', help='Input counts file from htseq-count.')
    parser.add_argument('gtf', help='Input GTF file.')
    parser.add_argument('output_fpkm', help='Output FPKM normalized counts file.')
    parser.add_argument('output_fpkm_uq', help='Output FPKM-UQ normalized counts file.')
    parser.add_argument('output_tpm', help='Output TPM normalized counts file.')
    parser.add_argument('--use_star_counts', action='store_true', help='Use gene counts from STAR.')
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
    print("From Output: " + str(stdout))
    eprint("From Error: " + str(stderr))

    return p.wait()


def eprint(*args, **kwargs):
    """
    http://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)


class GTFReader(object):
    """
    Input: GTF version
        #!genome-build GRCh37.p13
        #!genome-version GRCh37
        #!genome-date 2009-02
        #!genome-build-accession NCBI:GCA_000001405.14
        #!genebuild-last-updated 2013-09
        1	pseudogene	gene	11869	14412	.	+	.	gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source
         "ensembl_havana"; gene_biotype "pseudogene";
    For this particular GTF, we want gene and Selenocysteine gene_id entries.
    These are always the first entry of the last column
    """

    def __init__(self, filename):
        self.filename = filename
        self.gene_ids = self.find_gene_id()
        self.gene_coords = self.find_gene_coords()
        self.tx_coords = self.find_tx_coords()

    def find_gene_coords(self):
        """
        """
        gene_coords = {}
        with open(self.filename, 'rU') as my_gtf:
            for line in my_gtf:
                if not line.startswith('#'):
                    line = line.rstrip('\n').split('\t')
                    gene_type = line[2]
                    if (gene_type == "exon" or
                            gene_type == "start_codon" or
                            gene_type == "stop_codon" or
                            gene_type == "UTR"):
                        # dictionary of the column 9 because they can be in different order
                        dattr = {}
                        for item in [x for x in line[8].split(';') if x]:
                            attr = item.strip().split(' ')
                            dattr[attr[0]] = attr[1].strip('\"')
                        gene_id = str(dattr['gene_id'])
                        if gene_id not in gene_coords:
                            gene_coords[gene_id] = []
                        start = int(line[3])
                        stop = int(line[4])+1
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

    def find_tx_coords(self):
        """
        Return: ENST id coordinates
        """
        tx_coords = {}
        with open(self.filename, 'rU') as my_gtf:
            for line in my_gtf:
                if not line.startswith('#'):
                    line = line.rstrip('\n').split('\t')
                    gene_type = line[2]
                    if (gene_type == "exon" or
                            gene_type == "start_codon" or
                            gene_type == "stop_codon" or
                            gene_type == "UTR"):
                        tx_id = line[8].split(';')[1].split(' ')[2][1:-1]
                        if tx_id not in tx_coords:
                            tx_coords[tx_id] = []

                        start = int(line[3])
                        stop = int(line[4])
                        for i in range(start, stop):
                            tx_coords[tx_id].append(i)

        return tx_coords


class HtseqReader(object):
    """
    """

    def __init__(self, filename, gtf):
        self.filename = filename
        self.gtf = gtf
        self.total, self.rpk_total = self.count_total()
        self.counts = self.make_counts(count_type=None)

        try:
            self.counts_fpkm = self.make_counts('fpkm')
            self.counts_fpkm_uq = self.make_counts('fpkm_uq')
            self.counts_tpm = self.make_counts('tpm')
        except:
            sys.exit("All counts are zero")

    def coords(self, ident):
        """

        :param ident:
        :return:
        """
        if ident.startswith(("ENST", "rna")):
            coords = self.gtf.tx_coords[ident]
        else:
            coords = self.gtf.gene_coords[ident]
        return coords

    def count_total(self):
        """
        Total number of counts in the HTSeq table of counts.
        Also count the RPK total, for TPM calculation.
        """
        total = 0
        rpk_total = 0.0
        with open(self.filename, 'rU') as my_htseq:
            for line in my_htseq:
                if '_' not in line:
                    line = line.rstrip('\n').split('\t')
                    ident = line[0]
                    gene_len = len(set(self.coords(ident))) / 1000.0
                    count = int(line[1])
                    total += count
                    rpk_total += float(count / gene_len)
        return total, rpk_total

    def make_counts(self, count_type):
        """
        """
        counts_dict = {}
        with open(self.filename, 'rU') as my_htseq:
            for line in my_htseq:
                if '_' not in line:
                    line = line.rstrip('\n').split('\t')
                    ident = line[0]
                    count = float(line[1])

                    if count_type == 'fpkm':
                        count = self.calc_fpkm(ident, count)
                    elif count_type == 'fpkm_uq':
                        count = self.calc_fpkm_uq(ident, count)
                    elif count_type == 'tpm':
                        count = self.calc_tpm(ident, count)

                    counts_dict[id] = count

        return counts_dict

    def calc_fpkm(self, ident, count):
        """
        """
        gene_len = len(set(self.coords(ident)))
        return float((count * 1000000000.0) / (self.total * gene_len))

    def calc_fpkm_uq(self, ident, count):
        """
        """
        np_arr = numpy.array(self.counts.values())
        # upper quartile without zeroes
        total_uq = numpy.percentile(np_arr[np_arr > 0], 75)
        gene_len = len(set(self.coords(ident)))
        return float((count * 1000000000.0) / (total_uq * gene_len))

    def calc_tpm(self, ident, count):
        """
        Calculate Transcripts per Million normalized counts.
        :return:
        """
        gene_len = len(set(self.coords(ident))) / 1000.0
        return float((count / gene_len) / (self.rpk_total / 1000000.0))


class StarGeneCountReader(object):
    """
    TODO: This is mostly the same as HTSeq Reader.  Create superclass Counter, then these can go underneath.
    count_total and make_counts are different.
    """

    def __init__(self, filename, gtf):
        self.filename = filename
        self.gtf = gtf
        self.total, self.rpk_total = self.count_total()
        self.counts = self.make_counts(count_type=None)

        try:
            self.counts_fpkm = self.make_counts('fpkm')
            self.counts_fpkm_uq = self.make_counts('fpkm_uq')
            self.counts_tpm = self.make_counts('tpm')
        except:
            sys.exit("All counts are zero")

    def coords(self, ident):
        if id.startswith(("ENST", "rna")):
            coords = self.gtf.tx_coords[ident]
        else:
            coords = self.gtf.gene_coords[ident]
        return coords

    def count_total(self):
        """
        Total number of counts in the HTSeq table of counts.
        Also count the RPK total, for TPM calculation.
        """
        total = 0
        rpk_total = 0.0
        with open(self.filename, 'rU') as my_htseq:
            for line in my_htseq:
                if not line.startswith('N_'):
                    line = line.rstrip('\n').split('\t')
                    id = line[0]
                    gene_len = len(set(self.coords(id))) / 1000.0
                    count = int(line[3])
                    total += count
                    rpk_total += float(count / gene_len)
        return total, rpk_total

    def make_counts(self, count_type):
        """
        line[3] corresponds to the "reverse" method of counting, as htseq-count describes it.
        Shoould probably generalize this and add a parameter that lets you pick a counts column.
        """
        counts_dict = {}
        with open(self.filename, 'rU') as my_htseq:
            for line in my_htseq:
                if not line.startswith('N_'):
                    line = line.rstrip('\n').split('\t')
                    id = line[0]
                    count = float(line[3])

                    if count_type == 'fpkm':
                        count = self.calc_fpkm(id, count)
                    elif count_type == 'fpkm_uq':
                        count = self.calc_fpkm_uq(id, count)
                    elif count_type == 'tpm':
                        count = self.calc_tpm(id, count)

                    counts_dict[id] = count

        return counts_dict

    def calc_fpkm(self, id, count):
        """
        """
        gene_len = len(set(self.coords(id)))
        return float((count * 1000000000.0) / (self.total * gene_len))

    def calc_fpkm_uq(self, id, count):
        """
        """
        np_arr = numpy.array(self.counts.values())
        # upper quartile without zeroes
        total_uq = numpy.percentile(np_arr[np_arr > 0], 75)
        gene_len = len(set(self.coords(id)))
        return float((count * 1000000000.0) / (total_uq * gene_len))

    def calc_tpm(self, id, count):
        """
        Calculate Transcripts per Million normalized counts.
        :return:
        """
        gene_len = len(set(self.coords(id))) / 1000.0
        return float((count / gene_len) / (self.rpk_total / 1000000.0))


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
        for key, value in sorted(self.htseq.items()):
            value = '{:.3f}'.format(value)
            to_write = '\t'.join([key, value])
            self.filename.write(to_write)
            self.filename.write('\n')

        self.filename.close()


class StarReader(HtseqReader):
    """
    first 4 rows:
    N_unmapped      17266161        17266161        17266161
    N_multimapping  16713652        16713652        16713652
    N_noFeature     1093465 114459752       3227050
    N_ambiguous     8750904 263962  4478729

    column headings:
    column 1: gene ID
    column 2: counts for unstranded RNA-seq
    column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
    column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
    """
    def __init__(self, filename, gtf):
        super(StarReader, self).__init__(filename, gtf)
        self.total, self.rpk_total = self.count_total()
        self.counts = self.make_counts(None)

        print(self.total)

    def count_total(self):
        """
        Total number of counts in the HTSeq table of counts.
        Also count the RPK total, for TPM calculation.
        """
        total = 0
        rpk_total = 0.0
        with open(self.filename, 'r') as my_htseq:
            for line in my_htseq:
                if '_' not in line and not line.startswith('N_'):
                    line = line.rstrip('\n').split('\t')
                    id = line[0]
                    gene_len = len(set(self.coords(id))) / 1000.0
                    count = int(line[3])
                    total += count
                    rpk_total += float(count / gene_len)
        return total, rpk_total

    def make_counts(self, count_type):
        """
        """
        counts_dict = {}
        with open(self.filename, 'r') as my_htseq:
            for line in my_htseq:
                if '_' not in line and not line.startswith('N_'):
                    line = line.rstrip('\n').split('\t')
                    id = line[0]
                    count = float(line[3])

                    if count_type == 'fpkm':
                        count = self.calc_fpkm(id, count)
                    elif count_type == 'fpkm_uq':
                        count = self.calc_fpkm_uq(id, count)
                    elif count_type == 'tpm':
                        count = self.calc_tpm(id, count)

                    counts_dict[id] = count

        return counts_dict


def main():
    args = supply_args()
    my_gtf = GTFReader(args.gtf)

    if args.use_star_counts:
        my_htseq = StarGeneCountReader(args.input, my_gtf)
    else:
        my_htseq = HtseqReader(args.input, my_gtf)

    HtseqWriter(my_htseq.counts_fpkm, args.output_fpkm)
    HtseqWriter(my_htseq.counts_fpkm_uq, args.output_fpkm_uq)
    HtseqWriter(my_htseq.counts_tpm, args.output_tpm)


if __name__ == "__main__":
    main()


