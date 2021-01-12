#!/usr/bin/env python

import argparse

VERSION = '0.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input', help='Input QIAseq primer list (BED format).')
    parser.add_argument('genes', help='Input genes list.')
    parser.add_argument('output', help='Output TSV.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

def main():
    args = supply_args()
    writer = open(args.output, 'w')

    genes = []
    with open(args.genes, 'r') as mygenes:
        for gene in mygenes:
            gene = gene.rstrip('\n')
            if gene not in genes:
                genes.append(gene)

    # The STPv3_225 portion is hard-coded in to the Rscript that works on this file,
    # so that needs to be changed before I can use a reasonable header here.
    header = ['CONTIG', 'START', 'END', 'STRINGhg19', 'GENE', 'STPv3_225']
    writer.write('\t'.join(header))
    writer.write('\n')

    output = []
    with open(args.input, 'r') as myfile:
        for line in myfile:
            line = line.rstrip('\n').split('\t')
            gene = line[3].split(';')
            coord = 'chr' + line[0] + ':' + line[1] + '-' + line[2]
            to_write = [line[0], line[1], line[2], coord]

            for entry in gene:
                if entry in genes:
                    gene = [entry]

            to_write.extend(gene)

            if gene[0] in genes:
                to_write.append('1')
            else:
                to_write.append('0')

            output.append(to_write)

    coords = {}
    for entry in output:
        chrom = entry[0]
        start = int(entry[1])
        stop = int(entry[2])
        if chrom not in coords:
            coords[chrom] = []

        for pos in range(start, stop):
            if pos not in coords[chrom]:
                coords[chrom].append(pos)

    for chrom, pos in coords.items():
        for p in sorted(pos):
            pass


    # writer.write('\t'.join(to_write))
    # writer.write('\n')

    writer.close()

if __name__ == "__main__":
    main()


