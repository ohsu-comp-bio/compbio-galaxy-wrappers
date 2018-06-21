#!/usr/bin/env python

# DESCRIPTION:
# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
import argparse

from file_types.gff3 import GffReader


VERSION = '0.1.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('star_infile', help='Input STAR SJ.out.tab file.')
    parser.add_argument('gff3', help='Input GFF3 file.')
    parser.add_argument('output', help='Output filtered STAR SJ.out.tab file.')
    parser.add_argument('--refseq_genes', help='Input list of RefSeq gene '
                                               'IDs.')
    parser.add_argument('--min_count', type=int, help='Minimum count of '
                                                      'splice event '
                                            'to include in output.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def parse_star_sj(filename):
    """
    Turn the Sj.out.tab file in to a list per line.
    :return:
    """
    star_file = []
    with filename as star:
        for line in star:
            line = line.rstrip('\n').split('\t')
            star_file.append(line)
    return star_file


def compare_norm(star_file, norm_ints):
    """
    Compare the splicesites from the GFF to what is in the Sj.out.tab file.
    :param star_file:
    :param norm_ints:
    :return:
    """
    not_ref_coords = []
    for entry in star_file:
        coord = (entry[0], int(entry[1]), int(entry[2]))
        if coord not in norm_ints:
            not_ref_coords.append(entry)
    return not_ref_coords


def add_exon_gene(not_ref, gff):
    """

    :return:
    """
    chrom, start, stop = not_ref[0], int(not_ref[1]), int(not_ref[2])
    if chrom in gff.exon_parent:
        for parent in gff.exon_parent[chrom]:
            #print(gff.coords_parent[chrom][parent])
            if (start-1 >= int(gff.coords_parent[chrom][parent][0][0])) and \
                (stop+1 <= int(gff.coords_parent[chrom][parent][0][1])):

                temp = []
                temp.append(gff.hgnc_parent[chrom][parent])
                temp.append(gff.refseq_parent[chrom][parent])
                temp.extend(['NA', 'NA'])
                for coord in sorted(gff.exon_parent[chrom][parent]):
                    if int(coord[1]) == (start - 1):
#                        not_ref.append(gff.hgnc_parent[chrom][parent])
#                        not_ref.append(gff.refseq_parent[chrom][parent])
                        temp[2] = (str(gff.exon_parent[chrom][parent].index(
                        coord)+1))
                    elif int(coord[0]) == (stop + 1):
                        temp[3] = (str(gff.exon_parent[chrom][
                                               parent].index(coord)+1))
                not_ref.extend(temp)
        return not_ref


def write_filtered_splice(args, not_ref_coords, handle):
    """
    Write the final output.
    :return:
    """
    header = ["Chromosome",
              "First Base of Intron",
              "Last Base of Intron",
              "Strand: (0: undefined, 1: +, 2: -)",
              "Intron Motif",
              "0: unannotated, 1: annotated",
              "number of uniquely mapping reads crossing the junction",
              "number of multi-mapping reads crossing the junction",
              "maximum spliced alignment overhang",
              "HGNC Symbol",
              "RefSeq ID",
              "Left Exon",
              "Right Exon"]

    handle.write('\t'.join(header))
    handle.write('\n')

    for entry in sorted(not_ref_coords, key=lambda e: int(e[6]), reverse=True):
        if int(entry[6]) > args.min_count and not entry[0].startswith('GL'):
            entry[4] = motif_map(entry[4])
            handle.write('\t'.join(entry))
            handle.write('\n')


def motif_map(motif):
    """
    Map the coded motif values to actual base values, as shown in the header.
    :return:
    """
    mapper = {'0': 'non-canonical',
              '1': 'GT/AG',
              '2': 'CT/AC',
              '3': 'GC/AG',
              '4': 'CT/GC',
              '5': 'AT/AC',
              '6': 'GT/AT'}

    return mapper[str(motif)]


def main():

    args = supply_args()
    handle_out = open(args.output, 'w')

    if args.refseq_genes:
        gff = GffReader(args.gff3, args.refseq_genes)
    else:
        gff = GffReader(args.gff3)

    norm_ints = []
    star_splice = open(args.star_infile, 'rU')
    star_file = parse_star_sj(star_splice)

    for chrom in gff.exon_parent:
        for parent in gff.exon_parent[chrom]:
            coords = sorted(gff.exon_parent[chrom][parent])
            splice_ints = [(chrom, int(coords[i][1])+1, int(coords[i + 1][0])-1) for i in range(len(coords) - 1)]
            norm_ints.extend(splice_ints)

    not_ref_coords = compare_norm(star_file, norm_ints)

    for entry in not_ref_coords:
        entry = add_exon_gene(entry, gff)

    write_filtered_splice(args, not_ref_coords, handle_out)


if __name__ == "__main__":
    main()
