#!/usr/bin/env python

# DESCRIPTION:
# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
import argparse

VERSION = '0.1.5'

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


class GffReader(object):
    """
    Reads a GFF3 file and creates several data structures relating portions
    of the GFF info field to the parent id of the particular feature.

    Currently this takes a list of RefSeq transcript id's as well, to limit
    the search space for BED coordinates.  This will change when I have a
    chance to put more effort in to this.

    Input: filename, refseq_file

    * Describe seqnames
    region

    * No parent
    gene

    *Parent is "gene"
    primary_transcript
    transcript
    mRNA
    ncRNA
    rRNA
    tRNA

    *Parent is in above group.
    CDS
    exon

    *Singletons
    cDNA_match
    match
    C_gene_segment
    D_gene_segment
    J_gene_segment
    D_loop
    V_gene_segment

    """
    def __init__(self, filename, refseq_file=None, chrom_list=None):

        """
        Set following data structure:

        refseq - List of RefSeq id's, from file
        header - List of lines containing # as line.startswith()
        refseq_parent - Relates RefSeq NM values to the feature parent id
        hgnc_parent - Relates HGNC symbols to the feature parent id
        cds_parent - Relates CDS coordinates to the feature parent id
        exon_parent - Relates exon coordinates to the feature parent id
        coords_parent - Relates mRNA transcript coordinates to the feature
        parent id

        :param filename:
        :param refseq_file:
        """

        self.filename = filename
        self.refseq_file = refseq_file
        self.chrom_list = chrom_list

        if self.refseq_file is not None:
            self.refseq = self._parse_refseq(self.refseq_file)
            self.refseq_parent = self._refseq_parent(self.refseq,
                                                     feature='mRNA',
                                                     ident='Name')
        else:
            self.refseq_parent = self._refseq_parent()

        self.header = self._get_header()
        self.hgnc_parent = self._info_parent(self.refseq_parent,
                                               feature='mRNA',
                                               ident='gene')
        self.cds_parent = self._cds_parent(self.refseq_parent, feature='CDS')
        self.exon_parent = self._cds_parent(self.refseq_parent, feature='exon')
        self.coords_parent = self._cds_parent(self.refseq_parent,
                                              feature='mRNA', ident='ID')

    def _info_parent(self, refseq_parent, feature='mRNA', ident='Name'):

        """
        Connect stuff in info field to Parent ID's based on feature and the
        field in the info entry.

        :param refseq_parent:
        :param feature:
        :param ident:
        :return info_parent:
        """

        info_parent = {}
        with open(self.filename, 'rU') as gff3:
            for this_line in gff3:
                if not this_line.startswith('#'):
                    this_line = this_line.rstrip('\n').split('\t')
                    chrom = self._refseq_to_common_chrom(this_line[0])
                    if chrom is not None:
                        info = this_line[8]
                        this_feature = this_line[2]
                        if chrom not in info_parent:
                            info_parent[chrom] = {}
                        if this_feature == feature:
                            parent = self._find_info_fields(info, 'ID')
                            this_id = self._find_info_fields(info, ident)
                            if parent in refseq_parent[chrom]:
                                if parent not in info_parent[chrom]:
                                    info_parent[chrom][parent] = this_id

        return info_parent

    def _parse_refseq(self, filename):

        """
        Put RefSeq ID's from file in to a list.
        :param filename:
        :return refseq_list:
        """

        refseq_list = []
        with open(filename, 'rU') as refseq:
            for line in refseq:
                line = line.rstrip('\n')
                refseq_list.append(line)

        return refseq_list

    def _cds_parent(self, refseq_parent, feature='CDS', ident='Parent'):

        """
        Get 'feature' sequences, and put them in a dictionary with 'ident'.
        This is in position [8] and follows "Parent=(rna[0-9]+);"

        :param refseq_parent:
        :param feature:
        :param ident:
        :return cds_parent:
        """

        cds_parent = {}
        with open(self.filename, 'rU') as gff3:
            for this_line in gff3:
                if not this_line.startswith('#'):
                    this_line = this_line.rstrip('\n').split('\t')
                    chrom = self._refseq_to_common_chrom(this_line[0])
                    if chrom == None:
                        continue
                    info = this_line[8]
                    this_feature = this_line[2]
                    if chrom not in cds_parent:
                        cds_parent[chrom] = {}
                    if this_feature == feature:
                        parent = self._find_info_fields(info, ident)
                        start = this_line[3]
                        stop = this_line[4]
                        if parent in refseq_parent[chrom]:
                            if parent not in cds_parent[chrom]:
                                cds_parent[chrom][parent] = []
                            cds_parent[chrom][parent].append([start, stop])

        return cds_parent


    def _refseq_parent(self, refseq=None, feature='mRNA', ident='Name'):

        """
        Connect 'feature' id's to 'ident' id's.  Differs from _cds_parent in
        that we are not tracking coordinates here.
        NOTE: Come back to this and refactor based on overlap with _cds_parent.

        :param refseq:
        :param feature:
        :param ident:
        :return refseq_parent:
        """

        refseq_parent = {}
        with open(self.filename, 'rU') as gff3:
            for this_line in gff3:
                if not this_line.startswith('#'):
                    this_line = this_line.rstrip('\n').split('\t')
                    chrom = self._refseq_to_common_chrom(this_line[0])
                    if chrom is not None:
                        info = this_line[8]
                        this_feature = this_line[2]
                        if chrom not in refseq_parent:
                            refseq_parent[chrom] = {}
                        if this_feature == feature:
                            parent = self._find_info_fields(info, 'ID')
                            refseq_id = self._find_info_fields(info, ident)
                            if refseq is not None:
                                if refseq_id in refseq:
                                    if parent not in refseq_parent[chrom]:
                                        refseq_parent[chrom][parent] = refseq_id
                            else:
                                if parent not in refseq_parent[chrom]:
                                    refseq_parent[chrom][parent] = refseq_id

        return refseq_parent


    def _get_header(self):

        """
        Isolate the header section of the GFF3.  This is a list of entries
        beginning with the '#' symbol.
        :return header:
        """

        header = []
        with open(self.filename, 'rU') as gff3:
            for this_line in gff3:
                if this_line.startswith('#'):
                    header.append(this_line.rstrip('\n'))
        return header


    def _find_info_fields(self, info, field):

        """
        Find a field from the INFO column.  Return the value based on a
        particular field title.

        :return entry[1]:
        """

        for entry in info.split(';'):
            entry = entry.split('=')
            if entry[0] == field:
                return entry[1]

    def _refseq_to_common_chrom(self, chrom):

        """
        Set dictionary of RefSeq to chromosome identifiers.  Hard code this
        for the time being until I figure out a better way to deal with this
        situation.
        :return refseq_to_chrom[chrom] if in refseq_to_chrom:
        :else return None:
        """

        refseq_to_chrom = {'NC_012920.1': 'MT', 'NC_000001.10': '1',
                           'NC_000002.11': '2', 'NC_000003.11': '3',
                           'NC_000004.11': '4', 'NC_000005.9': '5',
                           'NC_000006.11': '6', 'NC_000007.13': '7',
                           'NC_000008.10': '8', 'NC_000009.11': '9',
                           'NC_000010.10': '10', 'NC_000011.9': '11',
                           'NC_000012.11': '12', 'NC_000013.10': '13',
                           'NC_000014.8': '14', 'NC_000015.9': '15',
                           'NC_000016.9': '16', 'NC_000017.10': '17',
                           'NC_000018.9': '18', 'NC_000019.9': '19',
                           'NC_000020.10': '20', 'NC_000021.8': '21',
                           'NC_000022.10': '22', 'NC_000023.10': 'X',
                           'NC_000024.9': 'Y'}

        try:
            if self.chrom_list is not None:
                if refseq_to_chrom[chrom] in self.chrom_list:
                    return refseq_to_chrom[chrom]
            else:
                return refseq_to_chrom[chrom]
        except:
            return None


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
        if entry[0].startswith('chr'):
            entry[0] = entry[0][3:]
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
            if (start-1 >= int(gff.coords_parent[chrom][parent][0][0])) and \
                (stop+1 <= int(gff.coords_parent[chrom][parent][0][1])):

                temp = []
                temp.append(gff.hgnc_parent[chrom][parent])
                temp.append(gff.refseq_parent[chrom][parent])
                temp.extend(['NA', 'NA'])
                for coord in sorted(gff.exon_parent[chrom][parent]):
                    if int(coord[1]) == (start - 1):
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
        # filtering MT here as well
        if int(entry[6]) > args.min_count and not entry[0].startswith('GL') and not entry[0].startswith('MT') and not entry[0] == 'M':
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
