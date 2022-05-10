#!/usr/bin/env python

# Sort BEDs based on contig list below.
# Only output those chromosomes contained in the following list.
# python sort_bed.py <in_bed> <out_bed>
# CODED BY: John Letaw

import sys

CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT"


def create_bed_dict(handle):
    """
    Create BED dictionary of form:
    chrom: [start, stop]
    """
    bed_dict = {}
    with handle as bed:
        for coords in bed:
            chrom = coords.rstrip('\n').split('\t')[0]
            if chrom[:3] == 'chr':
                chrom = chrom[3:]
            start = coords.rstrip('\n').split('\t')[1]
            stop = coords.rstrip('\n').split('\t')[2]

            if chrom not in bed_dict:
                bed_dict[chrom] = []

            bed_dict[chrom].append([start, stop])

    return bed_dict


def write_bed(handle, bed):
    """
    Write from a chrom: [start, stop] dict to file.
    """
    for chrom in CHROMS.split(' '):
        if chrom in bed:
            for coord in sorted(bed[chrom], key=lambda x: int(x[0])):
                to_write = [chrom, coord[0], coord[1], '\n']
                handle.write('\t'.join(to_write))


def main():

    handle = open(sys.argv[1], 'r')
    handle_out = open(sys.argv[2], 'w')
    bed_dict = create_bed_dict(handle)
    write_bed(handle_out, bed_dict)
    handle_out.close()


if __name__ == "__main__":
    main()
