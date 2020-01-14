#!/usr/bin/env python

# DESCRIPTION: Clip primers based on a BED file.
# USAGE: primer_clip.py --bed <BED> --bamfile <SAM> --outfile <OUTPUT SAM>
# CODED BY: John Letaw

from __future__ import print_function
from file_types import bed
import argparse
import pysam
import re

VERSION = '0.1.4'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--bed', help='Input BED file.')
    parser.add_argument('--bamfile', help='Input SAM file.')
    parser.add_argument('--outfile', help='Output SAM file.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def get_pcoords(start, stop, strand):
    """
    0 is positive
    1 is negative
    """
    if strand == 0:
        pstart = stop - 250 + 1
        pstop = start - 1
        return pstop
    else:
        pstart = stop + 1
        pstop = start + 250 - 1
        return pstart
#    return range(pstart, pstop)


def soft_clip_offset(entry):
    """
    Tell me if there are soft-clipped bases in the CIGAR string.
    :return:
    """
    if entry[0][0] == 4:
        return entry[0][1]
    return 0

def cigar_adj(change, cigar):
    """
    Op BAM Description ConsumesQuery ConsumesReference
    M 0 alignment match (can be a sequence match or mismatch) yes yes
    I 1 insertion to the reference yes no
    D 2 deletion from the reference no yes
    N 3 skipped region from the reference no yes
    S 4 soft clipping (clipped sequences present in SEQ) yes no
    H 5 hard clipping (clipped sequences NOT present in SEQ) no no
    P 6 padding (silent deletion from padded reference) no no
    = 7 sequence match yes yes
    X 8 sequence mismatch yes yes
    :return:
    """
    soft_clip = 0
    new_cigar = []
    for entry in cigar:
        field = entry[0]
        value = entry[1]
        if change != 0:
        # yes/yes
            if field == 0 or field == 7 or field == 8:
                if value >= change:
                    value -= change
                    soft_clip += change
                    change = 0
                    new_cigar.append([field, value])
                else:
                    change -= value
                    soft_clip += value
        # yes/no insertions
            elif field == 1:
                soft_clip += value
        # no/yes deletions
            elif field == 2 or field == 3:
                change -= value
            # soft clipping
            elif field == 4:
                new_cigar.append([field, value])
            # no/no (5 or 6)
            else:
                new_cigar.append([field, value])
        else:
            new_cigar.append([field, value])

    try:
        if new_cigar[0][0] == 5:
            if new_cigar[1][0] == 4:
                new_cigar[1][1] += soft_clip
            else:
                new_cigar.insert(1, (4, soft_clip))
        else:
            new_cigar.insert(0, (4, soft_clip))
    except:
        pass

    return new_cigar


def change_cigar_format(cigar):
    """
    Transform form 6S144M format to [(4, 6), (0, 144)]
    :return:
    """
    mapping = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5}
    r = re.compile(r"([0-9]+)([MIDNSHPX])")
    new_cigar = []
    for entry in re.findall(r, cigar):
        new_cigar.append((mapping[entry[1]], int(entry[0])))

    return new_cigar


def flip_cigar(cigar):
    """

    :param cigar:
    :return:
    """
    r = re.compile(r"([0-9]+[MIDNSHPX])")
    new_cigar = []
    for entry in re.findall(r, cigar):
        new_cigar.insert(0, entry)

    return ''.join(new_cigar)


def cigar_for_writing(cigar):
    """
    Transform change_cigar_format() output back to writable 6S144M format.
    :return:
    """
    mapping = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H'}
    new_cigar = ''
    for entry in cigar:
        to_add = str(entry[1]) + mapping[entry[0]]
        new_cigar += to_add

    return new_cigar


def check_len(cigar, minlen=10):
    """
    Check to make sure we aren't dealing with an all hard-coded sequence,
    or even one that is just too short to be useful.  Will remove both ends
    of the pair for the time being.
    :return:
    """
    for entry in cigar:
        field = entry[0]
        value = entry[1]
        if value <= 0:
            return False
        if field == 0:
            if value >= minlen:
                return True
    return False


def strip_hc(entry, strip, rev=False):
    """
    Take a BAM entry, remove the hardclipped sequence from the entry.seq
    and entry.qual.
    :param entry:
    :param strip:
    :return:
    """
    # new_seq = str(entry.qual)
    # entry.seq = entry.seq[strip:]
    # entry.query_qualities = pysam.qualitystring_to_array(new_seq[strip:])

    if rev:
        entry[9] = entry[9][:-strip]
        entry[10] = entry[10][:-strip]
    else:
        entry[9] = entry[9][strip:]
        entry[10] = entry[10][strip:]
        entry[3] = str(int(entry[3]) + strip)
    return entry


def val_r1(line1, primer_coords, samfile):
    """

    :param entry:
    :param primer_coords:
    :param samfile:
    :return:
    """
    flag = int(line1[1])
    chrom = str(line1[2])
    while ((flag & 0x800 != 0) or
               (flag & 0x8 != 0) or
               (flag & 0x2 == 0) or
               (flag & 0x40 == 0) or
               (chrom not in primer_coords)):
        try:
            line1 = next(samfile).rstrip('\n').split('\t')
        except StopIteration:
            break
        flag = int(line1[1])
        chrom = str(line1[2])

    return line1, samfile


def val_r2(line2, primer_coords, samfile):
    """

    :param entry:
    :param primer_coords:
    :param samfile:
    :return:
    """
    flag2 = int(line2[1])
    chrom = str(line2[2])
    while ((flag2 & 0x800 != 0) or
           (flag2 & 0x80 == 0) or
           (flag2 & 0x2 == 0) or
           (chrom not in primer_coords)):
        try:
            line2 = next(samfile).rstrip('\n').split('\t')
        except StopIteration:
            break
        flag2 = int(line2[1])
        chrom = str(line2[2])

    return line2, samfile


def fix_cigar(cigar):
    """
    This code is getting ugly.  If the CIGAR string contains two of the same type of field after another, merge.
    69M20S39S should be 69M59S
    :return:
    """
    new_cigar = []
    for i in range(len(cigar)):
        entry = cigar[i]
        field = entry[0]
        value = entry[1]
        if i == 0:
            new_cigar.append(list(entry))
        else:
            if field == new_cigar[-1][0]:
                new_cigar[-1][1] += value
            elif value != 0:
                new_cigar.append(list(entry))

    new_cigar = [tuple(x) for x in new_cigar]
    return new_cigar


def main():
    args = supply_args()
    samfile = open(args.bamfile, 'r')
    outfile = open(args.outfile, 'w')
    primer_coords = bed.ExtBedReader(args.bed, header=True, strand=5).get_primer_ends()
    mismatch = False

    i = 0
    with samfile as sam:
        for entry in sam:
            if i % 1000000 == 0:
                    print(str(i) + " read pairs processed.")
            to_remove = 0
            if entry.startswith('@'):
                outfile.write(entry)
            else:
                if not mismatch:
                    line1 = entry.rstrip('\n').split('\t')
                    line1, samfile = val_r1(line1, primer_coords, samfile)
                    id1 = line1[0]

                    try:
                        entry = next(samfile)
                    except StopIteration:
                        break

                    line2 = entry.rstrip('\n').split('\t')
                    line2, samfile = val_r2(line2, primer_coords, samfile)
                    id2 = line2[0]

                    if id1 != id2:
                        mismatch = True
                        line1 = line2
                        continue
                else:
                    line1, samfile = val_r1(line1, primer_coords, samfile)
                    line2 = entry.rstrip('\n').split('\t')
                    line2, samfile = val_r2(line2, primer_coords, samfile)

                    mismatch = False

                if line1[0] != line2[0]:
                    print(line1)
                    print(line2)
                    print("Reads 1 and 2 don't match, ouch!")
                    continue

                pos = int(line1[3])
                cigar = change_cigar_format(line1[5])
                pnext = int(line1[7])
                tlen = int(line1[8])
                sco = soft_clip_offset(cigar)

                l2_pos = int(line2[3])
                l2_cigar = change_cigar_format(line2[5])

                flag = int(line1[1])
                chrom = str(line1[2])

                # If first read is mapped to fwd strand.
                if not (flag & 0x10):
                    seqloc = (chrom, range(pos + sco, pos + tlen))
                    primer_choices = set(primer_coords[seqloc[0]][
                                '1']).intersection(set(seqloc[1]))

                    if primer_choices:
                        for coord in primer_choices:
                            check_remove = coord - pos
                            if check_remove < 50 and check_remove >= 5:
                                to_remove = check_remove
                        cigar = cigar_adj(to_remove, cigar)
                        line2_cigar_adj = to_remove - (l2_pos - pos)
                        line1[3] = str(int(line1[3]) + to_remove)
                        if line2_cigar_adj <= to_remove and line2_cigar_adj > 0:
                            l2_cigar = cigar_adj(line2_cigar_adj, l2_cigar)
                            line2[3] = str(int(line2[3]) + line2_cigar_adj)
                            line1[7] = line2[3]

                        line2[7] = line1[3]

                        line1[5] = cigar_for_writing(fix_cigar(cigar))
                        line2[5] = cigar_for_writing(fix_cigar(l2_cigar))
                        if check_len(cigar) and check_len(l2_cigar):
                            outfile.write('\t'.join(line1))
                            outfile.write('\n')
                            outfile.write('\t'.join(line2))
                            outfile.write('\n')
                else:
                    seqloc = (chrom, range(pos, (pnext - tlen) - sco))
                    primer_choices = set(primer_coords[seqloc[0]][
                                '0']).intersection(set(seqloc[1]))
                    if primer_choices:
                         for coord in primer_choices:
                             check_remove = (pnext - tlen) - coord
                             if check_remove < 50 and check_remove >= 5:
                                 to_remove = check_remove
                             else:
                                 to_remove = 0
                         if to_remove != 0:
                             cigar = cigar_adj(to_remove, cigar[::-1])
                             line2_cigar_adj = to_remove - (pos - l2_pos)
                             if line2_cigar_adj <= to_remove and line2_cigar_adj > 0:
                                 l2_cigar = cigar_adj(line2_cigar_adj, l2_cigar[::-1])
                                 line2[5] = cigar_for_writing(fix_cigar(l2_cigar[::-1]))
                             else:
                                 line2[5] = cigar_for_writing(fix_cigar(l2_cigar))
                             line1[5] = cigar_for_writing(fix_cigar(cigar[::-1]))
                             if check_len(cigar) and check_len(l2_cigar):
                                 outfile.write('\t'.join(line1))
                                 outfile.write('\n')
                                 outfile.write('\t'.join(line2))
                                 outfile.write('\n')

            i += 1

    outfile.close()


def main2():
    """
    This is the method that uses pysam.  I don't really like this method,
    and it is slower, but keeping this code for the hell of it.
    :return:
    """
    args = supply_args()
    samfile = pysam.AlignmentFile(args.bamfile, 'rb')
    outfile = pysam.AlignmentFile(args.outfile, 'wb', template=samfile)
    primer_coords = bed.ExtBedReader(args.bed, header=True, strand=5).get_primer_ends()

    i = 0
    for entry in samfile:

        if i % 100000 == 0:
            print(str(i) + " read pairs processed.")

        line1 = entry
        try:
            line2 = next(samfile)
        except:
            raise Exception("WHY AM I HERE?")

        chrom = samfile.header['SQ'][line1.rname]['SN']

        if (line1.flag & 0x800 == 0) and \
            (line1.flag & 0x8 == 0) and \
            (line1.flag & 0x2 != 0) and \
                (chrom in primer_coords):
            sco = soft_clip_offset(line1.cigar)
            if (line1.flag & 0x10 == 0):
                seqloc = (chrom, range(line1.pos + sco, line1.pos +
                                        line1.tlen))
                z = set(primer_coords[seqloc[0]]['1']).intersection(
                    set(seqloc[1]))
                if z:
                    for coord in list(z):
                        check_remove = coord-line1.pos
                        if check_remove < 50 and check_remove >= 0:
                            to_remove = check_remove
                    line1.cigar, line1_hc = cigar_adj(to_remove, line1.cigar)
                    line2_cigar_adj = to_remove - (line2.pos - line1.pos)
                    if line2_cigar_adj <= to_remove and line2_cigar_adj > 0:
                        line2.cigar, line2_hc = cigar_adj(line2_cigar_adj,
                                                   line2.cigar)
                        line2 = strip_hc(line2, line2_hc)

                    line1 = strip_hc(line1, line1_hc)
                    if check_len(line1.cigar) and check_len(line2.cigar):
                        outfile.write(line1)
                        outfile.write(line2)

            else:
                seqloc = (chrom, range(line1.pos, (line1.pnext - line1.tlen) - sco))
                z = set(primer_coords[seqloc[0]]['0']).intersection(
                    set(seqloc[1]))
                if z:
                    for coord in list(z):
                        check_remove = (line1.pnext-line1.tlen)-coord
                        if check_remove < 50 and check_remove >= 0:
                            to_remove = check_remove
                    line1.cigar, line1_hc = cigar_adj(to_remove, line1.cigar)
                    line2_cigar_adj = to_remove - (line1.pos - line2.pos)
                    if line2_cigar_adj <= to_remove and line2_cigar_adj > 0:
                        line2.cigar, line2_hc = cigar_adj(line2_cigar_adj,
                                                   line2.cigar)
                        line2 = strip_hc(line2, line2_hc)

                    line1 = strip_hc(line1, line1_hc)
                    if check_len(line1.cigar) and check_len(line2.cigar):
                        outfile.write(line1)
                        outfile.write(line2)
        i += 1

    outfile.close()

if __name__ == "__main__":
    main()


