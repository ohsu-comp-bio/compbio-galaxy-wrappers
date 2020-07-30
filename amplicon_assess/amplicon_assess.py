#!/usr/bin/env python

from collections import OrderedDict
from scipy import stats
import argparse
import numpy
import pysam
import vcfpy

VERSION = '0.0.1'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('bam',  help='Input BAM')
    parser.add_argument('bed',  help='Input BED listing primer regions of interest.')
    parser.add_argument('vcf',  help='Input VCF to pull primer information for.')
    parser.add_argument('--outfile', help='Output JSON')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class QiagenPrimers:
    def __init__(self, filename):
        self.filename = filename
        self.my_roi = self._parse_bed()
        self.my_primers = self._primer_regions()

    def _parse_bed(self):
        my_primers = OrderedDict()
        with open(self.filename, 'r') as myfile:
            for line in myfile:
                if not line.startswith('track'):
                    line = line.rstrip('\n').split('\t')
                    chrom = line[0][3:]
                    start = int(line[1])
                    stop = int(line[2]) - 1
                    uniq_key = (chrom, start, stop)
                    if line[5] == '+':
                        strand = 0
                    elif line[5] == '-':
                        strand = 1
                    else:
                        raise Exception("INVALID STRAND")
                    my_primers[uniq_key] = strand
        return my_primers

    @staticmethod
    def _get_pcoords(start, stop, strand):
        """
        0 is positive
        1 is negative
        If positive strand, just return the primer start position.
        If negative strand, return the primer end position.
        """
        if strand == 0:
            pstart = int(stop) - 250 + 1
            return pstart
        elif strand == 1:
            pstop = int(start) + 250 - 1
            return pstop
        else:
            raise ValueError("Invalid strand specification in get_pcoords.")

    def _primer_regions(self):
        """

        :return:
        """
        primer_coords = {}
        for coord, strand in self.my_roi.items():
            pcoord = (coord[0], self._get_pcoords(coord[1], coord[2], strand))
            primer_coords[pcoord] = strand
        return primer_coords

class BamEntry:
    """

    """
    def __init__(self, samline, primers):
        self.vrnt = samline
        self.primers = primers
        self.chrom = self.vrnt.reference_name
        self.pos = self.vrnt.reference_start
        self.name = self.vrnt.query_name
        self.direc = self._assign_direc()
        self.pstart = (self.chrom, self.pos)
        self.pend = (self.chrom, self.vrnt.reference_end - 1)
        self.unassigned = False
        self.primer = self._primer_assign()
        # self._params_stdout()

    def _params_stdout(self):
        """
        Print variable contents to stdout.
        :return:
        """

        print("Chromosome Name: {0}".format(self.chrom))
        print("Reference Position: {0}".format(self.pos))
        print("Query Name: {0}".format(self.name))
        print("Is Rev Comped: {0}".format(self.direc))
        print("Primer Start Option: {0}".format(self.pstart))
        print("Primer End Option: {0}".format(self.pend))
        print("Is Unassigned: {0}".format(self.unassigned))
        print("Primer Choice: {0}".format(self.primer))

    def _assign_direc(self):
        """
        Based on FLAG field, assign the directionality of the read.
        False = not rev comped
        True = rev comped
        :return:
        """
        if self.vrnt.flag & 0x10 == 0:
            return False
        else:
            return True

    def _check_all_primers(self, chrom, coord, dist=5):
        """
        When there is not an exact match between primer positions and the start position in the BAM, we need
        to check in the neighboring area, based on dist parameter.
        :return:
        """
        best_score = 6
        primer_choice = None
        for primer in self.primers:
            pchrom = primer[0]
            if pchrom == chrom:
                pcoord = primer[1]
                score = -dist
                for pos in range(pcoord - dist, pcoord + dist):
                    if coord == pos:
                        if abs(score) < best_score:
                            best_score = abs(score)
                            primer_choice = pcoord
                    else:
                        score += 1
        return primer_choice

    def _primer_assign(self):
        """
        Match up coords with a QiagenPrimers object.
        :param primers:
        :return:
        """
        if self.pstart in self.primers:
            return self.pstart
        elif self.pend in self.primers:
            return self.pend
        else:
            close_coord = self._check_all_primers(self.chrom, self.pstart[1])
            if close_coord:
                full_coord = (self.chrom, close_coord)
                return full_coord
            close_coord = self._check_all_primers(self.chrom, self.pend[1])
            if close_coord:
                full_coord = (self.chrom, close_coord)
                return full_coord
        self.unassigned = True


# class VrntMetrics:
#     """
#     Hold information for basecall metrics associated with each variant of interest.
#     vrnt = (chrom, pos)
#     samfile = pysam.AlignmentFile object
#     primers = QiagenPrimers.my_primers {int(POS) 0-based: int(STRAND 0 or 1)}
#     """
#     def __init__(self, samfile, vrnt, primers):
#         self.samfile = samfile
#         self.primers = primers
#         self.chrom = vrnt[0]
#         self.pos = vrnt[1]
#         self.ref = vrnt[2]
#         self.alt = [alt.value for alt in vrnt[3]][0]
#         self.raw_metrics, self.coords = self._collect_basecalls_region()
#         self.assigned = self._assign_primers()
#         # self.unassigned = self._assign_primers_rev()
#         self.total_cnt = self._get_total_cnt()
#         self.total_depth = self._get_total_depth()
#         self.total_vaf = self._get_total_vaf()
#         self.vaf = self._get_vafs()
#
#
#     def _get_total_depth(self):
#         """
#         Get the total depth for this variant.
#         :return:
#         """
#         return sum(self.total_cnt.values())
#
#     def _get_total_vaf(self):
#         """
#         Get the total VAF for this variant.
#         :return:
#         """
#         if self.alt in self.total_cnt:
#             return self.total_cnt[self.alt]/self.total_depth
#         else:
#             return None
#
#     def _get_total_cnt(self):
#         """
#         Find the total number of ref and alt counts for this variant.
#         total = {BASE: count, ...}
#         :return:
#         """
#         total = {}
#         for primer, bases in self.assigned.items():
#             for base in bases:
#                 if base not in total:
#                     total[base] = self.assigned[primer][base][0] + self.assigned[primer][base][1]
#                 else:
#                     total[base] += (self.assigned[primer][base][0] + self.assigned[primer][base][1])
#         return total
#
#     def _vaf_calc(self, counts, depth=50):
#         """
#         From this: {'T': [0, 3], 'C': [0, 17]}, calculate a VAF.
#         :return:
#         """
#         total = 0
#         for base in counts:
#             total += (counts[base][0] + counts[base][1])
#         if total < depth:
#             return None
#         if self.alt in counts:
#             vaf = (counts[self.alt][0] + counts[self.alt][1])/total
#         else:
#             vaf = 0
#         return vaf
#
#
#     def _get_vafs(self, depth=50):
#         """
#         For each primer group, estimate the VAF.  There will need to be a depth cutoff for this.
#         :return:
#         """
#         vafs = {}
#         for primer, bases in self.assigned.items():
#             vafs[primer] = self._vaf_calc(bases)
#         return vafs
#
#
#     def _assign_primers_rev(self):
#         """
#         Logic to assign counts from reverse oriented primers to the correct primer coordinate.
#         :return:
#         """
#         unassigned = []
#         # print(self.assigned)
#         # print(self.unassigned)
#         for coord in self.unassigned:
#             # print("COORD: " + str(coord))
#             ecoords = self.end_coords[coord]
#             # print("ECOORDS: " + str(ecoords))
#             if coord in self.end_coords:
#                 for ecoord in ecoords:
#                     full_ecoord = (coord[0], ecoord)
#                     # print("FULL ECOORD: " + str(full_ecoord))
#                     # print(self.primers)
#                     if ecoord in self.primers:
#                         if full_ecoord not in self.assigned:
#                             self.assigned[full_ecoord] = self.raw_metrics[coord]
#                         else:
#                             self.assigned[full_ecoord] = self._add_bases(self.assigned[full_ecoord],
#                                                                          self.raw_metrics[coord])
#                     else:
#                         close_coord = self._check_all_primers(ecoord)
#                         if close_coord:
#                             full_coord = (coord[0], close_coord)
#                             if full_coord not in self.assigned:
#                                 self.assigned[full_coord] = self.raw_metrics[coord]
#                             else:
#                                 self.assigned[full_coord] = self._add_bases(self.assigned[full_coord],
#                                                                             self.raw_metrics[coord])
#                         else:
#                             # This means we didn't assign this anywhere.
#                             unassigned.append(coord)
#             else:
#                 raise Exception("All coord tuples should be in end_coords.")
#         # print(self.assigned)
#         # print(unassigned)
#         return unassigned
#
#     def _assign_primers(self):
#         """
#         Logic to assign each count in basecalls to the appropriate primer coordinate.
#         :return:
#         """
#         primer_basecalls = {}
#         for coord, primer_choices in self.coords.items():
#             chrom = coord[0]
#             # If the coordinate is in the primer position list, assume this is the one we want.
#             # This will only match the forward-oriented primers.
#             for choice in primer_choices:
#                 pstart = (chrom, choice[0])
#                 pstop = (chrom, choice[1])
#                 idx = primer_choices.index(choice)
#                 if pstart in self.primers:
#                     if pstart not in primer_basecalls:
#                         primer_basecalls[pstart] = self.raw_metrics[coord][idx]
#                     else:
#                         primer_basecalls[pstart] = self._add_bases(primer_basecalls[pstart], self.raw_metrics[coord][idx])
#                 elif pstop in self.primers:
#                     if pstop not in primer_basecalls:
#                         primer_basecalls[pstop] = self.raw_metrics[coord][idx]
#                     else:
#                         primer_basecalls[pstop] = self._add_bases(primer_basecalls[pstop], self.raw_metrics[coord][idx])
#                 else:
#                     close_coord = self._check_all_primers(chrom, pstart[1])
#                     if close_coord:
#                         full_coord = (chrom, close_coord)
#                         if full_coord not in primer_basecalls:
#                             primer_basecalls[full_coord] = self.raw_metrics[coord][idx]
#                         else:
#                             primer_basecalls[full_coord] = self._add_bases(primer_basecalls[full_coord],
#                                                                        self.raw_metrics[coord][idx])
#                     close_coord = self._check_all_primers(chrom, pstop[1])
#                     if close_coord:
#                         full_coord = (chrom, close_coord)
#                         if full_coord not in primer_basecalls:
#                             primer_basecalls[full_coord] = self.raw_metrics[coord][idx]
#                         else:
#                             primer_basecalls[full_coord] = self._add_bases(primer_basecalls[full_coord],
#                                                                        self.raw_metrics[coord][idx])
#         print(primer_basecalls)
#         return primer_basecalls
#
#     def _check_all_primers(self, chrom, coord, dist=5):
#         """
#         When there is not an exact match between primer positions and the start position in the BAM, we need
#         to check in the neighboring area, based on dist parameter.
#         :return:
#         """
#         best_score = 6
#         primer_choice = None
#         for primer in self.primers:
#             pchrom = primer[0]
#             if pchrom == chrom:
#                 pcoord = primer[1]
#                 score = -dist
#                 for pos in range(pcoord - dist, pcoord + dist):
#                     if coord == pos:
#                         if abs(score) < best_score:
#                             best_score = abs(score)
#                             primer_choice = pcoord
#                     else:
#                         score += 1
#         return primer_choice
#
#     @staticmethod
#     def _add_bases(old, new):
#         """
#         This will handle combining dicts of {base: [fwdcount, revcount]} pairs.
#         OLD: {'G': (1, 0)}
#         NEW: {'G': (0, 1)} =
#             {'G': (1, 1)}
#         :return:
#         """
#         new_val = old
#         for base in new:
#             if base not in new_val:
#                 new_val[base] = new[base]
#             else:
#                 new_val[base][0] += new[base][0]
#                 new_val[base][1] += new[base][1]
#         return new_val
#
#     def _assign_direc(self, flag):
#         """
#         Based on FLAG field, assign the directionality of the read.
#         :return:
#         """
#         if flag & 0x10 == 0:
#             return 0
#         else:
#             return 1
#
#     def _del_from_cigar(self, cigar, pos):
#         """
#         Get the value of the deletion operation, if it exists.
#         [(0, 54), (2, 1), (0, 24)]
#         :return:
#         """
#         total_del = 0
#         cnt_opt = ['0', '3', '7', '8']
#         if 2 in [x[0] for x in cigar]:
#             scigar = ''.join([str(x[0]) * x[1] for x in cigar])
#             for op in scigar:
#                 if op in cnt_opt:
#                     cnt += 1
#
#             del_loc = scigar.count('2')
#             return del_loc
#         return total_del
#
#     def _collect_basecalls_region(self):
#         """
#         Get the base counts for each read in the BAM, based on given start position.
#         :return:
#         """
#         mysam = self.samfile.fetch(self.chrom, self.pos, self.pos+1)
#         basecalls = {}
#         coords = {}
#         for line in mysam:
#             # This ain't gonna work without a correction for deletions.
#             del_correct = self._del_from_cigar(line.cigartuples, self.pos-line.pos)
#             srch_idx = self.pos - line.pos - del_correct
#             basecall = line.query_sequence[srch_idx]
#             uniq_key = (line.reference_name, line.pos)
#             direc = self._assign_direc(line.flag)
#             primer_pair = (line.reference_start, line.reference_end-1)
#
#             # Fill the coords dict.
#             if uniq_key not in coords:
#                 coords[uniq_key] = [primer_pair]
#             else:
#                 coords[uniq_key].append(primer_pair)
#
#             # Fill the basecalls dict.
#             if uniq_key not in basecalls:
#                 basecalls[uniq_key] = []
#             if direc == 0:
#                 basecalls[uniq_key].append({basecall: [1, 0]})
#             else:
#                 basecalls[uniq_key].append({basecall: [0, 1]})
#         return basecalls, coords

class Vrnt:
    """

    """
    def __init__(self, vrnt, samfile):
        self.vrnt = vrnt
        self.samfile = samfile
        self.chrom = self.vrnt.CHROM
        self.pos = self.vrnt.POS - 1
        self.ref = self.vrnt.REF
        self.alt = [alt.value for alt in self.vrnt.ALT]
        if len(self.alt) == 1:
            self.alt = self.alt[0]
            self.type = self.vrnt.ALT[0].type
            self.altlen = self._set_vrnt_len()
        else:
            raise ValueError("Only one alternate allele from your VCF record is currently supported.")
        self.uniq_key = (self.chrom, self.pos, self.ref, self.alt)
        self.pileup = self._get_pileup()
        # self._params_stdout()

    def _set_vrnt_len(self):
        """

        :return:
        """
        if self.type == 'INS':
            return len(self.alt) - 1
        if self.type == 'DEL':
            return len(self.ref) - 1

    def _get_pileup(self):
        """
        Figure out how many counts there are in the BAM, at this location, for each base present.
        :return:
        """
        pile = {}
        start = self.pos
        stop = self.pos + 1

        for pileupColumn in self.samfile.pileup(self.chrom, start, stop, stepper='samtools',
                                                min_base_quality=0, ignore_overlaps=False, ignore_orphans=False,
                                                truncate=True):
            for pileupRead in pileupColumn.pileups:
                name = pileupRead.alignment.query_name
                base = pileupRead.alignment.query_sequence[pileupRead.query_position]
                if pileupRead.indel:
                    if abs(pileupRead.indel) == self.altlen:
                        base = self.type
                    else:
                        label = 'BAD_' + str(self.type)
                        base = label
                if name not in pile:
                    pile[name] = base
                else:
                    raise Exception("NOT POSSIBLE")
        return pile

    def _params_stdout(self):
        """
        Print variable contents to stdout.
        :return:
        """
        print("Variant Chromosome: {0}".format(self.chrom))
        print("Variant Position: {0}".format(self.pos))
        print("Variant Ref Allele: {0}".format(self.ref))
        print("Variant Alt Allele: {0}".format(self.alt))
        print("Variant Type: {0}".format(self.type))
        print("BAM Base Counts: {0}".format(self.pileup))


class VrntMetrics:
    """
    """
    def __init__(self, vrnts, samfile, primers):
        self.vrnts = vrnts
        self.primers = primers
        self.samfile = samfile

        self._gather()
        # self.raw_metrics, self.coords = self._collect_basecalls_region()
        # self.assigned = self._assign_primers()
        # # self.unassigned = self._assign_primers_rev()
        # self.total_cnt = self._get_total_cnt()
        # self.total_depth = self._get_total_depth()
        # self.total_vaf = self._get_total_vaf()
        # self.vaf = self._get_vafs()


    def _gather(self):
        """
        vrnt.pileup = {qname1: base, q   name2: base, etc.}
        """
        primer_counts = {}
        for vrnt in self.vrnts:
            primer_counts[vrnt.uniq_key] = {}
            mybam = self.samfile.fetch(vrnt.chrom, vrnt.pos, vrnt.pos + 1)
            for line in mybam:
                seq = BamEntry(line, self.primers)
                if seq.name in vrnt.pileup:
                    base = vrnt.pileup[seq.name]
                    if seq.primer not in primer_counts[vrnt.uniq_key]:
                        primer_counts[vrnt.uniq_key][seq.primer] = {}
                    if base not in primer_counts[vrnt.uniq_key][seq.primer]:
                        primer_counts[vrnt.uniq_key][seq.primer][base] = 1
                    else:
                        primer_counts[vrnt.uniq_key][seq.primer][base] += 1

        print(primer_counts)
            #     if seq.primer not in primer_counts[entry.uniq_key]:
            #         primer_counts[entry.uniq_key][seq.primer] = [line]
            #     else:
            #         primer_counts[entry.uniq_key][seq.primer].append(line)
            # for primer in primer_counts[entry.uniq_key].values():
            #     for read in primer:
            #         print(read.fetch(entry.chrom, entry.pos, entry.pos+1))
        return primer_counts


def main():
    """
    AMPBIAS metrics options:
    1. Use a z-score or something along those lines.  This may not work very well since we typically
    will not have more than 3 or 4 overlapping measurements of VAF.
    2. Note when one of the amplicons doesn't have any variant counts, this is certainly not good.
    Thresholding should include some minimum depth.
    :return:
    """
    args = supply_args()
    samfile = pysam.AlignmentFile(args.bam, "rb")
    primers = QiagenPrimers(args.bed).my_primers
    reader = vcfpy.Reader.from_path(args.vcf)
    vrnts = []
    for entry in reader:
        vrnts.append(Vrnt(entry, samfile))

    VrntMetrics(vrnts, samfile, primers)

    #
    # for entry in reader:
        # vrnt = (entry.CHROM, entry.POS-1, entry.REF, entry.ALT)
        # If this is too slow, we can prepare primers list ahead to be localized to region, or even chrom.
        # all_mets.append(VrntMetrics(samfile, vrnt, primers))


    samfile.close()
    reader.close()

    # for entry in all_mets:
    #     temp = []
    #     for line in entry.vaf.values():
    #         if line:
    #             temp.append(line)
    #         if len(temp) > 2:
    #             a = numpy.array(temp)
    #             print(stats.zscore(a))
    #             print(entry.total_vaf)
    #             print(entry.vaf)
    #             print(entry.total_cnt)
    #             print(entry.assigned)
    #             print('\n')




if __name__ == "__main__":
    main()
