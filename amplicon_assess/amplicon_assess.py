#!/usr/bin/env python

from _collections import OrderedDict
import argparse
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
        else:
            pstop = int(start) + 250 - 1
            return pstop

    def _primer_regions(self):
        """

        :return:
        """
        primer_coords = {}
        for coord, strand in self.my_roi.items():
            pcoord = (coord[0], self._get_pcoords(coord[1], coord[2], strand))
            primer_coords[pcoord] = strand
        return primer_coords


class VrntMetrics:
    """
    Hold information for basecall metrics associated with each variant of interest.
    vrnt = (chrom, pos)
    samfile = pysam.AlignmentFile object
    primers = QiagenPrimers.my_primers {int(POS) 0-based: int(STRAND 0 or 1)}
    """
    def __init__(self, samfile, vrnt, primers):
        self.samfile = samfile
        self.primers = primers
        self.chrom = vrnt[0]
        self.pos = vrnt[1]
        self.raw_metrics, self.end_coords = self._collect_basecalls_region()
        self.assigned, self.unassigned = self._assign_primers()
        self.unassigned = self._assign_primers_rev()

    def _assign_primers_rev(self):
        """
        Logic to assign counts from reverse oriented primers to the correct primer coordinate.
        :return:
        """
        unassigned = []
        for coord in self.unassigned:
            ecoords = self.end_coords[coord]
            if coord in self.end_coords:
                for ecoord in ecoords:
                    full_ecoord = (coord[0], ecoord)
                    if ecoord in self.primers:
                        if full_ecoord not in self.assigned:
                            self.assigned[full_ecoord] = self.raw_metrics[coord]
                        else:
                            self.assigned[full_ecoord] = self._add_bases(self.assigned[full_ecoord],
                                                                         self.raw_metrics[coord])
                    else:
                        close_coord = self._check_all_primers(ecoord)
                        if close_coord:
                            full_coord = (coord[0], close_coord)
                            if full_coord not in self.assigned:
                                self.assigned[full_coord] = self.raw_metrics[coord]
                            else:
                                self.assigned[full_coord] = self._add_bases(self.assigned[full_coord],
                                                                            self.raw_metrics[coord])
                        else:
                            # This means we didn't assign this anywhere.
                            unassigned.append(coord)
            else:
                raise Exception("All coord tuples should be in end_coords.")

        return unassigned

    def _assign_primers(self):
        """
        Logic to assign each count in basecalls to the appropriate primer coordinate.
        :return:
        """
        primer_basecalls = {}
        unassigned = []
        for coord in self.raw_metrics:
            # If the coordinate is in the primer position list, assume this is the one we want.
            # This will only match the forward-oriented primers.
            if coord in self.primers:
                if coord not in primer_basecalls:
                    primer_basecalls[coord] = self.raw_metrics[coord]
                else:
                    primer_basecalls[coord] = self._add_bases(primer_basecalls[coord], self.raw_metrics[coord])
            # Now check within a range.
            else:
                close_coord = self._check_all_primers(coord[1])
                full_coord = (coord[0], close_coord)
                if close_coord:
                    if full_coord not in primer_basecalls:
                        primer_basecalls[full_coord] = self.raw_metrics[coord]
                    else:
                        primer_basecalls[full_coord] = self._add_bases(primer_basecalls[full_coord],
                                                                       self.raw_metrics[coord])
                else:
                    # This means we didn't assign this anywhere.
                    unassigned.append(coord)

        return primer_basecalls, unassigned

    def _check_all_primers(self, coord, dist=5):
        """
        When there is not an exact match between primer positions and the start position in the BAM, we need
        to check in the neighboring area, based on dist parameter.
        :return:
        """
        best_score = 6
        primer_choice = None
        for primer in self.primers:
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

    @staticmethod
    def _add_bases(old, new):
        """
        This will handle combining dicts of {base: count} pairs.
        :return:
        """
        new_val = old
        for base in new:
            if base not in new_val:
                new_val[base] = new[base]
            else:
                new_val[base] += new[base]
        return new_val

    def _collect_basecalls_region(self):
        """
        Get the base counts for each read in the BAM, based on given start position.
        :return:
        """
        mysam = self.samfile.fetch(self.chrom, self.pos, self.pos+1)
        basecalls = {}
        end_coords = {}
        for line in mysam:
            srch_idx = self.pos - line.pos
            basecall = line.query_sequence[srch_idx]
            uniq_key = (line.rname, line.pos)

            if uniq_key not in end_coords:
                end_coords[uniq_key] = [line.reference_end-1]
            else:
                end_coords[uniq_key].append(line.reference_end-1)

            if uniq_key not in basecalls:
                basecalls[uniq_key] = {}
            if basecall not in basecalls[uniq_key]:
                basecalls[uniq_key][basecall] = 1
            else:
                basecalls[uniq_key][basecall] += 1
        return basecalls, end_coords


def main():
    args = supply_args()

    samfile = pysam.AlignmentFile(args.bam, "rb")
    primers = QiagenPrimers(args.bed).my_primers
    reader = vcfpy.Reader.from_path(args.vcf)
    all_mets = []
    for entry in reader:
        vrnt = (entry.CHROM, entry.POS-1)
        # If this is too slow, we can prepare primers list ahead to be localized to region, or even chrom.
        all_mets.append(VrntMetrics(samfile, vrnt, primers))

    reader.close()
    samfile.close()

    for entry in all_mets:
        print(entry.chrom)
        print(entry.pos)
        print(entry.assigned)
        print(entry.unassigned)


if __name__ == "__main__":
    main()
