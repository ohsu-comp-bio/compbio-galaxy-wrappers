class BedReader(object):
    """
    Simple class to ingest BED files and return a data structure as such:
    {chrom: [start1, stop1], [start2, stop2], ...}

    Input: filename
    """

    def __init__(self, filename):

        self.filename = open(filename, 'rU')
        self.bed_ints = self._create_bed()
        self.chrom_list = list(self.bed_ints)

    def _create_bed(self):

        """
        Create the structure of BED coordinates connected to chromosome
        identifiers.
        Return start and stop in 1-based coords.
        :return bed_ints:
        """
        bed_ints = {}
        with self.filename as bed:
            for interval in bed:
                interval = interval.rstrip('\n').split('\t')
                chrom = str(interval[0])
                # 0-based
                start = int(interval[1])+1
                # 1-based
                stop = int(interval[2])

                if chrom not in bed_ints:
                    bed_ints[chrom] = [[start, stop]]
                else:
                    bed_ints[chrom].append([start, stop])
        return bed_ints


    def split_coords(self):
        """
        Split out the intervals in to single coords.
        :return:
        """
        split_coords = {}
        for chrom in self.bed_ints:
            split_coords[chrom] = []
            for coords in self.bed_ints[chrom]:
                for coord in range(coords[0], coords[1]+1):
                    split_coords[chrom].append(coord)

        return split_coords


class ExtBedReader(object):
    """

    """
    def __init__(self, filename, header=False, hgnc=False, strand=False,
                 phase=False):
#        super(ExtBedReader, self).__init__(*args, **kwargs)
        self.filename = open(filename, 'rU')
        self.header = header
        self.hgnc = hgnc
        self.strand = strand
        self.phase = phase
        self.bed_ints = self._ext_bed_parse()

    def _ext_bed_parse(self):
        """
        Provide additional BED information for use.
        :return:
        """
        ext_bed_ints = {}
        with self.filename as bed:
            if self.header:
                next(bed)
            for interval in bed:
                interval = interval.rstrip('\n').split('\t')
                chrom = str(interval[0])
                if chrom.startswith('chr'):
                    chrom = chrom[3:]

                if chrom not in ext_bed_ints:
                    ext_bed_ints[chrom] = {}

                # 0-based
                start = int(interval[1]) + 1
                # 1-based
                stop = int(interval[2])

                new_key = (str(start), str(stop))
                ext_bed_ints[chrom][new_key] = {}
                ext_bed_ints[chrom][new_key]['start'] = start
                ext_bed_ints[chrom][new_key]['stop'] = stop

                if self.hgnc:
                    ext_bed_ints[chrom][new_key]['hgnc'] = interval[self.hgnc]
                if self.strand:
                    if interval[self.strand] == '-':
                        interval[self.strand] = 0
                    elif interval[self.strand] == '+':
                        interval[self.strand] = 1
                    ext_bed_ints[chrom][new_key]['strand'] = interval[
                            self.strand]
                if self.phase:
                    ext_bed_ints[chrom][new_key]['phase'] = interval[
                        self.phase]

        return ext_bed_ints

    def find_primer_coords(self, tsize=250):
        """
        Based on a given region length, and BED coordinates, find the
        primer coords.
        :return:
        """
        primer_coords = {}
        for chrom in self.bed_ints:
            primer_coords[chrom] = []
            for entry in self.bed_ints[chrom].values():
                rsize = entry['stop'] - entry['start'] + 1
                psize = tsize - rsize
                if entry['strand'] == 0:
                    # Theses are being flipped, so that the start > stop for
                    #  primers targeting reverse strand seqs.
                    pstop = entry['stop'] + 1
                    pstart = entry['stop'] + psize + 1
                elif entry['strand'] == 1:
                    pstart = entry['start'] - psize
                    pstop = entry['start'] - 1
                else:
                    raise ValueError("The strand should either be 0 or 1.")

                primer_coords[chrom].append([pstart, pstop])

        return primer_coords


    def get_primer_ends(self):
        """
        This is an extended version of find_primer_coords, where we pull
        down just the coordinates that are near the ends.
        {CHROM: {STRAND: (a, b, c, ...)}}
        :return:
        """
        primer_coords = {}
        for chrom in self.bed_ints.keys():
            primer_coords[chrom] = {}
            for coord in self.bed_ints[chrom].values():
                strand = str(coord['strand'])
                if strand not in primer_coords[chrom]:
                    primer_coords[chrom][strand] = []

                if coord['strand'] == 0:
                    pstop = coord['stop'] + 1
                elif coord['strand'] == 1:
                    pstop = coord['start'] - 1
                else:
                    raise ValueError("The strand should either be 0 or 1.")

                primer_coords[chrom][strand].append(pstop)

        return primer_coords
