import csv


class PerLocusRead(object):
    """
    Ingest a GATK DepthofCoverage per locus type file and prepare.
    Locus	Total_Depth	Average_Depth_sample	Depth_for_DNA-15-03359-1
    1:2496455	0	0.00	0
    """

    def __init__(self, infile):
        self.infile = csv.DictReader(open(infile, 'rU'), delimiter='\t')
        self.headers = self.infile.fieldnames
        self.perlocus = self._perlocus_fill()

    def __div__(self, other):

        div_perlocus = {}
        for coord in self.perlocus:
            try:
                div_perlocus[coord] = float(self.perlocus[coord]) / float(
                    other.perlocus[coord])
            except ZeroDivisionError:
                div_perlocus[coord] = 0.0

        return div_perlocus

    def _perlocus_fill(self):
        """
        Put information from per locus file in to dict of dicts.
        :return:
        """
        perlocus = {}
        for entry in self.infile:
            chrom = entry['Locus'].split(':')[0]
            coord = entry['Locus'].split(':')[1]
            perlocus[(chrom, coord)] = entry['Total_Depth']

        return perlocus


class ProbeQcRead(object):
    """
    Ingest a probe qc output file, as produced by the interval_qc_v2.py script.
    CHROM	START	STOP	REFSEQ	HGNC	AVGD	Q30	D2000	D500	D100	D20
    20	30946370	30946526			65.7	69.22	0.0	0.0	28.03	49.04
    NOTE: There is currently a total section at the bottom of this file.  We
    will be providing sample level metrics separately, so will need to account
    for this somehow.
    """

    def __init__(self, infile):

        if infile:
            self.infile = csv.DictReader(open(infile, 'rU'), delimiter='\t')
            self.headers = self.infile.fieldnames
            self.probeqc = self._probeqc_fill()

    def __div__(self, other):

        div_probeqc = {}
        for coord in self.probeqc:
            try:
                div_probeqc[coord] = float(self.probeqc[coord]['AVGD']) / \
                                 float(other.probeqc[coord]['AVGD'])
            except ZeroDivisionError:
                div_probeqc[coord] = 0.0

        return div_probeqc

    def _probeqc_fill(self):
        """
        Put information from ProbeQC file in to dict of dicts.
        :return:
        """
        probeqc = {}
        for entry in self.infile:
            chrom, start, stop = entry['CHROM'], entry['START'], entry['STOP']
            this_key = (chrom, start, stop)
            probeqc[this_key] = entry

        return probeqc
