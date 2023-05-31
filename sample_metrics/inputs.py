import csv
import os
import pysam
import vcfpy
import sys
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess


class GatkCollectRnaSeqMetrics:
    """
    Grab data from files in this format.
    ## denotes section headers
    # is more or less a comment?
    """
    def __init__(self, filename):
        self.filename = filename
        self.metrics = self._fill_metrics()

    def _fill_metrics(self):
        """
        This will only grab the METRICS section values.
        :return:
        """
        metrics = {}
        with open(self.filename, 'r') as outfile:
            for line in outfile:
                if line.startswith('PF_BASES'):
                    headers = line.rstrip('\n').split('\t')
                    raw_metrics = outfile.readline().rstrip('\n').split('\t')
                    for i in range(len(headers)):
                        metrics[headers[i]] = raw_metrics[i]
                    return metrics
        raise ImportError("Can't find METRICS CLASS data from the input CollectRnaSeqMetrics file.")


class SamReader:
    """
    Read a BAM file, get stuff we need from it, like total counts of all reads.
    """
    def __init__(self, filename, bedfile):
        self.filename = filename
        self.bedfile = bedfile
        try:
            self.count = self._pysam_get_target_count()
        except:
            self.count = self._run_cmd(self._get_target_count_cmd())

    @staticmethod
    def _run_cmd(cmd):
        """
        Run command.
        """
        print('Running the following command:')
        print('\t'.join(cmd))

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        return stdout

    def _pysam_get_target_count(self):
        """

        :return:
        """
        return int(pysam.view(self.filename, "-L", self.bedfile, "-c").strip())

    def _get_target_count_cmd(self):
        """
        From a BAM file, use samtools view to count the number of reads with targets in primers_bed.
        samtools-1.3.1 view -L primers_only.bed test.bam -c
        Core dumps occur trying to run below pysam command.  Command work great for several versions of samtools on command
        line, so we're doing subprocess until better solution.
        :return:
        """
        if self.bedfile:
            cmd = ['samtools', 'view', '-L', self.bedfile, '-c', self.filename]
        else:
            cmd = ['samtools', 'view', '-c', self.filename]
        return cmd


class MsiSensor:
    """
    Output comes from MSIsensor.  Looks like:

    """
    def __init__(self, filename):
        self.fh = open(filename, 'r')
        self.msi = self.msi_parse()
        self.fh.close()

    def msi_parse(self):
        """
        Get the MSI results from the TSV file.
        Total_Number_of_Sites	Number_of_Somatic_Sites	%
        230	60	26.09

        :return:
        """
        with self.fh as mymsi:
            for line in mymsi:
                line = line.rstrip('\n').split('\t')
                if not line[0].startswith('Total') and valid:
                    msi = {'total_sites': int(line[0]),
                           'somatic_sites': int(line[1]),
                           'somatic_pct': float(line[2])}
                    return msi
                elif line[0].startswith('Total') and len(line) == 3:
                    valid = True
                else:
                    raise Exception("Invalid formatting in MSIsensor input file, please check.")


class GatkCountReads:
    """
    Output files contain count at top of file.
    """
    def __init__(self, filename):
        self.fh = open(filename, 'r')
        self.count = self._get_count()
        self.fh.close()


    def _get_count(self):
        """
        File contains a string-representation of an integer.  Attempt to coerce
        to integer, to avoid most issues with bad values.
        :return:
        """
        return int(self.fh.readline().strip())


class PerLocusRead:
    """
    Ingest a GATK DepthofCoverage per locus type file and prepare.
    Locus	Total_Depth	Average_Depth_sample	Depth_for_DNA-15-03359-1
    1:2496455	0	0.00	0
    """

    def __init__(self, infile):
        self.infile = csv.DictReader(open(infile, 'r'), delimiter='\t')
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


class ProbeQcRead:
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
            self.infile = csv.DictReader(open(infile, 'r'), delimiter='\t')
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

        if probeqc != {}:
            return probeqc
        else:
            return None


class AlignSummaryMetrics:
    """
    Parse AlignmentSummaryMetrics output from Picard.
    """
    def __init__(self, infile):
        self.filename = open(infile, 'r')
        self.header, self.metrics = self._parse_metrics()

    def _parse_metrics(self):
        """
        Turn the text file of metrics in to a structure.
        :return:
        """
        header = []
        metrics = {}

        with self.filename as myfile:
            for line in myfile:
                if line.startswith('#'):
                    header.append(line)
                elif line.startswith('CATEGORY'):
                    titles = line.rstrip('\n').split('\t')
                elif line == '\n':
                    pass
                else:
                    line = line.rstrip('\n').split('\t')
                    metrics[line[0]] = {}
                    for i in range(len(line)):
                        if titles[i] not in metrics[line[0]]:
                            metrics[line[0]][titles[i]] = line[i]

        return header, metrics


class FastQcRead:
    """
    ##FastQC	0.11.8
    >>Basic Statistics	pass
    #Measure	Value
    Filename	19KD-004H8888-1_S1_R2_001_fastq_gz.gz
    File type	Conventional base calls
    Encoding	Sanger / Illumina 1.9
    Total Sequences	1000000
    Sequences flagged as poor quality	0
    Sequence length	32-151
    %GC	46
    >>END_MODULE
    """
    def __init__(self, infile):
        self.infile = infile
        self.gc_pct = self._get_header_val(metric='%GC')
        self.seq_cnt = self._get_header_val(metric='Total Sequences')

    def _get_header_val(self, metric):
        with open(self.infile, 'r') as myfile:
            for line in myfile:
                if line.startswith(metric):
                    val = line.rstrip('\n').split('\t')[1]
                    return val


class VcfRead:
    """
    """
    def __init__(self, infile):
        self.reader = vcfpy.Reader.from_path(infile)
        self.count = self._get_count()

    def _get_count(self):
        records = {}
        for record in self.reader:
            var_id = '{}:{}{}>{}'.format(record.CHROM, record.POS, record.REF, record.ALT)
            if var_id not in records:
                records[var_id] = record
        return len(records)


class GatkDepthOfCoverageRead:
    """
    Ingest a GATK DepthofCoverage sample type file
    """

    def __init__(self, infile):
        self.infile = csv.DictReader(open(infile, 'r'), delimiter='\t')
        self.headers = self.infile.fieldnames
        self.sample_mets = self._sample_fill()

    def _sample_fill(self):
        """
        Put information from per locus file in to dict of dicts.
        :return:
        """
        for entry in self.infile:
            for e in entry:
                if e == 'sample':
                    return entry
