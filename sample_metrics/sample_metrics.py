#!/usr/bin/env python3

# DESCRIPTION: Create sample level metrics to be passed to the CGD.  Metrics
#  are passed as a json dump.
# USAGE: sample_metrics.py -h
# CODED BY: John Letaw

import argparse
import json
# User libraries
from inputs import ProbeQcRead, AlignSummaryMetrics, GatkCountReads, MsiSensor, SamReader, GatkCollectRnaSeqMetrics
from inputs import FastQcRead

VERSION = '0.6.6'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    # Input files that will be parsed for data.
    parser.add_argument('--probeqc_after', type=ProbeQcRead,
                        help='Probe coverage QC after UMI deduplication metrics.')
    parser.add_argument('--probeqc_before', type=ProbeQcRead,
                        help='Probe coverage QC before UMI deduplication metrics.')
    parser.add_argument('--fastqc_r1', type=FastQcRead, help='FastQC stats for read 1.')
    parser.add_argument('--fastqc_r2', type=FastQcRead, help='FastQC stats for read 2.')
    parser.add_argument('--picard_summary', type=AlignSummaryMetrics, help='Picard alignment summary metrics file.')
    parser.add_argument('--gatk_coll_rnaseq_mets', type=GatkCollectRnaSeqMetrics,
                        help='GATK CollectRnaSeqMetrics file.')
    parser.add_argument('--gatk_count_reads_total', type=GatkCountReads,
                        help='Output from GATK4 CountReads with total read count.')
    parser.add_argument('--gatk_count_reads_ints', type=GatkCountReads,
                        help='Output from GATK4 CountReads with read count for reads overlapping targets.')
    parser.add_argument('--msi', type=MsiSensor, help='TSV file containing MSI results')

    parser.add_argument('--primers_bam', help='BAM file to calculate primer reads on target.')
    parser.add_argument('--primers_bed', help='BED file containing primer coordinates only.')

    # These just get attached to the final json output as-is.
    parser.add_argument('--json_in', nargs='*',
                        help='Arbitrary number of files to be included in sample metrics that are in json format.')

    parser.add_argument('--outfile', help='Output file with json string.')
    parser.add_argument('--outfile_new', help='Output file with new style json string.')
    parser.add_argument('--outfile_txt', help='Output file in human readable text format.')
    parser.add_argument('--workflow', help='Pass the Galaxy workflow name, if applicable.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    # Check to make sure if one gatk_count_reads option is set, both are set.
    if args.gatk_count_reads_total and not args.gatk_count_reads_ints:
        parser.error("Argument gatk_count_reads_total requires gatk_count_reads_ints.")
    if args.gatk_count_reads_ints and not args.gatk_count_reads_total:
        parser.error("Argument gatk_count_reads_ints requires gatk_count_reads_total.")

    # Check to make sure if one primers option is set, both are set.
    if args.primers_bam and not args.primers_bed:
        parser.error("Argument primers_bam requires primers_bed.")
    if args.primers_bed and not args.primers_bam:
        parser.error("Argument primers_bed requires primers_bam.")

    # Ensure at least one output option has been selected.
    if not args.outfile and not args.outfile_txt and not args.outfile_new:
        parser.error("You must specify one of outfile, outfile_new, or outfile_txt.")

    return args


class RawMetricCollector:
    """
    Gather all of the metrics from each input source, and allow it to be accessible from here.
    gatk_cr_total = total number of reads, as given by GATK4 CountReads
    gatk_cr_ints = total number of reads crossing target intervals, as given by GATK4 CountReads
    msi = from MSIsensor, TSV output file
    """
    def __init__(self, args):
        if args.gatk_count_reads_total:
            self.gatk_cr_total = args.gatk_count_reads_total.count
            self.gatk_cr_ints = args.gatk_count_reads_ints.count
        else:
            self.gatk_cr_total = None
            self.gatk_cr_ints = None

        if args.msi:
            self.msi = args.msi.msi
        else:
            self.msi = None

        if args.picard_summary:
            self.picard_summary = args.picard_summary.metrics
        else:
            self.picard_summary = None

        if args.gatk_coll_rnaseq_mets:
            self.gatk_coll_rnaseq_mets = args.gatk_coll_rnaseq_mets.metrics
        else:
            self.gatk_coll_rnaseq_mets = None

        if args.fastqc_r1:
            self.gc_pct_1 = args.fastqc_r1.gc_pct
        else:
            self.gc_pct_1 = None

        if args.fastqc_r2:
            self.gc_pct_2 = args.fastqc_r2.gc_pct
        else:
            self.gc_pct_2 = None

        if args.probeqc_before:
            self.probeqc_before = args.probeqc_before.probeqc
            self.probeqc_header_before = args.probeqc_before.headers
        else:
            self.probeqc_before = None
            self.probeqc_header_before = None

        if args.probeqc_after:
            self.probeqc_after = args.probeqc_after.probeqc
            self.probeqc_header_after = args.probeqc_after.headers
        else:
            self.probeqc_after = None
            self.probeqc_header_after = None

        if args.workflow:
            self.wf = args.workflow
        else:
            self.wf = None

        if args.primers_bam:
            self.primers_bam = SamReader(args.primers_bam, args.primers_bed).count
        else:
            self.primers_bam = None

        if args.json_in:
            self.json_in = args.json_in
            self.json_mets = self._json_in()
        else:
            self.json_mets = None

        # self._params_stdout()

    def _json_in(self):
        """

        :return:
        """
        json_mets = {}
        for filename in self.json_in:
            with open(filename, 'r') as myfile:
                for line in myfile:
                    for k, v in json.loads(line).items():
                        json_mets[str(k)] = v
        return json_mets

    def _params_stdout(self):
        """
        Print variable contents to stdout.
        :return:
        """
        print("MSI: {0}".format(self.msi))
        print("ProbeQC Before: {0}".format(self.probeqc_before))
        print("ProbeQC After: {0}".format(self.probeqc_after))
        print("FastQC R1 GC: {0}".format(self.gc_pct_1))
        print("FastQC R2 GC: {0}".format(self.gc_pct_2))
        print("GATK CountReads Total: {0}".format(self.gatk_cr_total))
        print("GATK CountReads Intervals: {0}".format(self.gatk_cr_ints))
        print("Picard: {0}".format(self.picard_summary))
        print("JSON Metrics: {0}".format(self.json_mets))
        print("Primers BAM: {0}".format(self.primers_bam))


class SampleMetrics:
    def __init__(self, raw_mets):
        self.raw_mets = raw_mets
        self.probeqc_before = self.raw_mets.probeqc_before
        if self.probeqc_before:
            self.total_cov_before = self._calc_cov(self.probeqc_before, 'AVGD')
            self.total_bp_before = self._calc_total_bp(self.probeqc_before)
            self.header_mets_before = self._metrics_from_probeqc_header(self.raw_mets.probeqc_header_before[5:],
                                                                        self.raw_mets.probeqc_before,
                                                                        self.total_bp_before)
        else:
            self.total_cov_before = None
            self.total_bp_before = None
            self.header_mets_before = None

        self.probeqc_after = self.raw_mets.probeqc_after
        if self.probeqc_after:
            self.total_cov_after = self._calc_cov(self.probeqc_after, 'AVGD')
            self.total_bp_after = self._calc_total_bp(self.probeqc_after)
            self.header_mets_after = self._metrics_from_probeqc_header(self.raw_mets.probeqc_header_after[5:],
                                                                       self.raw_mets.probeqc_after,
                                                                       self.total_bp_after)
        else:
            self.total_cov_after = None
            self.total_bp_after = None
            self.header_mets_before = None

        try:
            self.pumi = self._pumi()
        except:
            self.pumi = None
        try:
            self.on_target = self._add_on_target(self.raw_mets.picard_summary, self.total_cov_after)
        except:
            self.on_target = None
        try:
            self.on_primer_frag_count = self.raw_mets.primers_bam
        except:
            self.on_primer_frag_count = None
        try:
            self.on_primer_frag_count_pct = self._add_on_target(self.raw_mets.picard_summary,
                                                                self.on_primer_frag_count,
                                                                'PF_HQ_ALIGNED_READS')
        except:
            self.on_primer_frag_count_pct = None
        try:
            self.gatk_cr_on_target = self._gatk_cr_on_target()
        except:
            self.gatk_cr_on_target = None

        try:
            self.gatk_pct_mrna_bases = self.raw_mets.gatk_coll_rnaseq_mets['PCT_MRNA_BASES']
            self.gatk_pct_correct_strand_reads = self.raw_mets.gatk_coll_rnaseq_mets['PCT_CORRECT_STRAND_READS']
        except:
            self.gatk_pct_mrna_bases = None
            self.gatk_pct_correct_strand_reads = None

    def _gatk_cr_on_target(self):
        """
        Use GATK4 CountReads or CountReadsSpark to estimate on-target pct.
        :return:
        """
        if self.raw_mets.gatk_cr_total and self.raw_mets.gatk_cr_ints:
            return float(self.raw_mets.gatk_cr_ints) / (float(self.raw_mets.gatk_cr_total) + 0.0)

    def _pumi(self):
        """
        Create the PUMI metric, if applicable.
        :return:
        """
        if self.probeqc_before and self.probeqc_after:
            return self._calc_metric(self.total_cov_before, (self.total_cov_after * 100))

    def _metrics_from_probeqc_header(self, headers, probeqc, total_bp):
        """
        :param headers:
        :param probeqc:
        :param total_bp:
        :return:
        """
        sample_metrics = {}
        for label in headers:
            this_cov = self._calc_cov(probeqc, label)
            sample_metrics[label] = self._calc_metric(total_bp, this_cov)
        return sample_metrics

    @staticmethod
    def _calc_total_bp(probeqc):
        """
        Get the total number of base pairs covered by your targeted region set.
        :return:
        """
        total_bp = 0
        for line in probeqc.values():
            try:
                curr = int(line['STOP']) - int(line['START']) + 1.0
                total_bp += curr
            except ValueError:
                pass

        return total_bp

    @staticmethod
    def _calc_cov(probeqc, metric):
        """
        Calculate total coverage across sample.
        :param probeqc:
        :param metric:
        :return:
        """
        total_cov = 0
        for line in probeqc.values():
            try:
                curr_bp = int(line['STOP']) - int(line['START']) + 1.0
                curr_cov = curr_bp * float(line[metric])
                total_cov += curr_cov
            except ValueError:
                pass

        return total_cov

    @staticmethod
    def _calc_metric(bp, cov):
        """
        Calculate the sample level AVGD from total bp anc coverage.
        :param bp:
        :param cov:
        :return:
        """
        return '{:0.1f}'.format(cov / bp)

    @staticmethod
    def _add_on_target(picard, total_cov, total_lbl='PF_ALIGNED_BASES'):
        """
        Include the percent on target reads metric, mainly for amplicon assays.
        Also include the on primer frag count percentage.  Rename variables...
        :return:
        """
        if 'PAIR' in picard:
            pf_bases_aligned = int(picard['PAIR'][total_lbl])
        elif 'UNPAIRED' in picard:
            pf_bases_aligned = int(picard['UNPAIRED'][total_lbl])
        else:
            pf_bases_aligned = None

        if pf_bases_aligned and total_cov:
            on_target = str("{:.4}".format((total_cov * 100.0) / pf_bases_aligned))
        else:
            on_target = None

        return on_target


class MetricPrep(SampleMetrics):
    """
    Get the metrics we want from SampleMetrics, and prepare them for being written.
    """
    def __init__(self, raw_mets):
        super(MetricPrep, self).__init__(raw_mets)
        self.mets = self._gather_metrics()
        self.req_old = self._req_old()
        self.req_new = self._req_new()

    def _gather_metrics(self):
        """
        Old style metrics.  Will need to handle special json metrics separately.
        :return:
        'parentage_sites', 'parentage_disc', 'parentage_binom', 'parentage_confirmed'
        {"COUNT_0": 21916, "COUNT_1": 15287, "COUNT_10": 11012, "COUNT_100": 3505, "COUNT_1000": 187,
        "COUNT_10000": 9, "COUNT_100000": 0, "TPM_0": 22296, "TPM_0.01": 15366, "TPM_0.1": 15360, "TPM_1": 14323,
        "TPM_10": 9791, "TPM_100": 1505, "TPM_1000": 98}
        """
        mets = {'qthirty': self._get_avg_probeqc('Q30'),
                'averageDepth': self._get_avg_probeqc('AVGD'),
                'depthTen': self._get_avg_probeqc('D10'),
                'depthTwenty': self._get_avg_probeqc('D20'),
                'depthFifty': self._get_avg_probeqc('D50'),
                'depthOneHundred': self._get_avg_probeqc('D100'),
                'depthTwoHundredFifty': self._get_avg_probeqc('D250'),
                'depthFiveHundred': self._get_avg_probeqc('D500'),
                'depthSevenHundred': self._get_avg_probeqc('D700'),
                'depthOneThousand': self._get_avg_probeqc('D1000'),
                'depthTwelveHundredFifty': self._get_avg_probeqc('D1250'),
                'depthTwoThousand': self._get_avg_probeqc('D2000'),
                'allele_balance': self._reduce_sig(self._add_json_mets(lookin=self.raw_mets.json_mets,
                                                                       metric='allele_balance')),
                'allele_balance_het_count': self._add_json_mets(lookin=self.raw_mets.json_mets,
                                                                metric='allele_balance_het_count'),
                'gatk_cr_on_target': self._reduce_sig_pct(self.gatk_cr_on_target),
                'gatk_cr_total': self.raw_mets.gatk_cr_total,
                'gatk_cr_ints': self.raw_mets.gatk_cr_ints,
                'gatk_pct_mrna_bases': self._reduce_sig_pct(self.gatk_pct_mrna_bases),
                'gatk_pct_correct_strand_reads': self._reduce_sig_pct(self.gatk_pct_correct_strand_reads),
                'gc_pct_r1': self.raw_mets.gc_pct_1,
                'gc_pct_r2': self.raw_mets.gc_pct_2,
                'gender_check': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='gender_check'),
                'parentage_binom': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='parentage_binom'),
                'parentage_disc': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='parentage_disc'),
                'parentage_confirmed': self._add_json_mets(lookin=self.raw_mets.json_mets,
                                                           metric='parentage_confirmed'),
                'parentage_sites': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='parentage_sites'),
                'percentOnTarget': self._reduce_sig(self.on_target),
                'percentUmi': self.pumi,
                'rna_count_zero': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='COUNT_0'),
                'rna_count_one': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='COUNT_1'),
                'rna_count_ten': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='COUNT_10'),
                'rna_count_onehundred': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='COUNT_100'),
                'rna_count_onethousand': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='COUNT_1000'),
                'rna_count_tenthousand': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='COUNT_10000'),
                'rna_count_hundredthousand': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='COUNT_100000'),
                'rna_tpm_zero': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='TPM_0'),
                'rna_tpm_hundredth': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='TPM_0.01'),
                'rna_tpm_tenth': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='TPM_0.1'),
                'rna_tpm_one': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='TPM_1'),
                'rna_tpm_ten': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='TPM_10'),
                'rna_tpm_onehundred': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='TPM_100'),
                'rna_tpm_onethousand': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='TPM_1000'),
                'tmb': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='tmb'),
                'msi_pct': self._add_json_mets(lookin=self.raw_mets.msi, metric='somatic_pct'),
                'msi_sites': self._add_json_mets(lookin=self.raw_mets.msi, metric='total_sites'),
                'msi_somatic_sites': self._add_json_mets(lookin=self.raw_mets.msi, metric='somatic_sites'),

                'total_on_target_transcripts': self.on_primer_frag_count,
                'total_on_target_transcripts_pct': self.on_primer_frag_count_pct
                }

        return mets

    @staticmethod
    def _add_json_mets(lookin, metric):
        """
        Currently the following metrics are sent in standalone.  Need to think about ways to do this better.
        allele_balance
        tmb
        :return:
        """
        if lookin:
            if metric in lookin:
                return str(lookin[metric])
        return None

    @staticmethod
    def _reduce_sig(metric):
        """
        Get rid of digits that are not significant.
        :return:
        """
        if metric:
            return "{:.2f}".format(float(metric))

    @staticmethod
    def _reduce_sig_pct(metric):
        """
        Get rid of digits that are not significant, and provide in percent format.
        :return:
        """
        if metric:
            return "{:.2f}".format(float(metric)*100)

    def _get_avg_probeqc(self, label):
        """
        Retrieve average metric from the probeqc_after dict.  All ProbeQCs are after, unless they
        are being run on UMI-containing assays.  If this is the case, the probeqc_before dict
        captures metrics before UMI deduplication occurs.
        :return:
        """
        try:
            this_cov = self._calc_cov(self.probeqc_after, label)
            return self._calc_metric(self.total_bp_after, this_cov)
        except:
            return None

    @staticmethod
    def _req_old():
        """
        Based on test name, list which metrics should be provided.
        :return:
        """
        return {'QIAseq_V3_RNA': ['qthirty', 'averageDepth', 'percentUmi'],
                'TruSightOne': ['qthirty', 'averageDepth', 'depthTwoHundredFifty', 'depthTwenty',
                                'depthOneHundred', 'percentOnTarget', 'depthTen', 'depthFifty'],
                'TruSightOneV2_5': ['qthirty', 'averageDepth', 'depthTwoHundredFifty', 'depthTwenty',
                                    'depthOneHundred', 'percentOnTarget', 'depthTen', 'depthFifty'],
                'AgilentCRE_V1': ['qthirty', 'averageDepth', 'depthTwoHundredFifty', 'depthTwenty',
                                  'depthOneHundred', 'percentOnTarget', 'depthTen', 'depthFifty'],
                'QIAseq_V3_HEME2': ['qthirty', 'averageDepth', 'depthTwoHundredFifty', 'depthTwelveHundredFifty',
                                    'depthOneHundred', 'percentOnTarget', 'depthSevenHundred', 'percentUmi'],
                'QIAseq_V3_STP3': ['qthirty', 'averageDepth', 'depthTwoHundredFifty', 'depthTwelveHundredFifty',
                                   'depthOneHundred', 'percentOnTarget', 'depthSevenHundred', 'percentUmi'],
                'TruSeq_RNA_Exome_V1-2': ['qthirty'],
                'QIAseq_V3_HOP': ['qthirty', 'averageDepth', 'depthTwoHundredFifty', 'depthTwenty',
                                  'depthOneHundred', 'percentOnTarget', 'depthTen', 'depthFifty', 'percentUmi'],
                'QIAseq_V3_HOP2': ['qthirty', 'averageDepth', 'depthTwoHundredFifty', 'depthTwenty',
                                   'depthOneHundred', 'percentOnTarget', 'depthTen', 'depthFifty', 'percentUmi']
                }

    @staticmethod
    def _req_new():
        """
        Based on test name, list which metrics should be provided.
        :return:
        """
        return {'QIAseq_V3_RNA': ['total_on_target_transcripts', 'total_on_target_transcripts_pct'],
                'TruSightOne': ['gc_pct_r1', 'gc_pct_r2', 'gender_check'],
                'TruSightOneV2_5': ['gc_pct_r1', 'gc_pct_r2', 'gender_check'],
                'AgilentCRE_V1': ['parentage_sites', 'parentage_disc', 'parentage_binom', 'parentage_confirmed',
                                  'gc_pct_r1', 'gc_pct_r2', 'gender_check'],
                'QIAseq_V3_HEME2': [],
                'QIAseq_V3_STP3': ['msi_sites', 'msi_somatic_sites', 'msi_pct', 'tmb'],
                'TruSeq_RNA_Exome_V1-2': ['total_on_target_transcripts', 'gatk_pct_mrna_bases',
                                          'gatk_pct_correct_strand_reads', 'rna_count_zero', 'rna_count_one',
                                          'rna_count_ten', 'rna_count_onehundred', 'rna_count_onethousand',
                                          'rna_count_tenthousand', 'rna_count_hundredthousand', 'rna_tpm_zero',
                                          'rna_tpm_hundredth', 'rna_tpm_tenth', 'rna_tpm_one', 'rna_tpm_ten',
                                          'rna_tpm_onehundred', 'rna_tpm_onethousand'],
                'QIAseq_V3_HOP': ['allele_balance', 'allele_balance_het_count'],
                'QIAseq_V3_HOP2': ['allele_balance', 'allele_balance_het_count']
                }


class Writer:
    """

    """
    def __init__(self, mets):
        self.mets = mets
        self.wf = self.mets.raw_mets.wf

    def write_cgd_new(self, filename):
        """
        Provide new metric style for CGD import.
        {
        "sampleRunMetrics": [
            {
                "metric": "total_on_target_reads",
                "value": 1230411
            },
            {
                "metric": "percent_on_target_reads",
                "value": 0.91
            }
        ],
        "geneMetrics": [
            {
                "gene": "ASXL1",
                "metric": "total_on_target_reads",
                "value": 469012
            },
            {
                "gene": "BRCA1",
                "metric": "total_on_target_reads",
                "value": 362330
            }
        ]
        }
        :return:
        """
        to_write = {'sampleRunMetrics': [], 'geneMetrics': []}
        for metric, val in self.mets.mets.items():
            if metric in self.mets.req_new[self.wf]:
                if val:
                    metric_dict = {'metric': str(metric), 'value': str(val)}
                    to_write['sampleRunMetrics'].append(metric_dict)

        with open(filename, 'w') as jwrite:
            json.dump(to_write, jwrite)

    def write_cgd_old(self, filename):
        """
        Provide old metric style for CGD import.
        {
        "depthSevenHundred": "95.6",
        "depthTwelveHundredFifty": "66.8",
        "percentOnTarget": "61.03"
        }
        """
        to_write = {}
        for metric, val in self.mets.mets.items():
            if metric in self.mets.req_old[self.wf]:
                to_write[str(metric)] = str(val)

        with open(filename, 'w') as jwrite:
            json.dump(to_write, jwrite)

    def write_to_text(self, filename):
        """
        Write metrics to a text file, mainly to be viewed in Galaxy.
        :return:
        """
        with open(filename, 'w') as to_write:
            for metric, val in self.mets.mets.items():
                to_write.write("{}: {}\n".format(metric, val))


def main():
    args = supply_args()
    raw_mets = RawMetricCollector(args)
    samp_mets = MetricPrep(raw_mets)
    writer = Writer(samp_mets)
    if args.outfile_txt:
        writer.write_to_text(args.outfile_txt)
    if args.outfile:
        writer.write_cgd_old(args.outfile)
    if args.outfile_new:
        writer.write_cgd_new(args.outfile_new)


if __name__ == "__main__":
    main()
