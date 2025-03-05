"""
Create sample level metrics to be passed to the CGD.  Metrics are passed as a json dump.

"""

import argparse
import json
# User libraries
from inputs import (ProbeQcRead, PerLocusRead, AlignSummaryMetrics, GatkDepthOfCoverageRead, GatkCountReads,
                    MsiSensor, SamReader, GatkCollectRnaSeqMetrics)
from inputs import FastQcRead


# TODO: Update version number in ./README.md
VERSION = '0.9.3'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--metrics', required=False,
                        default="qthirty, averageDepth",
                        help='List of metrics to be provided.')
    # xml inputs params type=text will pass metrics args
    # as 'qthirty averageDepth depthTwoHundredFifty' (even with nargs=*)
    # best to split to leave as string and split later

    # Input files that will be parsed for data.
    parser.add_argument('--probeqc_after', required=False, type=ProbeQcRead,
                        help='Probe coverage QC after UMI deduplication metrics.')
    parser.add_argument('--probeqc_before', required=False, type=ProbeQcRead,
                        help='Probe coverage QC before UMI deduplication metrics.')

    parser.add_argument('--fastqc_r1', type=FastQcRead, help='FastQC stats for read 1.')
    parser.add_argument('--fastqc_r2', type=FastQcRead, help='FastQC stats for read 2.')

    parser.add_argument('--picard_summary', type=AlignSummaryMetrics,
                        help='Picard alignment summary metrics file.')
    parser.add_argument('--picard_summary_umi', type=AlignSummaryMetrics,
                        help='Picard alignment summary metrics file taken after umi deduplication.')

    parser.add_argument('--gatk_depth_cov_prop', type=GatkDepthOfCoverageRead,
                        help='GATK DepthOfCoverage file.')
    parser.add_argument('--gatk_depth_cov_cnts', type=PerLocusRead,
                        help='GATK DepthOfCoverage locus file.')
    parser.add_argument('--gatk_coll_rnaseq_mets', type=GatkCollectRnaSeqMetrics,
                        help='GATK CollectRnaSeqMetrics file.')
    parser.add_argument('--gatk_count_reads_total', type=GatkCountReads,
                        help='Output from GATK4 CountReads with total read count.')
    parser.add_argument('--gatk_count_reads_ints', type=GatkCountReads,
                        help='Output from GATK4 CountReads with read count for reads overlapping targets.')

    parser.add_argument('--msi', type=MsiSensor, help='TSV file containing MSI results')

    parser.add_argument('--primers_bam', help='BAM file to calculate primer reads on target.')
    parser.add_argument('--primers_bed', help='BED file containing primer coordinates only.')

    parser.add_argument('--blia_pre', help='JSON from Ding correlation subtyping, pre-normalization.')
    parser.add_argument('--blia_post', help='JSON from Ding correlation subtyping, post-normalization.')

    parser.add_argument('--dragen_metrics', help='JSON file containing metrics produced by DRAGEN')
    parser.add_argument('--dragen_qc', help='CSV file containing fastqc data produced by DRAGEN')

    # These just get attached to the final json output as-is.
    parser.add_argument('--json_in', nargs='*',
                        help='Arbitrary number of files to be included in sample metrics that are in json format.')

    parser.add_argument('--outfile', help='Output file with json string.')
    parser.add_argument('--outfile_txt', help='Output file in human readable text format.')
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
    if not args.outfile and not args.outfile_txt:
        parser.error("You must specify outfile or outfile_txt.")

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

        if args.picard_summary_umi:
            self.picard_summary_umi = args.picard_summary_umi.metrics
        else:
            self.picard_summary_umi = None

        if args.gatk_depth_cov_prop:
            self.gatk_depth_cov_prop = args.gatk_depth_cov_prop
        else:
            self.gatk_depth_cov_prop = None

        if args.gatk_depth_cov_cnts:
            self.gatk_depth_cov_cnts = args.gatk_depth_cov_cnts
        else:
            self.gatk_depth_cov_cnts = None

        if args.gatk_coll_rnaseq_mets:
            self.gatk_coll_rnaseq_mets = args.gatk_coll_rnaseq_mets.metrics
        else:
            self.gatk_coll_rnaseq_mets = None

        if args.fastqc_r1:
            self.gc_pct_1 = args.fastqc_r1.gc_pct
            self.fastqc_seq_count = args.fastqc_r1.seq_cnt
        else:
            self.gc_pct_1 = None
            self.fastqc_seq_count = None

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

        if args.primers_bam:
            self.primers_bam = SamReader(args.primers_bam, args.primers_bed).count
        else:
            self.primers_bam = None

        if args.blia_pre:
            self.blia_pre = self._json_in([args.blia_pre])
        else:
            self.blia_pre = {'blia': None, 'blis': None, 'lar': None, 'mes': None}

        if args.blia_post:
            self.blia_post = self._json_in([args.blia_post])
        else:
            self.blia_post = {'blia': None, 'blis': None, 'lar': None, 'mes': None}

        if args.json_in:
            self.json_mets = self._json_in(args.json_in)
        else:
            self.json_mets = None

        if args.dragen_metrics:
            self.dragen_metrics = self._dragen_json(args.dragen_metrics)
        else:
            self.dragen_metrics = None

        if args.dragen_qc:
            self.dragen_qc = args.dragen_qc
        else:
            self.dragen_qc = None

    def _dragen_json(self, json_in):
        """
         Collect metrics from DRAGEN metrics JSON
         :return:
        """
        with open(json_in, "r") as metric_handle:
            dragen_mets = {}
            for line in metric_handle:
                line = line.strip().strip(",").replace("\"", "").replace(":", "")
                if "{" not in line and "}" not in line:
                    dragen_mets[line.split()[0]] = line.split()[1]
        return dragen_mets
    @staticmethod
    def _json_in(json_in):
        """
        Get all of the incoming JSON metrics.
        :return:
        """
        json_mets = {}
        for filename in json_in:
            with open(filename, 'r') as myfile:
                for line in myfile:
                    try:
                        for k, v in json.loads(line).items():
                            json_mets[str(k)] = v
                    except AttributeError:
                        for k, v in json.loads(line)[0].items():
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

        if self.raw_mets.dragen_qc:
            self.gc_r1, self.gc_r2 = self._dragen_gc_calc(self.raw_mets.dragen_qc)

        try:
            self.gatk_depth_cov_prop = self.raw_mets.gatk_depth_cov_prop
        except:
            self.gatk_depth_cov_prop = None

        try:
            self.gatk_depth_cov_cnts = self.raw_mets.gatk_depth_cov_cnts
        except:
            self.gatk_depth_cov_cnts = None

        try:
            self.pumi = self._pumi()
        except:
            self.pumi = None
        try:
            self.on_target = self._add_on_target(self.raw_mets.picard_summary, self.raw_mets.fastqc_seq_count,
                                                 total_lbl='PF_READS_ALIGNED')
            self.on_target_after = self._add_on_target(self.raw_mets.picard_summary_umi, self.raw_mets.fastqc_seq_count,
                                                 total_lbl='PF_READS_ALIGNED')
        except:
            self.on_target = None
            self.on_target_after = None
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

        try:
            self.blia_pre_mets = self._top_two_diff(self.raw_mets.blia_pre)
            self.blia_post_mets = self._top_two_diff(self.raw_mets.blia_post)
            self.assess_blia = self._assess_blia()
        except:
            self.blia_pre_mets = {'best': None, 'second': None, 'diff': None}
            self.blia_post_mets = {'best': None, 'second': None, 'diff': None}
            self.assess_blia = None

    def _dragen_gc_calc(self, qc_data):
        """
        Calculate GC content for R1 and R2 from DRAGEN QC data
        """
        r1_total, r2_total = 0, 0
        r1_gc, r2_gc = 0, 0

        # sum bases from DRAGEN FastQC output file
        with open(qc_data) as qc_handle:
            for line in qc_handle:
                if "POSITIONAL BASE CONTENT" in line:
                    num = int(line.split(',')[3])
                    if "Read1" in line:
                        r1_total += num
                        if "G Bases" in line or "C Bases" in line:
                            r1_gc += num
                    if "Read2" in line:
                        r2_total += num
                        if "G Bases" in line or "C Bases" in line:
                            r2_gc += num

        # compute GC content for R1 and R2
        gc_r1_pct = round(r1_gc * 100 / r1_total)
        gc_r2_pct = round(r2_gc * 100 / r2_total)

        return gc_r1_pct, gc_r2_pct

    def _blia_blis_map(self, name):
        """
        Map between numbers and strings.
        :return:
        """
        mapping = {'blia': '0',
                   'blis': '1',
                   'lar': '2',
                   'mes': '3'}
        return mapping[name]

    def _assess_blia(self, thresh=0.1):
        """
        If the difference between the top two scores is >=0.1 for both pre and post, and the two
        top types are the same, this is a reportable subtyping.
        :return:
        """
        if (self.blia_pre_mets['diff'] >= thresh
                and self.blia_post_mets['diff'] >= thresh
                and self.blia_pre_mets['best'] == self.blia_post_mets['best']
                and self.blia_pre_mets['second'] == self.blia_post_mets['second']
                and (self.blia_post_mets['best'] == '0' or self.blia_post_mets['best']) == '1'):
            return '1'
        return '0'

    def _top_two_diff(self, scores):
        """
        Get the difference between the top two correlation scores.
        :return:
        """
        sort_scores = dict(sorted((val, key) for (key, val) in scores.items()))
        diff = list(sort_scores.keys())[-1] - list(sort_scores.keys())[-2]
        best = list(sort_scores.values())[-1]
        second = list(sort_scores.values())[-2]
        return {'diff': diff,
                'best': self._blia_blis_map(best),
                'second': self._blia_blis_map(second)
                }

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
        if self.raw_mets.picard_summary and self.raw_mets.picard_summary_umi:
            before = int(self.raw_mets.picard_summary['FIRST_OF_PAIR']['PF_ALIGNED_BASES'])
            after = int(self.raw_mets.picard_summary_umi['FIRST_OF_PAIR']['PF_ALIGNED_BASES'])
            return self._calc_metric(before, (after * 100))

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
        if 'FIRST_OF_PAIR' in picard:
            pf_bases_aligned = int(picard['FIRST_OF_PAIR'][total_lbl])
        elif 'UNPAIRED' in picard:
            pf_bases_aligned = int(picard['UNPAIRED'][total_lbl])
        else:
            pf_bases_aligned = None

        if pf_bases_aligned and total_cov:
            if total_lbl == 'PF_HQ_ALIGNED_READS':
                on_target = str("{:.4}".format((int(total_cov) / pf_bases_aligned) * 100.0))
            else:
                on_target = str("{:.4}".format(pf_bases_aligned / (int(total_cov)) * 100.0))
        else:
            on_target = None

        return on_target


class MetricPrep(SampleMetrics):
    """
    Get the metrics we want from SampleMetrics, and prepare them for being written.
    """
    def __init__(self, raw_mets, metrics):
        super(MetricPrep, self).__init__(raw_mets)
        self.metrics_req = metrics
        self.mets = self._gather_metrics()

    def _gather_metrics(self):
        """
        Old style metrics.  Will need to handle special json metrics separately.
        :return:
        'parentage_sites', 'parentage_disc', 'parentage_binom', 'parentage_confirmed'
        {"COUNT_0": 21916, "COUNT_1": 15287, "COUNT_10": 11012, "COUNT_100": 3505, "COUNT_1000": 187,
        "COUNT_10000": 9, "COUNT_100000": 0, "TPM_0": 22296, "TPM_0.01": 15366, "TPM_0.1": 15360, "TPM_1": 14323,
        "TPM_10": 9791, "TPM_100": 1505, "TPM_1000": 98}
        """
        if self.raw_mets.dragen_metrics:
            mets = self.raw_mets.dragen_metrics
            if self.raw_mets.dragen_qc:
                mets['dragen_gc_pct_r1'] = self.gc_r1
                mets['dragen_gc_pct_r2'] = self.gc_r2
        else:
            mets = {'qthirty': self._get_avg_probeqc('Q30', self.probeqc_after, self.total_bp_after),
                    'qthirty_before': self._get_avg_probeqc('Q30', self.probeqc_before, self.total_bp_before),
                    'averageDepth': self._get_avg_probeqc('AVGD', self.probeqc_after, self.total_bp_after),
                    'averageDepth_before': self._get_avg_probeqc('AVGD', self.probeqc_before, self.total_bp_before),
                    'depthTen': self._get_avg_probeqc('D10', self.probeqc_after, self.total_bp_after),
                    'depthTen_before': self._get_avg_probeqc('D10', self.probeqc_before, self.total_bp_before),
                    'depthTwenty': self._get_avg_probeqc('D20', self.probeqc_after, self.total_bp_after),
                    'depthTwenty_before': self._get_avg_probeqc('D20', self.probeqc_before, self.total_bp_before),
                    'depthFifty': self._get_avg_probeqc('D50', self.probeqc_after, self.total_bp_after),
                    'depthFifty_before': self._get_avg_probeqc('D50', self.probeqc_before, self.total_bp_before),
                    'depthOneHundred': self._get_avg_probeqc('D100', self.probeqc_after, self.total_bp_after),
                    'depthOneHundred_before': self._get_avg_probeqc('D100', self.probeqc_before, self.total_bp_before),
                    'depthTwoHundredFifty': self._get_avg_probeqc('D250', self.probeqc_after, self.total_bp_after),
                    'depthTwoHundredFifty_before':
                        self._get_avg_probeqc('D250', self.probeqc_before, self.total_bp_before),
                    'depthFiveHundred': self._get_avg_probeqc('D500', self.probeqc_after, self.total_bp_after),
                    'depthFiveHundred_before':
                        self._get_avg_probeqc('D500', self.probeqc_before, self.total_bp_before),
                    'depthSevenHundred':
                        self._get_avg_probeqc('D700', self.probeqc_after, self.total_bp_after),
                    'depthSevenHundred_before':
                        self._get_avg_probeqc('D700', self.probeqc_before, self.total_bp_before),
                    'depthOneThousand':
                        self._get_avg_probeqc('D1000', self.probeqc_after, self.total_bp_after),
                    'depthOneThousand_before':
                        self._get_avg_probeqc('D1000', self.probeqc_before, self.total_bp_before),
                    'depthTwelveHundredFifty':
                        self._get_avg_probeqc('D1250', self.probeqc_after, self.total_bp_after),
                    'depthTwelveHundredFifty_before':
                        self._get_avg_probeqc('D1250', self.probeqc_before, self.total_bp_before),
                    'depthTwoThousand': self._get_avg_probeqc('D2000', self.probeqc_after, self.total_bp_after),
                    'depthTwoThousand_before':
                        self._get_avg_probeqc('D2000', self.probeqc_before, self.total_bp_before),
                    'uniformity_of_coverage':
                        self._get_cov_uniformity(self._get_avg_probeqc('AVGD', self.probeqc_after,
                                                                       self.total_bp_after)),
                    'pf_reads': self.raw_mets.picard_summary['PAIR']['PF_READS'] if self.raw_mets.picard_summary else None,
                    'pf_reads_after': self.raw_mets.picard_summary_umi['PAIR']['PF_READS'] if self.raw_mets.picard_summary_umi else None,
                    'pct_pf_reads_aligned': self.raw_mets.picard_summary['PAIR']['PCT_PF_READS_ALIGNED'] if self.raw_mets.picard_summary else None,
                    'pct_pf_reads_aligned_after': self.raw_mets.picard_summary_umi['PAIR']['PCT_PF_READS_ALIGNED'] if self.raw_mets.picard_summary_umi else None,
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
                    'bio_sex_check': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='bio_sex_check'),
                    'homozygosity_flag': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='homozygosity_flag'),
                    'parentage_binom': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='parentage_binom'),
                    'parentage_disc': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='parentage_disc'),
                    'parentage_confirmed': self._add_json_mets(lookin=self.raw_mets.json_mets,
                                                               metric='parentage_confirmed'),
                    'parentage_sites': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='parentage_sites'),
                    'percentOnTarget': self._reduce_sig(self.on_target),
                    'percentOnTarget_after': self._reduce_sig(self.on_target_after),
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
                    'tmb': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='tmb') if
                    self._add_json_mets(lookin=self.raw_mets.json_mets, metric='tmb') else 0.0001,
                    'msi_pct': self._add_json_mets(lookin=self.raw_mets.msi, metric='somatic_pct'),
                    'msi_sites': self._add_json_mets(lookin=self.raw_mets.msi, metric='total_sites'),
                    'msi_somatic_sites': self._add_json_mets(lookin=self.raw_mets.msi, metric='somatic_sites'),
                    'blia_pre_best': self.blia_pre_mets['best'],
                    'blia_pre_second': self.blia_pre_mets['second'],
                    'blia_pre_diff': self._reduce_sig(self.blia_pre_mets['diff']),
                    'blia_post_best': self.blia_post_mets['best'],
                    'blia_post_second': self.blia_post_mets['second'],
                    'blia_post_diff': self._reduce_sig(self.blia_post_mets['diff']),
                    'blia_reportable': self.assess_blia,
                    'blia_raw_pre': self.raw_mets.blia_pre['blia'],
                    'blis_raw_pre': self.raw_mets.blia_pre['blis'],
                    'lar_raw_pre': self.raw_mets.blia_pre['lar'],
                    'mes_raw_pre': self.raw_mets.blia_pre['mes'],
                    'blia_raw_post': self.raw_mets.blia_post['blia'],
                    'blis_raw_post': self.raw_mets.blia_post['blis'],
                    'lar_raw_post': self.raw_mets.blia_post['lar'],
                    'mes_raw_post': self.raw_mets.blia_post['mes'],
                    'total_on_target_transcripts': self.on_primer_frag_count,
                    'total_on_target_transcripts_pct': self.on_primer_frag_count_pct,
                    'forced_calls_above': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='forced_calls_above'),
                    'forced_calls_below': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='forced_calls_below'),
                    'y_ploidy_check': self._add_json_mets(lookin=self.raw_mets.json_mets, metric='y_ploidy_check'),
                    'cnv_median_segment_mad_cn':
                        self._add_json_mets(lookin=self.raw_mets.json_mets, metric='cnv_median_segment_mad_cn'),
                    'q30_bases_pct': self.raw_mets.dragen_metrics.q30 if self.raw_mets.dragen_metrics else None,
                    'average_alignment_coverage_over_target_region':
                        self.raw_mets.dragen_metrics.avg_depth if self.raw_mets.dragen_metrics else None,
                    'pct_of_target_region_with_coverage_10x_inf':
                        self.raw_mets.dragen_metrics.depth10 if self.raw_mets.dragen_metrics else None,
                    'pct_of_target_region_with_coverage_20x_inf':
                        self.raw_mets.dragen_metrics.depth20 if self.raw_mets.dragen_metrics else None,
                    'pct_of_target_region_with_coverage_50x_inf':
                        self.raw_mets.dragen_metrics.depth50 if self.raw_mets.dragen_metrics else None,
                    'pct_of_target_region_with_coverage_100x_inf':
                        self.raw_mets.dragen_metrics.depth100 if self.raw_mets.dragen_metrics else None,
                    'aligned_reads_in_target_region_pct':
                        self.raw_mets.dragen_metrics.pct_on_target if self.raw_mets.dragen_metrics else None,
                    'dragen_gc_pct_r1': self.raw_mets.dragen_qc.gc_r1 if self.raw_mets.dragen_metrics else None,
                    'dragen_gc_pct_r2': self.raw_mets.dragen_qc.gc_r2 if self.raw_mets.dragen_metrics else None,
                    'number_of_large_roh_gt_eq_3000000':
                        self.raw_mets.dragen_metrics.roh if self.raw_mets.dragen_metrics else None,
                    'ploidy_estimation': self.raw_mets.dragen_metrics.ploidy_est if self.raw_mets.dragen_metrics else None
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

    def _get_avg_probeqc(self, label, probeqc, bp):
        """
        Retrieve average metric from the probeqc_after dict.  All ProbeQCs are after, unless they
        are being run on UMI-containing assays.  If this is the case, the probeqc_before dict
        captures metrics before UMI deduplication occurs.
        :return:
        """
        try:
            this_cov = self._calc_cov(probeqc, label)
            return self._calc_metric(bp, this_cov)
        except:
            return None

    def _get_cov_uniformity(self, average_depth):
        """
        Calculate the uniformity of coverage
        The percentage of targeted base positions in which the read depth is greater than 0.2 times the mean region
        target coverage depth.
        :return:
        """
        try:
            # calculate using perlocusread instead to get decimal places
            return "{:.1f}".format(float(sum(float(i) > float(average_depth)*0.2 for i in
                                             list(self.gatk_depth_cov_cnts.perlocus.values())) /
                                         len(list(self.gatk_depth_cov_cnts.perlocus.values())) * 100))
        except:
            return None


class Writer:
    def __init__(self, mets):
        self.mets = mets

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
            if metric in self.mets.metrics_req:
                if val:
                    metric_dict = {'metric': str(metric), 'value': str(val)}
                    to_write['sampleRunMetrics'].append(metric_dict)

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

    metrics_req = args.metrics.split(" ")

    raw_mets = RawMetricCollector(args)
    samp_mets = MetricPrep(raw_mets, metrics_req)

    writer = Writer(samp_mets)
    if args.outfile_txt:
        writer.write_to_text(args.outfile_txt)
    if args.outfile:
        writer.write_cgd_new(args.outfile)


if __name__ == "__main__":
    main()
