"""
This script is specifically written to further annotate data coming from the GATK Germline CNV workflow.

# Iniitialize DB
#db = gffutils.create_db('test/ref_GRCh37.p13_top_level.gff3', 'ref_GRCh37.p13_top_level.db',
merge_strategy="create_unique")
"""

import argparse
from collections import OrderedDict
import gffutils
import logging
import os
import re
import vcfpy

logger = logging.getLogger(__name__)

chrom_map = {'1': 'NC_000001.10', '2': 'NC_000002.11', '3': 'NC_000003.11', '4': 'NC_000004.11', '5': 'NC_000005.9',
             '6': 'NC_000006.11', '7': 'NC_000007.13', '8': 'NC_000008.10', '9': 'NC_000009.11', '10': 'NC_000010.10',
             '11': 'NC_000011.9', '12': 'NC_000012.11', '13': 'NC_000013.10', '14': 'NC_000014.8', '15': 'NC_000015.9',
             '16': 'NC_000016.9', '17': 'NC_000017.10', '18': 'NC_000018.9', '19': 'NC_000019.9', '20': 'NC_000020.10',
             '21': 'NC_000021.8', '22': 'NC_000022.10', 'X': 'NC_000023.10', 'Y': 'NC_000024.9', 'MT': 'NC_012920.1'}

VERSION = '0.0.1'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', help='Input VCF, from PostProcessGermlineCNVCalls.')
    parser.add_argument('gffutils_db', help='GFFutils library sqlite3 db.')
    parser.add_argument('gatk_counts', help='Counts file produced by GATK CollectReadCounts.')
    parser.add_argument('segments_dir', help='Directory of segments files, contained VCFs with CNV '
                                             'calls for each historical sample.')
    parser.add_argument('outfile', help='Output TSV file, for CGD import.')
    parser.add_argument('outfile_nogene', help='Output TSV file for entries without genes.')
    parser.add_argument('logfile', help='Output log file.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class VcfFile:
    """
    Provide dictionary (ordered) of input CNV VCFs:
    {(chrom, start, stop): VCF Record}
    """
    def __init__(self, vcf_path: str):
        self.vcf_path = vcf_path

    def parse(self):
        vrnts = OrderedDict()
        reader = vcfpy.Reader.from_path(self.vcf_path)
        for vrnt in reader:
            uniq_key = (vrnt.CHROM, vrnt.POS, vrnt.INFO['END'])
            vrnts[uniq_key] = vrnt
        return vrnts


class VcfRegion:
    """
    Connect genomic regions with gene names, size of region.
    """
    def __init__(self, chrom: str, start: int, end: int, dbfile: str):
        self.db = gffutils.FeatureDB(dbfile)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.size_bp = self.end - self.start + 1
        self.genes = self._get_genes()

    def _get_genes(self):
        """
        Get the list of genes from gffutils db associated with a particular genomic region.
        :return:
        """
        genes = []
        for entry in self.db.region(seqid=chrom_map[self.chrom], start=self.start, end=self.end):
            if entry.featuretype == 'gene':
                this_gene = entry.attributes['gene']
                genes.extend(this_gene)
        return list(set(genes))


class Writer:
    def __init__(self, outfile, outfile_nogene, counts):
        self.outfile = outfile
        self.outfile_nogene = outfile_nogene
        self.counts = counts

    def write(self, myvcf, stats, pop_freqs, vcf_count, gffutils_db):
        write_nogenes = open(self.outfile_nogene, 'w')
        with open(self.outfile, 'w') as f:
            f.write('\t'.join(self._write_header()))
            f.write('\n')
            write_nogenes.write('\t'.join(self._write_header()))
            write_nogenes.write('\n')
            for vrnt in myvcf:
                chrom = vrnt[0]
                start = vrnt[1]
                end = vrnt[2]
                region = VcfRegion(chrom=chrom, start=start, end=end, dbfile=gffutils_db)
                size_bp = str(region.size_bp)
                genes = ', '.join(region.genes)
                geno = str(myvcf[vrnt].calls[0].data['GT'])
                geno_lbl = self._assign_geno_label(geno)
                copy_number = str(myvcf[vrnt].calls[0].data['CN'])
                probes = str(myvcf[vrnt].calls[0].data['NP'])
                avg_count = str(int(stats[vrnt]['total'] / stats[vrnt]['regions']))
                min_count = str(stats[vrnt]['min'])
                max_count = str(stats[vrnt]['max'])
                qc_all = str(myvcf[vrnt].calls[0].data['QA'])
                qc_one = str(myvcf[vrnt].calls[0].data['QS'])
                qc_start = str(myvcf[vrnt].calls[0].data['QSS'])
                qc_end = str(myvcf[vrnt].calls[0].data['QSE'])
                try:
                    local_freq = "%.3f" % float((pop_freqs[vrnt]['1'] + pop_freqs[vrnt]['2']) / vcf_count)
                    del_freq = "%.3f" % float(pop_freqs[vrnt]['1'] / vcf_count)
                    dup_freq = "%.3f" % float(pop_freqs[vrnt]['2'] / vcf_count)
                except KeyError:
                    local_freq = '0.000'
                    del_freq = '0.000'
                    dup_freq = '0.000'
                    logger.info(f"{vrnt} not found in segments directory.")
                if geno != '0':
                    if not genes:
                        to_write = [chrom, str(start), str(end), genes, geno_lbl, copy_number, probes, avg_count,
                                    max_count, min_count, size_bp, local_freq, del_freq, dup_freq, qc_all, qc_one,
                                    qc_start, qc_end]
                        write_nogenes.write('\t'.join(to_write))
                        write_nogenes.write('\n')
                    else:
                        to_write = [chrom, str(start), str(end), genes, geno_lbl, copy_number, probes, avg_count,
                                    max_count, min_count, size_bp, local_freq, del_freq, dup_freq, qc_all, qc_one,
                                    qc_start, qc_end]
                        f.write('\t'.join(to_write))
                        f.write('\n')
        write_nogenes.close()
    @staticmethod
    def _assign_geno_label(geno):
        """
        Assign labels based on called genotype.
        0 - None
        1 - Del
        2 - Dup
        :return:
        """
        geno_label = {'0': 'Dip',
                      '1': 'Del',
                      '2': 'Dup'}
        return geno_label[geno]

    @staticmethod
    def _write_header():
        header = ['REGION_CHROM', 'REGION_START', 'REGION_STOP', 'GENE', 'CALL', 'COPY_NUMBER', 'PROBES',
                  'AVG_READ_COUNT', 'MAX_READ_COUNT', 'MIN_READ_COUNT', 'SIZE', 'LOCAL_FREQ', 'DEL_FREQ', 'DUP_FREQ',
                  'QC_ALL', 'QC_ONE', 'QC_START', 'QC_END']
        return header


class GatkCounts:
    """
    Collect all of the counts in to dictionary format.
    {<CHROM>: {(<START>, <END>): <COUNT>}}
    """
    def __init__(self, filename: str):
        self.filename = filename
        self.counts = self._gather_counts()

    def _gather_counts(self):
        counts = {}
        with open(self.filename) as f:
            for entry in f:
                if not entry.startswith('@') and not entry.startswith('CONTIG'):
                    entry = entry.rstrip('\n').split('\t')
                    chrom = entry[0]
                    start = entry[1]
                    end = entry[2]
                    count = entry[3]
                    if chrom not in counts:
                        counts[chrom] = OrderedDict()
                    counts[chrom][(start, end)] = count
        return counts


def get_stats_region(counts, myvcf):
    """
    Given coordinates, collect depth information across the region.
    :param counts:
    :param myvcf:
    :return:
    """
    stats = {}
    for vrnt in myvcf:
        chrom = vrnt[0]
        start = vrnt[1]
        end = vrnt[2]
        stat = get_counts_match(counts, chrom, start, end)
        stats[vrnt] = stat
    return stats


def get_counts_match(counts, chrom, start, end):
    stat = {}
    for entry in counts[chrom]:
        counts_start = int(entry[0])
        counts_end = int(entry[1])
        total = int(counts[chrom][entry])
        min_cnt = int(counts[chrom][entry])
        max_cnt = int(counts[chrom][entry])
        if start == counts_start:
            stat = {'total': total,
                    'regions': 1,
                    'min': min_cnt,
                    'max': max_cnt}
            if end == counts_end:
                return stat
        elif end == counts_end:
            stat['total'] += total
            stat['regions'] += 1
            stat['min'] = min(stat['min'], min_cnt)
            stat['max'] = max(stat['max'], max_cnt)
            return stat
        else:
            if 'total' in stat:
                stat['total'] += total
                stat['regions'] += 1
                stat['min'] = min(stat['min'], min_cnt)
                stat['max'] = max(stat['max'], max_cnt)


class VcfCollect:
    """
    From a directory structure, locate VCFs.
    """
    def __init__(self, seg_dir):
        self.seg_dir = seg_dir
        self.big_list = self._get_vcfs()
        self.total_vcfs = len(self.big_list)

    def _get_vcfs(self):
        """
        Get the list of all VCFs that we need to process.
        :return:
        """
        big_list = []
        vcf_dir = os.listdir(self.seg_dir)
        for d in vcf_dir:
            this_dir = self.seg_dir + d
            if self.find_runid_history(d):
                vcf_list = [this_dir + '/' + x for x in os.listdir(this_dir) if x.endswith('.vcf')]
                big_list.extend(vcf_list)
        return big_list

    @staticmethod
    def find_runid_history(runid):
        """
        Look through all histories, then find those that match a given runid.
        :param runid:
        :return:
        """
        my_regex = "([0-9]{6}_N[BS][0-9]{6}_[0-9]{4}_[A-Z0-9]{10})_*"
        my_match = re.findall(my_regex, runid)
        return my_match


def collect_cnv_counts(vcf_list: list):
    """
    Go through the list of VCF files and collect number of del/dup counts assocaited with cnv regions.
    :param vcf_list:
    :return:
    """
    all_cnvs = {}
    for entry in vcf_list:
        this_vcf = VcfFile(entry).parse()
        for vrnt in this_vcf:
            gt = this_vcf[vrnt].calls[0].data['GT']
            if vrnt not in all_cnvs:
                all_cnvs[vrnt] = {'0': 0,
                                  '1': 0,
                                  '2': 0}
            all_cnvs[vrnt][str(gt)] += 1
    return all_cnvs


def main():
    args = supply_args()
    logging.basicConfig(filename=args.logfile, level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    logger.info('Begin annotating.')
    logger.info("Collecting VCFs from segments directory and establishing local cohort of calls.")
    vcfs = VcfCollect(args.segments_dir)
    vcf_list = vcfs.big_list
    vcf_count = vcfs.total_vcfs
    all_cnvs = collect_cnv_counts(vcf_list)
    logger.info("Parsing input VCF.")
    myvcf = VcfFile(args.input_vcf)
    myvcf_calls = myvcf.parse()
    logger.info("Parsing counts file.")
    counts = GatkCounts(args.gatk_counts)
    my_stats = get_stats_region(counts.counts, myvcf_calls)
    logger.info("Writing results.")
    Writer(args.outfile, args.outfile_nogene, counts).write(myvcf_calls, my_stats, all_cnvs, vcf_count, args.gffutils_db)
    logger.info('End annotating.')

if __name__ == '__main__':
    main()
