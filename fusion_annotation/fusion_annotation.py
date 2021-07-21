#!/usr/bin/env python
# USAGE: python fusion_annotation.py [-f] <starfusion> <sample_level_metrics.txt> <genome_refpath>

import argparse
import json
import pysam
import re
from collections import OrderedDict
from ensembldb import EnsemblDbImport
from gtf import Gtf

VERSION = '0.5.1'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='Annotate starfusion output with oncotator and nts sequences.')
    parser.add_argument('--starfusion', help='STAR-Fusion output star-fusion.fusion_candidates.final.abridged.FFPM')
    parser.add_argument('--json_sample_metrics', help="Sample level metrics file")
    parser.add_argument('--ensembl_mapping', help="TSV mapping of ENST to RefSeq transcript IDs.")
    parser.add_argument('--path_to_fasta', help='Full path to fasta.fa file, with index file')
    parser.add_argument('--gencode_gtf', help='GENCODE GTF, as used in CTAT resource package.')
    parser.add_argument('--filt', action='store_true', help='filter out specific ribosomal and mitochondrial fusions')
    parser.add_argument('--cr_rem', action='store_true', help='this is a v1.10.0 star-fusion output, get rid of the crs')
    parser.add_argument('--output', default='starfusion_output.bedpe', help='output file')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class StarFusionOutput:
    def __init__(self, filename, cr_rem=False):
        self.filename = filename
        if cr_rem:
            self.header, self.fusions = self.out_create_cr_rm(self.remove_cr())
        else:
            self.header, self.fusions = self.out_create()
        self.cds_left = self._get_tx_list(left=True)
        self.cds_right = self._get_tx_list(left=False)

    def _get_tx_list(self, left=True):
        """
        Get a list of all ENST IDs in output.
        :return:
        """
        tx_list = []
        for fuse in self.fusions:
            if left:
                tx_list.append(fuse.cds_left_id)
            else:
                tx_list.append(fuse.cds_right_id)
        return tx_list

    def remove_cr(self):
        """
        In v1.10.0 ^Ms are being added all over the file.
        :return:
        """
        no_crs = []
        first = True
        with open(self.filename, 'r') as my_sf:
            for line in my_sf:
                if first:
                    temp = line.replace('\r', '').rstrip()
                    first = False
                else:
                    temp += line
                    first = True
                    no_crs.append(temp)
        return no_crs

    def out_create_cr_rm(self, my_sf):
        """
        If we have to use remove_cr, use this function.
        :return:
        """
        fusions = []
        header = []
        for fusion in my_sf:
            if fusion.startswith('#'):
                header = fusion.rstrip('\n').lstrip('#').split('\t')
            else:
                fusions.append(StarFusionOutputLine(fusion, header))
        return header, fusions

    def out_create(self):
        fusions = []
        header = []
        with open(self.filename, 'r') as my_sf:
            for fusion in my_sf:
                if fusion.startswith('#'):
                    header = fusion.rstrip('\n').lstrip('#').split('\t')
                else:
                    fusions.append(StarFusionOutputLine(fusion, header))
        return header, fusions


class StarFusionOutputLine:
    """
    COLUMN HEADERS FROM V1.9.1 STAR-Fusion output, include all extra headings as well.
    #FusionName
    JunctionReadCount
    SpanningFragCount
    LeftGene
    LeftLocalBreakpoint
    LeftBreakpoint
    RightGene
    RightLocalBreakpoint
    RightBreakpoint
    SpliceType
    LargeAnchorSupport
    NumCounterFusionLeft
    NumCounterFusionRight
    FAR_left
    FAR_right
    TrinGG_Fusion
    LeftBreakDinuc
    LeftBreakEntropy
    RightBreakDinuc
    RightBreakEntropy
    FFPM
    annots
    CDS_LEFT_ID
    CDS_LEFT_RANGE
    CDS_RIGHT_ID
    CDS_RIGHT_RANGE
    PROT_FUSION_TYPE
    FUSION_MODEL
    FUSION_CDS
    FUSION_TRANSL
    PFAM_LEFT
    PFAM_RIGHT
    """
    def __init__(self, line, header):
        self.line = line.rstrip('\n').split('\t')
        self.header = header
        self.fusion = self._fusion_dict_create()
        if len(self.line) != len(self.header):
            raise Exception("Data line length does not match header length.")
        self.cds_left_id = self.fusion['CDS_LEFT_ID']
        self.cds_right_id = self.fusion['CDS_RIGHT_ID']

    def _fusion_dict_create(self):
        """
        Create a dictionary that has column headings as keys, data points as values.
        :return:
        """
        fusion = {}
        for entry in self.header:
            fusion[entry] = self.line[self.header.index(entry)]
        return fusion


class FusionAnnot:
    """

    """
    def __init__(self, line):
        self.line = line.fusion
        self.output_line = self._convert_starfusion_to_bedpe()

    @staticmethod
    def _kinase_from_pfam(pfam):
        """
        Determine whether the Pkinase annotation exists within PFAM strings.
        NOTE: Temporarily returning dots, until this functionality has been improved.
        :return:
        """
        return '.'
        # pfam = pfam.rstrip('\n').split('^')
        # for entry in pfam:
        #     if entry.startswith('Pkinase'):
        #         return 'YES'
        # return 'NO'

    def _convert_starfusion_to_bedpe(self):
        """
        Takes starfusion run output as the input and reorders a line from file to bedpe format
        :return:
        """
        annot = OrderedDict()
        valsL = self.line['LeftBreakpoint'].split(':')
        chrL = valsL[0]
        posL = valsL[1]
        strandL = valsL[2]
        valsR = self.line['RightBreakpoint'].split(':')
        chrR = valsR[0]
        posR = valsR[1]
        strandR = valsR[2]
        geneL = self.line['LeftGene'].split('^')[0]
        geneR = self.line['RightGene'].split('^')[0]

        annot['chrom1'] = chrL
        annot['start1'] = str(int(posL) - 1)
        annot['end1'] = posL
        annot['chrom2'] = chrR
        annot['start2'] = str(int(posR) - 1)
        annot['end2'] = posR
        annot['name'] = '--'.join([geneL, geneR])
        annot['score'] = '0'
        annot['strand1'] = strandL
        annot['strand2'] = strandR
        annot['JunctionReadCount'] = self.line['JunctionReadCount']
        annot['SpanningFragCount'] = self.line['SpanningFragCount']
        annot['SpliceType'] = self.line['SpliceType']
        annot['HGVSGene1'] = geneL
        annot['EnsGene1'] = self.line['LeftGene'].split('^')[1]
        annot['HGVSGene2'] = geneR
        annot['EnsGene2'] = self.line['RightGene'].split('^')[1]
        annot['LargeAnchorSupport'] = self.line['LargeAnchorSupport']
        annot['LeftBreakDinuc'] = self.line['LeftBreakDinuc']
        annot['LeftBreakEntropy'] = self.line['LeftBreakEntropy']
        annot['RightBreakDinuc'] = self.line['RightBreakDinuc']
        annot['RightBreakEntropy'] = self.line['RightBreakEntropy']
        if 'PROT_FUSION_TYPE' in self.line:
            annot['PROT_FUSION_TYPE'] = self.line['PROT_FUSION_TYPE']
            annot['FUSION_CDS'] = self.line['FUSION_CDS']
            annot['FAR_left'] = self.line['FAR_left']
            annot['FAR_right'] = self.line['FAR_right']
            annot['PFAM_LEFT'] = self.line['PFAM_LEFT']
            annot['PFAM_RIGHT'] = self.line['PFAM_RIGHT']
            annot['KINASE_IN_PFAM_LEFT'] = self._kinase_from_pfam(annot['PFAM_LEFT'])
            annot['KINASE_IN_PFAM_RIGHT'] = self._kinase_from_pfam(annot['PFAM_RIGHT'])

        # This means J_FFPM and S_FFPM are not listed.  Currently, for compatibility with CGD, placeholders
        # for these values need to be sent.
        if 'J_FFPM' in self.line:
            annot['J_FFPM'] = self.line['J_FFPM']
        else:
            annot['J_FFPM'] = '-1'

        if 'S_FFPM' in self.line:
            annot['S_FFPM'] = self.line['S_FFPM']
        else:
            annot['S_FFPM'] = '-1'

        if 'FFPM' in self.line:
            annot['FFPM'] = self.line['FFPM']

        annot['ENST_LEFT_ID'] = self.line['CDS_LEFT_ID']
        annot['ENST_RIGHT_ID'] = self.line['CDS_RIGHT_ID']
        annot['REFSEQ_LEFT_ID'] = '.'
        annot['REFSEQ_RIGHT_ID'] = '.'
        annot['CCDS_LEFT_ID'] = '.'
        annot['CCDS_RIGHT_ID'] = '.'
        annot['EXON_LEFT'] = '.'
        annot['EXON_RIGHT'] = '.'

        return annot


def calc_on_target(ontarget_count, *args):
    """
    calculate on-target cpm from junction & spanning frag count and sample level metrics
    :return:
    """
    num = 0.0
    for arg in args:
        num += float(arg)
    j_s_cpm = (num / ontarget_count) * 1e6
    return round(j_s_cpm, 3)


def get_sample_met(filename):
    """

    :param filename:
    :return:
    """
    with open(filename, 'r') as samp_metrics_fh:
        metrics = json.load(samp_metrics_fh)
        for entry in metrics['sampleRunMetrics']:
            if entry['metric'] == 'total_on_target_transcripts':
                return float(entry['value'])
    return None


def reverse_complement(seq):
    """
    Takes the reverse complement of genomic sequence. Sequence should only contain nts letters or exception error
    will come up
    Args:
        seq (str): sequence of nts
    Returns:
        seq (str): sequence of nts
    Examples:
        seq = "TCGGGCCC"
        print(reverse_complement(seq))
        > GGGCCCGA
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    return bases


def get_nucleotides_with_samtools(mafline, genome_refpath, fusionside):
    """
    Use samtools faidx to find genomic sequences 25 nts upstream (if maflite_left), downatream (if maflite_right).
    Args:
        mafline(list):  [chrom, start, stop, nt/-, nt/-]
        genome_refpath (str): string containing the fai reference path
        fusionside
    Returns:
        maf_seqs (list of str): sequence of nts for each row of coordinates in the maflite file
    Examples:
        original bash command: samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa 2:112615802-112615805
    """
    interval = None
    if fusionside == 'left':
        if mafline[4] == "-":
            interval = mafline[0] + ':' + mafline[1] + '-' + str(int(mafline[1])+25)
        else:
            interval = mafline[0] + ':' + str(int(mafline[2])-25) + '-' + mafline[2]

    if fusionside == 'right':
        if mafline[4] == "-":
            interval = mafline[0] + ':' + str(int(mafline[2])-25) + '-' + mafline[2]
        else:
            interval = mafline[0] + ':' + mafline[1] + '-' + str(int(mafline[1])+25)

    seq = pysam.FastaFile(filename=genome_refpath).fetch(region=interval)
    return [interval, seq]


def collect_exon_coords(txs, gencode, left):
    """
    Get the exon-coordinate mappings.
    :return:
    """
    exons = {}
    mygtf = Gtf(gencode, txs=txs)
    for entry in txs:
        tx_only = mygtf.filt_tx(mygtf.raw, entry)
        feat_only = mygtf.filt_feature(tx_only)
        exons[entry] = mygtf.exon_to_coord(feat_only, left)
    return exons


def collect_ccds(txs, gencode):
    """
    Get CCDS values from the GENCODE GTF.
    :param txs:
    :param gencode:
    :return:
    """
    ccds = {}
    mygtf = Gtf(gencode, txs=txs)
    for entry in txs:
        tx_only = mygtf.filt_tx(mygtf.raw, entry)
        feat_only = mygtf.filt_feature(tx_only, feat='transcript')
        ccds[entry] = mygtf.exon_to_ccds(feat_only)
    return ccds


def main():
    args = supply_args()
    sf_out = open(args.output, 'w')
    if args.ensembl_mapping:
        enst_refseq = EnsemblDbImport(args.ensembl_mapping).enst_refseq
    else:
        enst_refseq = None

    on_target = get_sample_met(args.json_sample_metrics)
    if args.filt:
        hard_filter = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1",
                       "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2",
                       "TRNA", "TRNC", "TRND", "TRNE", "TRNF", "TRNG", "TRNH",
                       "TRNI", "TRNK", "TRNL1", "TRNL2", "TRNM", "TRNN", "TRNP",
                       "TRNQ", "TRNR", "TRNS1", "TRNS2", "TRNT", "TRNV", "TRNW", "TRNY",
                       "IGH-@-ext", "IGH@", "IGH-@", "IGL-@", "IGL@", "ABC7-42389800N19.1", "IGH@-ext"]
        regex_filt = [r'A[FPLC][0-9]{6}\.[0-9]{1,2}',
                      r'RP[0-9]{1,2}-[0-9]{1,4}[A-Z]{1}[0-9]{1,2}\.[0-9]{1,2}',
                      r'CT[CDAB]{1}-[0-9]{3,4}[A-Z]{1}[0-9]+.[0-9]+',
                      r'XX[ybacos]{0,3}-[0-9A-Z_]+\.[0-9]{1,2}',
                      r'GS1-[0-9]{2,3}[A-Z]{1}[0-9]{1,2}\.[0-9]{1,2}',
                      r'hsa-mir-[0-9]+',
                      r'U[0-9]{5}\.[0-9]{1,2}',
                      r'MT-[A-Z0-9]{1,4}',
                      r'KB-[0-9]{1,4}[A-Z]{1}[0-9]\.[0-9]{1,2}',
                      r'FLJ[0-9]{5}']

    else:
        hard_filter = None
        regex_filt = None

    sf = StarFusionOutput(args.starfusion, args.cr_rem)
    my_sf = sf.fusions

    # Get the transcript list from the star-fusion file, so we can then create exon-coord mappings from gencode.
    if args.gencode_gtf:
        left_exons = collect_exon_coords(sf.cds_left, args.gencode_gtf, left=True)
        right_exons = collect_exon_coords(sf.cds_right, args.gencode_gtf, left=False)
        ccds = collect_ccds(sf.cds_left + sf.cds_right, args.gencode_gtf)
    else:
        left_exons, right_exons, ccds = None, None, None

    cols = False
    # If there are no results, we need at least a header to send to the CGD...
    hard_header = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2",
                   "JunctionReadCount", "SpanningFragCount", "SpliceType", "HGVSGene1", "EnsGene1", "HGVSGene2",
                   "EnsGene2", "LargeAnchorSupport", "LeftBreakDinuc", "LeftBreakEntropy", "RightBreakDinuc",
                   "RightBreakEntropy", "PROT_FUSION_TYPE", "FUSION_CDS", "FAR_left", "FAR_right", "PFAM_LEFT",
                   "PFAM_RIGHT", "KINASE_IN_PFAM_LEFT", "KINASE_IN_PFAM_RIGHT", "J_FFPM", "S_FFPM", "FFPM",
                   "NormalizedFrags", "leftgene", "leftseq", "rightgene", "rightseq", "combinedseq", "ENST_LEFT_ID",
                   "ENST_RIGHT_ID", "REFSEQ_LEFT_ID", "REFSEQ_RIGHT_ID", "CCDS_LEFT_ID", "CCDS_RIGHT_ID",
                   "EXON_LEFT", "EXON_RIGHT"]
    if len(my_sf) == 0:
        sf_out.write('\t'.join(hard_header))
        sf_out.write('\n')

    header_set = False
    for line in my_sf:
        fusion = FusionAnnot(line)
        linebedpe = fusion.output_line

        if not header_set:
            if 'PROT_FUSION_TYPE' not in linebedpe:
                hard_header = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1",
                               "strand2", "JunctionReadCount", "SpanningFragCount", "SpliceType", "HGVSGene1",
                               "EnsGene1", "HGVSGene2", "EnsGene2", "LargeAnchorSupport", "LeftBreakDinuc",
                               "LeftBreakEntropy", "RightBreakDinuc", "RightBreakEntropy", "J_FFPM", "S_FFPM", "FFPM",
                               "NormalizedFrags", "leftgene", "leftseq", "rightgene", "rightseq", "combinedseq",
                               "ENST_LEFT_ID", "ENST_RIGHT_ID", "REFSEQ_LEFT_ID", "REFSEQ_RIGHT_ID", "CCDS_LEFT_ID",
                               "CCDS_RIGHT_ID", "EXON_LEFT", "EXON_RIGHT"]
            header_set = True

        # get seqs left and right
        mafline_left = [linebedpe['chrom1'], linebedpe['start1'], linebedpe['end1'],
                        linebedpe['LeftBreakDinuc'][0], linebedpe['strand1']]
        mafline_right = [linebedpe['chrom2'], linebedpe['start2'], linebedpe['end2'],
                         linebedpe['RightBreakDinuc'][0], linebedpe['strand2']]
        seq_left = get_nucleotides_with_samtools(mafline_left, args.path_to_fasta, "left")
        seq_right = get_nucleotides_with_samtools(mafline_right, args.path_to_fasta, "right")

        if on_target:
            linebedpe['NormalizedFrags'] = calc_on_target(on_target, linebedpe['JunctionReadCount'],
                                                          linebedpe['SpanningFragCount'])
        else:
            linebedpe['NormalizedFrags'] = '-1'

        linebedpe['leftgene'] = seq_left[0]
        linebedpe['leftseq'] = seq_left[1]
        linebedpe['rightgene'] = seq_right[0]
        linebedpe['rightseq'] = seq_right[1]

        if linebedpe['strand1'] == "-":
            combined_left = reverse_complement(linebedpe['leftseq'])
        else:
            combined_left = linebedpe['leftseq']
        if linebedpe['strand2'] == "-":
            combined_right = reverse_complement(linebedpe['rightseq'])
        else:
            combined_right = linebedpe['rightseq']
        combinedseq = combined_left + combined_right
        linebedpe['combinedseq'] = combinedseq

        if enst_refseq:
            enst_l = linebedpe['ENST_LEFT_ID']
            enst_r = linebedpe['ENST_RIGHT_ID']
            if enst_l in enst_refseq:
                linebedpe['REFSEQ_LEFT_ID'] = enst_refseq[enst_l]
            if enst_r in enst_refseq:
                linebedpe['REFSEQ_RIGHT_ID'] = enst_refseq[enst_r]

        if left_exons:
            enst_l = linebedpe['ENST_LEFT_ID']
            coord = linebedpe['end1']
            if enst_l in left_exons:
                if coord in left_exons[enst_l]:
                    linebedpe['EXON_LEFT'] = left_exons[enst_l][coord]

        if right_exons:
            enst_r = linebedpe['ENST_RIGHT_ID']
            coord = linebedpe['end2']
            if enst_r in right_exons:
                if coord in right_exons[enst_r]:
                    linebedpe['EXON_RIGHT'] = right_exons[enst_r][coord]

        if ccds:
            enst_l = linebedpe['ENST_LEFT_ID']
            enst_r = linebedpe['ENST_RIGHT_ID']
            if enst_l in ccds:
                if ccds[enst_l]:
                    linebedpe['CCDS_LEFT_ID'] = ccds[enst_l]
            if enst_r in ccds:
                if ccds[enst_r]:
                    linebedpe['CCDS_RIGHT_ID'] = ccds[enst_r]

        # Writing section.
        filter_me = False
        if not cols:
            sfoutcolumns = [str(x) for x in linebedpe]
            if len(sfoutcolumns) == 0:
                sfoutcolumns = hard_header
            sf_out.write('\t'.join(sfoutcolumns))
            sf_out.write('\n')
            cols = True
        if args.filt:
            gene1 = linebedpe['HGVSGene1']
            gene2 = linebedpe['HGVSGene2']
            if (gene1 in hard_filter) or (gene2 in hard_filter):
                filter_me = True
            for filt in regex_filt:
                if re.match(filt, gene1) or re.match(filt, gene2):
                    filter_me = True
        if not filter_me:
            to_write = [str(x) for x in linebedpe.values()]
            sf_out.write('\t'.join(to_write))
            sf_out.write('\n')


if __name__ == "__main__":
    main()
