#!/usr/bin/env python
# USAGE: python fusion_annotation.py [-f] <starfusion> <sample_level_metrics.txt> <genome_refpath>

import argparse
import json
import re
from collections import OrderedDict

import pysam

VERSION = '0.2.6'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='Annotate starfusion output with oncotator and nts sequences.')
    parser.add_argument('--starfusion', help='STAR-Fusion output star-fusion.fusion_candidates.final.abridged.FFPM')
    parser.add_argument('--json_sample_metrics', help="Sample level metrics file")
    parser.add_argument('--path_to_fasta', help='Full path to fasta.fa file, with index file')
    parser.add_argument('--filt', action='store_true', help='filter out specific ribosomal and mitochondrial fusions')
    parser.add_argument('--output', default='starfusion_output.bedpe',
                        help='Full path to fasta.fa file, with index file')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class StarFusionOutput:
    def __init__(self, filename):
        self.filename = filename
        self.header, self.fusions = self.out_create()

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
    COLUMN HEADERS FROM V1.8.1 STAR-Fusion output, include all extra headings as well.
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

    def _kinase_from_pfam(self, pfam):
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


def main():
    args = supply_args()
    sf_out = open(args.output, 'w')
    on_target = get_sample_met(args.json_sample_metrics)
    if args.filt:
        hard_filter = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1",
                       "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2",
                       "TRNA", "TRNC", "TRND", "TRNE", "TRNF", "TRNG", "TRNH",
                       "TRNI", "TRNK", "TRNL1", "TRNL2", "TRNM", "TRNN", "TRNP",
                       "TRNQ", "TRNR", "TRNS1", "TRNS2", "TRNT", "TRNV", "TRNW", "TRNY",
                       "IGH-@-ext", "IGH@-ext", "IGH@", "IGL-@"]
        regex_filt = [r'A[LC][0-9]{6}\.[0-9]{1}',
                      r'RP[0-9]{1,2}-[0-9]{1,4}[A-Z]{1}[0-9]{1,2}\.[0-9]{1,2}',
                      r'CT[CDAB]{1}-[0-9]{3,4}[A-Z]{1}[0-9]+.[0-9]+',
                      r'XX[ybacos]{0,3}-[0-9A-Z_]+\.[0-9]{1,2}',
                      r'GS1-[0-9]{2,3}[A-Z]{1}[0-9]{1,2}\.[0-9]{1,2}',
                      r'hsa-mir-[0-9]+']

    else:
        hard_filter = None
        regex_filt = None

    my_sf = StarFusionOutput(args.starfusion).fusions
    cols = False
    # If there are no results, we need at least a header to send to the CGD...
    hard_header = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2",
                   "JunctionReadCount", "SpanningFragCount", "SpliceType", "HGVSGene1", "EnsGene1", "HGVSGene2",
                   "EnsGene2", "LargeAnchorSupport", "LeftBreakDinuc", "LeftBreakEntropy", "RightBreakDinuc",
                   "RightBreakEntropy", "PROT_FUSION_TYPE", "FUSION_CDS", "FAR_left", "FAR_right", "PFAM_LEFT",
                   "PFAM_RIGHT", "KINASE_IN_PFAM_LEFT", "KINASE_IN_PFAM_RIGHT", "J_FFPM", "S_FFPM", "FFPM",
                   "NormalizedFrags", "leftgene", "leftseq", "rightgene", "rightseq", "combinedseq"]
    if len(my_sf) == 0:
        sf_out.write('\t'.join(hard_header))
        sf_out.write('\n')

    header_set = False
    for line in my_sf:
        fusion = FusionAnnot(line)
        linebedpe = fusion.output_line

        if not header_set:
            if 'PROT_FUSION_TYPE' in linebedpe:
                hard_header = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1",
                               "strand2", "JunctionReadCount", "SpanningFragCount", "SpliceType", "HGVSGene1", "EnsGene1",
                               "HGVSGene2", "EnsGene2", "LargeAnchorSupport", "LeftBreakDinuc", "LeftBreakEntropy",
                               "RightBreakDinuc", "RightBreakEntropy", "PROT_FUSION_TYPE", "FUSION_CDS", "FAR_left", "FAR_right",
                               "PFAM_LEFT", "PFAM_RIGHT", "KINASE_IN_PFAM_LEFT", "KINASE_IN_PFAM_RIGHT", "J_FFPM", "S_FFPM", "FFPM",
                               "NormalizedFrags", "leftgene", "leftseq", "rightgene", "rightseq", "combinedseq"]
            else:
                hard_header = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1",
                               "strand2", "JunctionReadCount", "SpanningFragCount", "SpliceType", "HGVSGene1", "EnsGene1",
                               "HGVSGene2", "EnsGene2", "LargeAnchorSupport", "LeftBreakDinuc", "LeftBreakEntropy",
                               "RightBreakDinuc", "RightBreakEntropy", "J_FFPM", "S_FFPM", "FFPM",
                               "NormalizedFrags", "leftgene", "leftseq", "rightgene", "rightseq", "combinedseq"]
            header_set = True

        # get seqs left and right
        mafline_left = [linebedpe['chrom1'], linebedpe['start1'], linebedpe['end1'],
                        linebedpe['LeftBreakDinuc'][0], linebedpe['strand1']]
        mafline_right = [linebedpe['chrom2'], linebedpe['start2'], linebedpe['end2'],
                         linebedpe['RightBreakDinuc'][0], linebedpe['strand2']]
        seq_left = get_nucleotides_with_samtools(mafline_left, args.path_to_fasta, "left")
        seq_right = get_nucleotides_with_samtools(mafline_right, args.path_to_fasta, "right")

        if on_target:
            linebedpe['NormalizedFrags'] = calc_on_target(on_target, linebedpe['JunctionReadCount'], linebedpe['SpanningFragCount'])
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
