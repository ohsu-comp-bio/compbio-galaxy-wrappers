#!/usr/bin/env python
# USAGE: python fusion_annotation.py [-f] <starfusion> <sample_level_metrics.txt> <genome_refpath>


import subprocess
import argparse
import re
import json

VERSION = '0.1.3'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='Annotate starfusion output with oncotator and nts sequences.')
    parser.add_argument('-f', action='store_true', help='filter out specific ribosomal and mitochondrial fusions')
    parser.add_argument('starfusion', help='STAR-Fusion output star-fusion.fusion_candidates.final.abridged.FFPM')
    parser.add_argument('json_sample_metrics', help="Sample level metrics file")
    parser.add_argument('path_to_fasta', help='Full path to fasta.fa file, with index file')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

def convert_starfusion_to_bedpe(line):
    """
    Takes starfusion run output as the input and reorders a line from file to bedpe format
    Args:
        line (string): unconcatenated or stripped. Just parsing line by line
    Returns:
        linebedpe (list): item per column
    Examples:
    """
    vals = line.strip().split()
    valsL = vals[5].split(':')
    valsR = vals[7].split(':')
    chrL = valsL[0]
    posL = valsL[1]
    geneL = vals[4].split('^')[0]
    strandL = valsL[2]
    chrR = valsR[0]
    posR = valsR[1]
    geneR = vals[6].split('^')[0]
    strandR = valsR[2]
    linebedpe=[chrL, str(int(posL) - 1), posL, chrR, str(int(posR) - 1), posR, '--'.join([geneL, geneR]), '0', strandL, strandR, vals[1], vals[2], vals[3], geneL, vals[4].split('^')[1], geneR, vals[4].split('^')[1]]
    linebedpe.extend(vals[8:15])
    return linebedpe

def reverse_complement(seq):
    """
    Takes the reverse complement of genomic sequence. Sequence should only contain nts letters or exception error will come up
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
    Returns:
        maf_seqs (list of str): sequence of nts for each row of coordinates in the maflite file
    Examples:
        original bash command: samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa 2:112615802-112615805
    """
    #sequence_list=[]
    #for item in mafline:
    if fusionside=='left':
        if mafline[4]=="-":
            interval = mafline[0] + ':' + mafline[1] + '-' + str(int(mafline[1])+25)
        else:
            interval = mafline[0] + ':' + str(int(mafline[2])-25) +'-' + mafline[2]

    if fusionside=='right':
        if mafline[4]=="-":
            interval = mafline[0] + ':' + str(int(mafline[2])-25) +'-' + mafline[2]
        else:
            interval = mafline[0] + ':' + mafline[1] + '-' + str(int(mafline[1])+25)

    cmd = ['samtools', 'faidx', genome_refpath, interval]
    seq = subprocess.check_output(cmd)
    seq_strip = seq.strip().split()
    return seq_strip

def main():
    args = supply_args()
    bedpe = []
    #took this file parsing out of the convert_starfusion_to_bedpe function because of regex matching to filter
    with open(args.json_sample_metrics, 'r') as samp_metrics_fh:
        metrics = json.load(samp_metrics_fh)
        for entry in metrics['sampleRunMetrics']:
            if entry['metric'] == 'total_on_target_transcripts':
                ontarget_count = float(entry['value'])

    if args.f:
        # filterout = ""
        hard_filter = ["ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", 
                   "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "RNR1", "RNR2", 
                   "TRNA", "TRNC", "TRND", "TRNE", "TRNF", "TRNG", "TRNH", 
                   "TRNI", "TRNK", "TRNL1", "TRNL2", "TRNM", "TRNN", "TRNP", 
                   "TRNQ", "TRNR", "TRNS1", "TRNS2", "TRNT", "TRNV", "TRNW", "TRNY"]
        with open(args.starfusion, 'r') as sf_handle:
            for line in sf_handle:
                if not line.startswith('#'):
                    # Removing this statement since filterout is empty.
                    # not re.search(filterout, line)
                    linebedpe = convert_starfusion_to_bedpe(line)

                    #calculate on-target cpm from junction & spanning frag count and sample level metrics
                    # JL: Lab has requested junction and spanning be summed for this calculation.
                    j_s_cpm = ((float(linebedpe[10]) + float(linebedpe[11])) / ontarget_count) * 1e6

                    #get seqs left and right
                    mafline_left = [linebedpe[0], linebedpe[1], linebedpe[2], linebedpe[18][0], linebedpe[8]]
                    mafline_right = [linebedpe[3], linebedpe[4], linebedpe[5], linebedpe[20][0], linebedpe[9]]
                    seq_left = get_nucleotides_with_samtools(mafline_left, args.path_to_fasta, "left")
                    seq_right = get_nucleotides_with_samtools(mafline_right, args.path_to_fasta, "right")

                    linebedpe.extend([round(j_s_cpm, 3)])
                    linebedpe.extend(seq_left)
                    linebedpe.extend(seq_right)
                    if linebedpe[8] == "-":
                        combined_left = reverse_complement(seq_left[1])
                    else:
                        combined_left = seq_left[1]
                    if linebedpe[9] == "-":
                        combined_right = reverse_complement(seq_right[1])
                    else:
                        combined_right = seq_right[1]
                    combinedseq = combined_left + combined_right
                    linebedpe.append(combinedseq)
                    if args.f:
                        if (linebedpe[13] not in hard_filter) and (linebedpe[15] not in hard_filter): 
                            bedpe.append(linebedpe)
                            
     
    sfoutcolumns = ["chrom1","start1","end1","chrom2","start2","end2","name","score","strand1","strand2","JunctionReadCount","SpanningFragCount","SpliceType","HGVSGene1","EnsGene1","HGVSGene2","EnsGene2","LargeAnchorSupport","LeftBreakDinuc","LeftBreakEntropy","RightBreakDinuc","RightBreakEntropy","J_FFPM","S_FFPM", "NormalizedFrags", "leftgene","leftseq","rightgene","rightseq", "combinedseq"]
    
    with open("starfusion_output.bedpe", 'w') as sf_out:
        sf_out.writelines('\t'.join(sfoutcolumns))
        sf_out.write('\n')
        for item in bedpe:
            sf_out.writelines('\t'.join(str(v) for v in item))
            sf_out.write('\n')


if __name__ == "__main__":
    main()
