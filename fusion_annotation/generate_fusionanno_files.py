#!/usr/bin/env python
# USAGE: python fusion_annotation.py <starfusionOutput> <genome_refpath>


import sys
import subprocess
import os
import argparse
from getSeq_revcomplement import *

VERSION = '0.4.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='Annotate starfusion output with oncotator and nts sequences.')
    parser.add_argument('starfusionOutput', help='STAR-Fusion output star-fusion.fusion_candidates.final.abridged.FFPM')
    parser.add_argument('genome_refPath', help='Full path to fasta.fa file, with index file')
    #parser.add_argument('oncotator_dbdir', help='Path to Oncotator database directory')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

def convert_starfusion_to_bedpe(starfusionOutput):
    #with open(starfusionOutput, 'r') as starfusion_ffpm:
    bedpe =[]
    for line in starfusionOutput:
        if not line.startswith('#'):
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
            bedpe.append(linebedpe)
    return bedpe

def main():
    args = supply_args()
    with open(args.starfusionOutput, 'r') as starfusion_handle:
        bedpe=convert_starfusion_to_bedpe(starfusion_handle)

    maf_left=[]
    maf_right=[]
    for i in range(len(bedpe)):
        l=[bedpe[i][0], bedpe[i][1], bedpe[i][2], bedpe[i][18][0], bedpe[i][8]]
        r=[bedpe[i][3], bedpe[i][4], bedpe[i][5], bedpe[i][20][0], bedpe[i][9]]
        maf_left.append(l)
        maf_right.append(r)

    
    seq_left = get_nucleotides_with_samtools(maf_left, args.genome_refPath)
    seq_right = get_nucleotides_with_samtools(maf_right, args.genome_refPath)
    lr_seq_list = map(list.__add__, seq_left, seq_right)

    bedpe_wseq = map(list.__add__, bedpe , lr_seq_list)
    sfoutcolumns = ["chrom1","start1","end1","chrom2","start2","end2","name","score","strand1","strand2","JunctionReadCount","SpanningFragCount","SpliceType","HGVSGene1","EnsGene1","HGVSGene2","EnsGene2","LargeAnchorSupport","LeftBreakDinuc","LeftBreakEntropy","RightBreakDinuc","RightBreakEntropy","J_FFPM","S_FFPM","leftgene","leftseq","rightgene","rightseq"]

    with open("starfusion_output.bedpe", 'w') as sf_out:
        sf_out.writelines('\t'.join(sfoutcolumns))
        sf_out.write('\n')
        for item in bedpe_wseq:
            sf_out.writelines('\t'.join(item))
            sf_out.write('\n')


if __name__ == "__main__":
    main()
