#!/usr/bin/env python
# USAGE: python fusion_annotation.py <starfusionOutput> <genome_refpath>


import sys
import subprocess
import os
import argparse
from getSeq_revcomplement import *
from starfusion_oncotator import *

VERSION = '0.3.0'

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

def main():
    args = supply_args()


    mafleft_filename = "left_brkpnt.maf"
    mafright_filename = "right_brkpnt.maf"

    with open(args.starfusionOutput, 'r') as starfusion_ffpm:
        maflite_leftls, maflite_rightls = convert_starfusion_to_maflites(starfusion_ffpm)

    # Oncotator Input
    with open(mafleft_filename, 'w') as outfile_mafleft:
        write_mafs(maflite_leftls, outfile_mafleft)

    with open(mafright_filename, 'w') as outfile_mafright:
        write_mafs(maflite_rightls, outfile_mafright)

    # Run Oncotator, output is automatic

    #runOncotator(args.oncotator_dbdir, mafleft_filename, "oncotated_left_output")
    #runOncotator(args.oncotator_dbdir, mafright_filename, "oncotated_right_output")

    # oncotator_dbdir = '.../BioCoders/DataResources/AnnotationSources/Oncotator/oncotator_v1_ds_April052016'
    # genome_refpath = '.../BioCoders/DataResources/Genomes/hg19/release-75/genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'


    seq_list_left = get_nucleotides_with_samtools(maflite_leftls, args.genome_refPath)
    seq_list_right = get_nucleotides_with_samtools(maflite_rightls, args.genome_refPath)
    lr_seq_list = map(list.__add__, seq_list_left, seq_list_right)


    with open(args.starfusionOutput, 'r') as starfusion_ffpm:
        bedpe=convert_starfusion_to_bedpe(starfusion_ffpm)

    bedpe_wseq = map(list.__add__, bedpe , lr_seq_list)
    sfoutcolumns = ["chrom1","start1","end1","chrom2","start2","end2","name","score","strand1","strand2","JunctionReadCount","SpanningFragCount","SpliceType","EnsGene1","EnsGene2","LargeAnchorSupport","LeftBreakDinuc","LeftBreakEntropy","RightBreakDinuc","RightBreakEntropy","J_FFPM","S_FFPM","leftcoord","leftseq","rightgene","rightseq"]

    with open("starfusion_output.bedpe", 'w') as sf_out:
        sf_out.writelines('\t'.join(sfoutcolumns))
        sf_out.write('\n')
        for item in bedpe_wseq:
            sf_out.writelines('\t'.join(item))
            sf_out.write('\n')


if __name__ == "__main__":
    main()
