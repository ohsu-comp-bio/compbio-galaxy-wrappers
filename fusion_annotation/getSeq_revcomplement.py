#!/usr/bin/env python
# USAGE: python getSeq_revcomplement.py <starfusionOutput> <genome_refpath> <left> <seq_out_prefix>
# Args:
#   Input: 
#       1) starfusionOutput 
#       2) genome fasta file
#       3) "left" or "right"
#       4) outfile name
#   intermediate file is maflite:
#   maflite
#       chr	start	end	ref_allele	alt_allele
#       2	89521179	89521180	T	-
#       2	89521179	89521180	T	-
#   genome_refpath= '/home/exacloud/lustre1/BioCoders/DataResources/Genomes/hg19/broad_variant/genome/Homo_sapiens_assembly19.fasta'
#   seq_out_prefix = <prefix of outfile>
# Returns:
#   list of sequence strings as txt file

import sys
import argparse
import subprocess

VERSION="0.2.0"
#April 12, 2018

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='Input starfusion results and retrieve the 25 nts sequences.')
    parser.add_argument('starfusionOutput', help='STAR-Fusion output star-fusion.fusion_candidates.final.abridged.FFPM')
    parser.add_argument('genome_refpath', help='Full path to fasta.fa file, with index file')
    parser.add_argument('fusion_side', help="'left' or 'right'")
    parser.add_argument('seq_out_prefix', help="outfile name")
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args



def convert_starfusion_to_maflites(starfusion_handle, fusion_side):
    """
    Convert starfusion results into maflite format. Simplest input for oncotator and std input for seq retriever. Two lists, left & right breakpoint
    Args:
        starfusionOutput (file): opened for read file type
	fusion_side (flag): "left" or "right"
    :Returns:
        maflite (list of list)
    """
    maflite = []
    if fusion_side == "left":
        for line in starfusion_handle:
            if not line.startswith('#'):
                    vals = line.strip().split()
                    valsL = vals[5].split(':')
                    valsR = vals[7].split(':')
                    chrL = valsL[0]
                    posL = valsL[1]
                    dinucL = vals[9][0]
                    chrR = valsR[0]
                    posR = valsR[1]
                    dinucR = vals[11][0]
                    maflite.append([chrL, str(int(posL)-1), posL, dinucL, '-'])
    elif fusion_side == "right":
        for line in starfusion_handle:
            if not line.startswith('#'):
                    vals = line.strip().split()
                    valsL = vals[5].split(':')
                    valsR = vals[7].split(':')
                    chrL = valsL[0]
                    posL = valsL[1]
                    dinucL = vals[9][0]
                    chrR = valsR[0]
                    posR = valsR[1]
                    dinucR = vals[11][0]
                    maflite.append([chrR, str(int(posR)-1), posR, dinucR, '-'])
    return maflite


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

def get_nucleotides_with_samtools(maflite, genome_refpath):
    """
    Use samtools faidx to find genomic sequences 25 nts upstream (if maflite_left), downatream (if maflite_right).
    Args:
        maflite (list of list): simple maflite file [chrom, start, stop, nt/-, nt/-]
        genome_refpath (str): string containing the fai reference path
    Returns:
        maf_seqs (list of str): sequence of nts for each row of coordinates in the maflite file
    Examples:
        original bash command: samtools faidx Homo_sapiens.GRCh37.75.dna.primary_assembly.fa 2:112615802-112615805
    """
    sequence_list=[]
    for item in maflite:
        interval = item[0] + ':' + item[1] + '-' + str(int(item[1])+25)
        cmd = ['samtools', 'faidx', genome_refpath, interval]
        seq = subprocess.check_output(cmd)
        print seq
        seq_strip=seq.strip().split()
        if item[4] == '-':
            print "reverse complementing"
            revseq = reverse_complement(seq_strip[1])
            sequence_list.append([seq_strip[0],revseq])
        else:
            sequence_list.append(seq_strip)
    return sequence_list


def main():
    args = supply_args()

    with open(args.starfusionOutput, 'r') as starfusion_handle:
	   maf = convert_starfusion_to_maflites(starfusion_handle, args.fusion_side)
       seq_list = get_nucleotides_with_samtools(maf, args.genome_refpath)

    #hg19
    #genome_refpath = '/home/users/patterja/BioCoders/DataResources/Genomes/hg19/release-75/genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
    #see fusion_annotation.py
    if args.fusion_side =="left":
        seq_out_filename = args.seq_out_prefix + "_left.txt"
    elif args.fusion_side =="right":
        seq_out_filename = args.seq_out_prefix + "_right.txt"

    with open(seq_out_filename, 'w') as outfile:
        for i in seq_list:
            outfile.writelines('\t'.join(i))
            outfile.write('\n')


if __name__ == "__main__":
    main()
