#!/usr/bin/env python
# USAGE: python getSeq_revcomplement.py <maflite> <genome_refpath> <seq_out_filename>
# Args:
#   maflite
#       chr	start	end	ref_allele	alt_allele
#       2	89521179	89521180	T	-
#       2	89521179	89521180	T	-
#   genome_refpath= '/home/users/patterja/BioCoders/DataResources/AnnotationSources/Oncotator/oncotator_v1_ds_April052016'
#   seq_out_filename = <name of output file>
# Returns:
#   list of sequence strings as txt file

import sys
import subprocess

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
    genome_refpath = sys.argv[2]
    maflite=[]
    with open(sys.argv[1], 'r') as maf_file:
        maf_file.next()
        for line in maf_file:
            vals = line.strip().split()
            maflite.append(vals)
    #hg19
    #genome_refpath = '/home/users/patterja/BioCoders/DataResources/Genomes/hg19/release-75/genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
    #see fusion_annotation.py


    seq_list = get_nucleotides_with_samtools(maflite, genome_refpath)


    with open(seq_out_filename, 'w') as outfile:
        for i in seq_list:
            outfile.writelines('\t'.join(i))
            outfile.write('\n')


if __name__ == "__main__":
    main()
