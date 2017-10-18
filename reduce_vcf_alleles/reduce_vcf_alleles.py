#!/usr/bin/env python

# USAGE: python reduce_vcf_alleles.py <input_vcf> <output_vcf>
# Correct the variant callers that like to break the rules.
# John Letaw

# Can use pyvcf package as well to deal with VCF files.  For this
# script, I found it overkill, and just stuck to vanilla.
import argparse
import itertools

VERSION = '0.1.1'

def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """

    parser = argparse.ArgumentParser(description='VCF\'s created by FreeBayes, among others, do not consistently write REF and ALT alleles to the output file.  This script will reduce these representations down in a consistent manner, and will allow annotations to be applied to variants correctly.')
    parser.add_argument('input_vcf', type=file, help='Input VCF to reduce alleles.')
    parser.add_argument('output_vcf', help='Output VCF to reduce alleles.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    args = parser.parse_args()
    return args

def compare_begin(ref, alt):
    """
    Work from beginning of the string, remove matching bases.
    See perform_begin().
    """
    i = 0
    while ref[0] == alt[0] and (len(ref) > 1) and (len(alt) > 1):
        ref = ref[1:]
        alt = alt[1:]
        i += 1

    return i


def compare_end(ref, alt):
    """
    Same as compare_begin, but this time we will work from the end of the string 
    in towards the beginning.
    """
    i = 0
    while ref[-1] == alt[-1] and (len(ref) > 1) and (len(alt) > 1):
        ref = ref[:-1]
        alt = alt[:-1]
        i -= 1

    return i
    


def perform_begin(ref, alt):
    """
    Decide whether we will reduce ref and alt alleles starting at the beginning
    or the end of the string.
    [TAG, TAT] = beginning (no match at end)
    [CCCTG, CG] = end (end is matching)
    [CCCTG, CCG] = end then beginning (we will always perform a start check after an end check)
    [AG, AGATCG] = end (ambiguity here with second base equal to last base, try to maintain start position)
    [TGGGT, TGGT] = end (will left-align indels with homopolymeric insertions/deletions)
    return True if begin
    return False if end
    """
    
    if len(ref) == 1 or len(alt) == 1:
        return 0
    else:
        # compare_end()
        if ref[-1] == alt[-1]:
            return 1
        # compare_begin()
        elif ref[0] == alt[0]:
            return 2


def main():

    args = supply_args()
    handle_out = open(args.output_vcf, 'w')

    with args.input_vcf as vcf:
        for variant in vcf:
            if variant[0] == '#':
                handle_out.write(variant)
            else:
                variant = variant.rstrip('\n').split('\t')
                chrom = variant[0]
                coord = int(variant[1])
                ref = variant[3]
                alt = variant[4]
                temp_offset = []

                for allele in alt.split(','):
                    if perform_begin(ref, allele) == 1:
                        offset_end = compare_end(ref, allele)
                        offset_begin = compare_begin(ref[:offset_end], allele[:offset_end])
                    elif perform_begin(ref, allele) == 2:
                        offset_begin = compare_begin(ref, allele)
                        offset_end = None
                    else:
                        offset_begin = None
                        offset_end = None

                    temp_offset.append([offset_begin, offset_end])

                temp_offset.sort()
                temp_offset = list(temp_offset for temp_offset,_ in itertools.groupby(temp_offset))

                if len(temp_offset) == 1:
                    ref = ref[offset_begin:offset_end]
                    if offset_begin != None:
                        coord = str(coord + offset_begin)
                    if ',' in alt:
                        alt = ','.join([x[temp_offset[0][0]:temp_offset[0][1]] for x in alt.split(',')])
                    else:
                        alt = alt[offset_begin:offset_end]

                to_write = [chrom, str(coord), variant[2], ref, alt]
                for entry in variant[5:]:
                    to_write.append(entry)

                handle_out.write('\t'.join(to_write))
                handle_out.write('\n')

    handle_out.close()

if __name__ == "__main__":
    main()


