#!/usr/bin/env python

import sys

def writeFakeDoc(chrom, coord):
    return (chrom + ':' + coord + '\t100\t100.00\t100\n')

def main():

    handle_doc = open(sys.argv[1], 'w')

    handle_doc.write("Locus\tTotal_Depth\tAverage_Depth_sample\tDepth_for_normal\n")    

    handle_doc.write(writeFakeDoc("7", "6026775"))
    handle_doc.write(writeFakeDoc("13", "32929387"))
    handle_doc.write(writeFakeDoc("1", "169519049"))
    handle_doc.write(writeFakeDoc("14", "75513883"))

main()
