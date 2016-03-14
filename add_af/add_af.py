#!/usr/bin/env python

### Add an allele frequency entry, useful for Pindel VCF's.
### USAGE: python add_af.py <input VCF> <output VCF>

import sys

def createHeaderEntry():
    """
    This goes in the header of the VCF, and looks like this:
    ##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth, how many reads support this allele">
    Modify to allow for different Number, Type, and Description field to be passed.
    """

    af_header = "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Variant allele frequency\">\n"
    return af_header


def findIndex(splitting, delim, to_find):
    """
    Find the index of the field in the FORMAT column you care about.  Usually will be AD.
    """

    for entry in splitting.split(delim):
        if entry == to_find:
            curr_index = splitting.split(delim).index(to_find)

    return curr_index


def calcAF(ref, alt):
    """
    Calculate the allele frequency and format the output.
    """
    
    if ref + alt != 0:
        af = alt / (ref + alt + 0.0)
        if af != 0.0:
            af = "{:.7f}".format(alt / (ref + alt + 0.0)).rstrip('0')
    else:
        af = 0.0

    return af


def main():
    
    handle_vcf = open(sys.argv[1], 'rU')
#    handle_vcf = open("/home/exacloud/clinical/Galaxy/database/files/019/dataset_19126.dat", 'rU')
    handle_out = open(sys.argv[2], 'w')
    curr_index = None
    delim = ':'  # We will usually target ':' delimeters, as that is what the VCF uses.
    to_find = 'AD'  # Generally should be looking for AD's (allele depth) here.

    with handle_vcf as myvcf:
        for variant in myvcf:
            if variant[0] != "#":

                split_variant = variant.rstrip('\n').split('\t')
                format = split_variant[8]
                normal = split_variant[9]
                tumor = split_variant[10]
                normal_ref = int(normal.split(':')[1].split(',')[0])
                normal_alt = int(normal.split(':')[1].split(',')[1])
                tumor_ref = int(tumor.split(':')[1].split(',')[0])
                tumor_alt = int(tumor.split(':')[1].split(',')[1])

                if curr_index == None:
                    curr_index = findIndex(format, delim, to_find)

                normal_af = calcAF(normal_ref, normal_alt)
                tumor_af = calcAF(tumor_ref, tumor_alt)

                new_format = format + delim + "AF"
                new_normal = normal + delim + str(normal_af)
                new_tumor = tumor + delim + str(tumor_af)
                handle_out.write('\t'.join(['\t'.join(split_variant[:8]), new_format, new_normal, new_tumor, '\n']))
                
            else:
                if "ID=AD," in variant:
                    handle_out.write(variant)
                    handle_out.write(createHeaderEntry())
                else:
                    handle_out.write(variant)

main()
