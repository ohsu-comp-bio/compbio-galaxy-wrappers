#!/usr/bin/env python

# DESCRIPTION: Given a panel of normals VCF, as created by
# panel_of_normals.py, annotate or remove entries in a sample VCF.
# USAGE: use_mutect2_calls.py <pon> <infile> <outfile>
# CODED BY: John Letaw

from __future__ import print_function
import argparse

VERSION = '0.3.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('mutect2', help='MuTect2 VCF.')
    parser.add_argument('infile', help='Input VCF.')
    parser.add_argument('hotspots', help='Hotspots VCF.')
    parser.add_argument('outfile', help='Output VCF.')
    parser.add_argument('outfile_m2', help='Output VCF M2 calls.')
    parser.add_argument('outfile_indels', help='Output VCF with Mutect2 '
                                               'indels.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def parse_m2(args, mutect_out):
    """
    Run through the Mutect2 VCF and grab what we need.
    :return:
    """
    annot = {}
    indels = []
    filtered = True
    check = True
    with open(args.mutect2, 'rU') as pon:
        for line in pon:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                coord = line[0]
                pos = line[1]
                ref = line[3]
                alt = line[4]
                format = line[7]
                tlod = find_info_val(format, 'TLOD')
                annot[(coord, pos, ref, alt)] = tlod

                if float(tlod) > 40.0 and (len(ref) > 1 or len(alt) > 1):
                    indels.append(line)
            else:
                mutect_out.write(line)
                if check:
                    mutect_out.write("##reference=/opt/installed/galaxy_genomes/hg19/Homo_sapiens_assembly19.fasta")
                    check = False
                if line.startswith("##FILTER") and filtered:
                    mutect_out.write("##FILTER=<ID=m2,Description=\"Variant found in MuTect2.\">\n")
                    filtered = False

    return annot, indels


def find_info_val(splitting, to_find, delim=';'):
    """
    Find the value of the field in the INFO column you care about.  Usually
    will be AD.
    """
    for entry in splitting.split(delim):
        if to_find in entry:
            try:
                value = entry.split('=')[1]
                return value
            except:
                raise Exception(to_find + " is in the field, but is not "
                                          "formatted as expected.")


def replace_filter(filt, anno):
    """
    Remove the dot if it's there, otherwise just append something with a
    semicolon.
    :return:
    """
    if filt == '.' or filt == 'PASS':
        return anno
    else:
        return ';'.join([filt, anno])


def write_new_vcf(args, annot, indels, outfile_indels):
    """
    Write a new VCF with the m2 tag in FILTER section.
    :return:
    """
    outfile = open(args.outfile, 'w')
    just_wrote = []
    filtered = True
    with open(args.infile, 'rU') as infile:
        for line in infile:
            if not line.startswith('#'):
                sline = line.rstrip('\n').split('\t')
                uniq_key = (sline[0], sline[1], sline[3], sline[4])
                af = float(sline[7].split(';')[0].split('=')[1])
                if len(sline[3]) == 1 and len(sline[4]) == 1:
                    if uniq_key in annot:
                        sline[6] = replace_filter(sline[6], 'm2')
                        tlod_str = 'TLOD=' + annot[uniq_key]
                        sline[7] = ';'.join([sline[7], tlod_str])
                        outfile.write('\t'.join(sline))
                        outfile.write('\n')
                        just_wrote.append(uniq_key)
                    # No, make this an argument.
                    elif af >= 0.2:
                        outfile.write(line)
                        just_wrote.append(uniq_key)
                    elif af >= 0.08 and (len(sline[3]) > 1 or len(sline[4]) > 1):
                        outfile.write(line)
                        just_wrote.append(uniq_key)
                        print(line)
                    elif 'hotspot' in sline[6]:
                        outfile.write(line)
                        just_wrote.append(uniq_key)
                    else:
                        pass
    #                    outfile.write(line)
            else:
                outfile.write(line)
                if line.startswith("##FILTER") and filtered:
                    outfile.write("##FILTER=<ID=m2,Description=\"Variant found in MuTect2.\">\n")
                    outfile.write("##INFO=<ID=TLOD,Number=A,Type=Float,Description=\"Tumor LOD score\">\n")
#                    outfile.write("##FILTER=<ID=m2_hotspot,Description=\"Variant would have been filtered but appears in the hotspot list.\">\n")
                    filtered = False

    for entry in indels:
        if (entry[0], entry[1], entry[3], entry[4]) not in just_wrote:
            entry[6] = 'm2'
            outfile_indels.write('\t'.join(entry))
            outfile_indels.write('\n')

    outfile.close()
    outfile_indels.close()

    return just_wrote


def hotspots_grab(filename):
    """

    :param filename:
    :return:
    """
    hotspots = []
    with open(filename, 'rU') as hots:
        for line in hots:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                chrom = line[0]
                pos = line[1]
                ref = line[3]
                alt = line[4]
                uniq_key = (chrom, pos, ref, alt)
                if uniq_key not in hotspots:
                    hotspots.append(uniq_key)
    return hotspots


def write_m2_vcf(filename, outfile, just_wrote, hotspots):
    """
    """
    handle_out = open(outfile, 'w')
    check = True
    filtered = True
    with open(filename, 'rU') as infile:
        for line in infile:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                uniq_key = (line[0], line[1], line[3], line[4])
                if uniq_key in hotspots and uniq_key not in just_wrote:
                    pass
                    # line[6] = replace_filter(line[6], 'm2_hotspot')
                    # handle_out.write('\t'.join(line))
                    # handle_out.write('\n')
            else:
                handle_out.write(line)
                if check:
                    handle_out.write("##reference=/opt/installed/galaxy_genomes/hg19/Homo_sapiens_assembly19.fasta")
                    check = False
                elif line.startswith("##FILTER") and filtered:
                    handle_out.write("##FILTER=<ID=m2,Description=\"Variant found in MuTect2.\">\n")
                    # handle_out.write("##FILTER=<ID=m2_hotspot,Description=\"Variant would have been filtered but appears in the hotspot list.\">\n")
                    filtered = False

    handle_out.close()


def main():

    args = supply_args()
    hotspots = hotspots_grab(args.hotspots)
    mutect_out = open(args.outfile_indels, 'w')
    annot, indels = parse_m2(args, mutect_out)
    just_wrote = write_new_vcf(args, annot, indels, mutect_out)
    write_m2_vcf(args.mutect2, args.outfile_m2, just_wrote, hotspots)

if __name__ == "__main__":
    main()
