#!/usr/bin/env python

# DESCRIPTION: Remove calls when forward and reverse read counts do not make
#  sense.
# USAGE:
# CODED BY: John Letaw

from __future__ import print_function
import argparse

VERSION = '0.4.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('infile', help='Input VCF.')
    parser.add_argument('cosmic', help='COSMIC VCF.')
    parser.add_argument('hotspots', help='Hotspots VCF.')
    parser.add_argument('outfile', help='Output VCF.')
    parser.add_argument('outfile_removed', help='Output VCF of removed '
                                                'sites.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def adj_alts(saf, sar, srf, srr, thresh=0.015):
    """
    Adjust alternate read strand counts to reflect ref proportion.
    For 0/0 and 0/1 only.
    """
    sref_af = (srf*1.0/srr)
    salt_af = (saf*1.0/sar)

    if sref_af > 1:
        if salt_af <= 1:
            new_saf = saf
            new_sar = saf/sref_af 
        elif salt_af > 1 and salt_af <= sref_af:
            new_saf = saf
            new_sar = saf/sref_af 
        elif salt_af > 1 and salt_af > sref_af:
            new_saf = sar*sref_af
            new_sar = sar
        else:
            raise Exception("This shouldn't be possible.")

    elif sref_af < 1:
        if salt_af >= 1:
            new_saf = sar*sref_af
            new_sar = sar
        elif salt_af < 1 and salt_af >= sref_af:
            new_saf = sar*sref_af
            new_sar = sar
        elif salt_af < 1 and salt_af < sref_af:
            new_saf = saf
            new_sar = saf/sref_af 
        else:
            print('\t'.join([str(saf), str(sar), str(srf), str(srr)]))
            raise Exception("This shouldn't be possible.")

    else:
        if salt_af <= 1:
            new_saf = saf
            new_sar = saf/sref_af 
        else:
            new_saf = sar*sref_af
            new_sar = sar
            
        
    new_af = (new_saf + new_sar) / (new_saf + new_sar + srf + srr)

    if new_af >= thresh:
        return True
    else:
        return False


def cosmic_grab(filename):
    """

    :return:
    """
    cosmic_vars = {}
    with open(filename, 'rU') as cosmic:
        for line in cosmic:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                chrom = line[0]
                pos = line[1]
                cosm = line[2]
                ref = line[3]
                alt = line[4]
                if cosm not in cosmic_vars:
                    cosmic_vars[(chrom, pos, ref, alt)] = cosm
    return cosmic_vars


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


def main():

    args = supply_args()

    handle_vcf = open(args.infile, 'rU')
    handle_out = open(args.outfile, 'w')
    handle_outr = open(args.outfile_removed, 'w')

    hotspots = hotspots_grab(args.hotspots)
    cosmic_vars = cosmic_grab(args.cosmic)
    check = True

    with handle_vcf as vcf:
        for line in vcf:
            if not line.startswith('#'):
                line = line.rstrip('\n').split('\t')
                chrom = line[0]
                pos = line[1]
                ref = line[3]
                alt = line[4]
                info = line[7]
                geno = line[9].split(':')[0]
                saf = info.split(';')[34].split('=')[1]
                sar = info.split(';')[36].split('=')[1]
                srf = info.split(';')[37].split('=')[1]
                srr = info.split(';')[39].split('=')[1]
                all_bal = info.split(';')[0].split('=')[1]

                if ',' not in saf and ',' not in sar:
                    saf = int(saf)
                    sar = int(sar)
                    srf = int(srf)
                    srr = int(srr)

                    if geno != '1/1':
                        if saf > 10 and sar > 10 and srf > 0 and srr > 0:
#                            if not ((saf*1.0/sar) > 1.2 and (srf*1.0/srr) < 0.8): 
#                                if not ((saf*1.0/sar) < 0.8 and (srf*1.0/srr) > 1.2): 
                            if all_bal =='0' and adj_alts(saf, sar, srf, srr):
                                handle_out.write('\t'.join(line))
                                handle_out.write('\n')
                            elif all_bal != '0':
                                handle_out.write('\t'.join(line))
                                handle_out.write('\n')
                        elif (chrom, pos, ref, alt) in hotspots:
                            line[6] = replace_filter(line[6], 'strand_hotspot')
                            handle_out.write('\t'.join(line))
                            handle_out.write('\n')
                        elif (chrom, pos, ref, alt) in cosmic_vars:
                            line[2] = cosmic_vars[(chrom, pos, ref, alt)]
                            handle_outr.write('\t'.join(line))
                            handle_outr.write('\n')
                        else:
                            pass

                    else:
                        handle_out.write('\t'.join(line))
                        handle_out.write('\n')
            else:
                handle_out.write(line)
                handle_outr.write(line)
                if line.startswith("##FILTER") and check:
                    handle_out.write("##FILTER=<ID=strand_hotspot,Description=\"Variant would be filtered but is seen in the hotspots file.\">\n")
                    check = False

    handle_out.close()
    handle_outr.close()

main()
