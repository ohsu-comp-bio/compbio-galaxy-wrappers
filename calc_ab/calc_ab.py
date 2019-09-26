#!/usr/bin/env python

from file_types import vcfreader
import argparse
import json

VERSION = '0.1.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--infile',  help='Input VCF')
    parser.add_argument('--outfile', help='Output JSON')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def write_json_out(metric, outfile):
    """
    Prepare output json file.
    :return:
    """
    outfile = open(outfile, 'w')
    out_metric = {'allele_balance': metric}
    json.dump(out_metric, outfile)
    outfile.close()


def calc_af(ad):
    """
    Given a comma sep list of ref and alt depths, return allele freq.
    :param ad:
    :return:
    """
    ref = int(ad.split(',')[0])
    alt = int(ad.split(',')[1])

    try:
        return alt / (alt + ref + 0.0)
    except ZeroDivisionError:
        return None

def af_list(samps):
    """
    Given a VcfSamples object, return a list of VAFs.
    :return:
    """
    afs = []
    for samp in samps:
        if 'AF' in samp:
            afs.append(float(samp['AF']))
        elif 'AD' in samp:
            afs.append(calc_af(samp['AD']))
        else:
            raise Exception('FORMAT must currently contain either an AD or AF field to capture VAF.')
    return afs


def create_af_counts(sample_list, my_vcf, min_vaf = 0.35, max_vaf = 0.65):
    """
    Create a dictionary {SAMPLE: [BAD VAF, TOTAL]
    :return:
    """
    af_dict = {}
    for vrnt in my_vcf.myvcf.values():

        i = 0
        for samp in vrnt.samples:
            if sample_list[i] not in af_dict:
                af_dict[sample_list[i]] = {'TOTAL': 0,
                                           'BAD_VAF': 0}
            i += 1

        if len(vrnt.ALT) == 1 and len(vrnt.REF) == 1:
            if len(vrnt.ALT[0]) == 1 and len(vrnt.REF[0]) == 1:
                if 'DP' in vrnt.INFO:
                    if int(vrnt.INFO['DP']) >= 50:
                        i = 0
                        for samp in vrnt.samples:
                            if samp['GT'] == '0/1':
                                af_dict[sample_list[i]]['TOTAL'] += 1
                                if af_list(vrnt.samples)[i] < min_vaf or af_list(vrnt.samples)[i] > max_vaf:
                                    af_dict[sample_list[i]]['BAD_VAF'] += 1
                            i += 1

    return af_dict

def main():
    args = supply_args()
    my_vcf = vcfreader.VcfReader(args.infile)
    sample_list = my_vcf.header.sample_list
    af_dict = create_af_counts(sample_list, my_vcf)
    if len(sample_list) == 1:
        try:
            bad_vaf = (af_dict[sample_list[0]]['BAD_VAF'] / (af_dict[sample_list[0]]['TOTAL'] + 0.0)) * 100
        except ZeroDivisionError:
            bad_vaf = '100'
        write_json_out(bad_vaf, args.outfile)
    else:
        print('Multi sample functionality not completely implemented...')
        for entry in af_dict:
            try:
                bad_vaf = (af_dict[entry]['BAD_VAF'] / (af_dict[entry]['TOTAL'] + 0.0)) * 100
            except ZeroDivisionError:
                bad_vaf = '100'
            print('\t'.join([entry, str(bad_vaf)]))

if __name__ == "__main__":
    main()
