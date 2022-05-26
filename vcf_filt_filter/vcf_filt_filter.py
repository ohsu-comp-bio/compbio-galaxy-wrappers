#!/usr/bin/env python

"""
Perform basic filtering on a VCF based on annotations in the FILTER column.
"""
import argparse
import vcfpy

VERSION = '0.0.2'


def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_vcf', help='Input VCF')
    parser.add_argument('output_vcf', help='Output VCF')
    parser.add_argument('output_vcf_filt', help='Output VCF Filtered Variants')
    parser.add_argument('--callers', help='FILTER annotations corresponding to variant callers being used.')
    parser.add_argument('--sing_rm', help='Remove variants containing only one of these terms '
                                                     'in the FILTER column.')
    parser.add_argument('--sing_vaf', type=float, help='AF to use when determining whether to drop a sing_rm record.')
    parser.add_argument('--only_sing', help='Only include this FILTER label if it is the only variant caller label.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def check_sing_filt(sing_rm, filt, vaf):
    """
    Return True if the sing_rm term is the same as the filt term.
    :param sing_rm:
    :param filt:
    :return:
    """
    for s in sing_rm:
        if [s] == filt:
            if vaf:
                return True
    return False


def check_only_sing(only_sing, filt):
    """
    Check to see whether a variant caller label specified in only_sing is listed by itself.
    :param only_sing:
    :param filt:
    :return:
    """
    if sorted(only_sing.split(' ')) == sorted(filt):
        return filt
    for f in only_sing.split(' '):
        if f in filt and len(filt) > 1:
            filt.remove(f)
    return filt


def main():
    args = supply_args()
    reader = vcfpy.Reader.from_path(args.input_vcf)
    header = reader.header
    writer = vcfpy.Writer.from_path(args.output_vcf, header=header)
    writer_bad = vcfpy.Writer.from_path(args.output_vcf_filt, header=header)
    for vrnt in reader:
        filt = vrnt.FILTER
        vaf = float(vrnt.calls[0].data['AF'])
        low_vaf = False
        if args.sing_vaf:
            low_vaf = (vaf < args.sing_vaf)
        for f in filt:
            if f not in args.callers.split(' '):
                filt.remove(f)
        res = check_sing_filt(args.sing_rm.split(' '), filt, low_vaf)
        if args.only_sing:
            vrnt.FILTER = check_only_sing(args.only_sing, filt)
        if not res:
            writer.write_record(vrnt)
        else:
            writer_bad.write_record(vrnt)
    writer.close()
    writer_bad.close()


if __name__ == "__main__":
    main()
