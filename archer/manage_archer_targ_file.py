#!/usr/bin/env python

from archer import ArcherOps
from cgd import CgdOps, CgdPropertiesParser
from samplesheet_class import SampleSheetReader
from vrnts import ArcherVrnts, ArcherWriteVrnts
import argparse


VERSION = '0.1.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('archer_creds', help='Archer API token.')
    parser.add_argument('cgd_props', help='CGD Client properties file.')
    parser.add_argument('samplesheet', help='Illumina SampleSheet.')
    parser.add_argument('runid', help='Run ID for corresponding to the SampleSheet you are parsing.')
    parser.add_argument('--outloc', default='.', help='Output VCF write location')
    parser.add_argument('--archer_url', default='https://ohsukdldev.analysis.archerdx.com/rest_api',
                        help='Archer URL Base.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def collect_reported(barcodes, runid, cgd):
    """
    Place all reported variants from CGD in same structure, removing duplicates.
    :return:
    """
    all_report_vars = []
    for barcode in barcodes:
        terms = [runid, barcode]
        reported = cgd.get_url(terms)
        for report in reported:
            rrep = reformat_reported(report)
            if rrep not in all_report_vars:
                all_report_vars.append(rrep)
    return all_report_vars


def reformat_reported(vrnt):
    """
    Take something that looks like this:
    {'genomicVariantId': 2218832, 'chromosome': 'chr1', 'positionStart': 216363636,
    'referenceBase': 'A', 'variantBase': 'G'}
    Return this:
    (chromosome, positionStart, referenceBase, variantBase)
    :return:
    """
    return vrnt['chromosome'], str(vrnt['positionStart']), vrnt['referenceBase'], vrnt['variantBase']


def archer_version_incr(name):
    """
    If we need to send a new file, increment the version number.
    File name something like this:
    Archer Comprehensive Targets v1.5
    Take last split piece and increment.
    :param name:
    :return:
    """
    if '.' not in name:
        raise Exception("Archer targets file malformed.")
    base = name.split('.')[0]
    vers = int(name.split('.')[-1])
    new_vers = str(vers + 1)
    return  '.'.join([base, new_vers])

def archer_name_reformat(name):
    """
    Now, make the file look like this:
    Archer_Comprehensive_Targets-v1.5.vcf
    :param name:
    :return:
    """
    name = name.split(' ')
    return f"{name[0]}_{name[1]}_{name[2]}-{name[3]}.vcf"

def write_string(vrnt):
    """
    Produce the write string to be appended to the end of the Archer targeted mutations file.
    Goes from this:
    ('chr10', '112356246', 'CAGA', 'C')
    Should look like this:
    chr7	148543637	EZH2_p.N57fs*25	G	GT	.	.	Archer_Gene=EZH2
    Though, we will not be including info in col3, and Archer_Gene will be set to NA
    :param vrnt:
    :return:
    """
    new_vrnt = [vrnt[0], vrnt[1], '.', vrnt[2], vrnt[3], '.', '.', 'Archer_Gene=.']
    return '\t'.join(new_vrnt) + '\n'


def outfile_loc(base, outfile):
    """
    Take the location as defined in the argument and create the outfile write location.
    :return:
    """
    return '/'.join([base.rstrip('/'), outfile])


def main():
    args = supply_args()
    # Set send to False, until we have found new variants to add.
    send = False
    # Get the oauth consumer_secret from properties file.
    props = CgdPropertiesParser(args.cgd_props)
    cgd_token = props.prop_vals['oauthsecret']
    # Get the archer token from file.
    archer_props = CgdPropertiesParser(args.archer_creds)
    archer_token = archer_props.prop_vals['token']
    # Load the SampleSheet
    ss = SampleSheetReader(args.samplesheet)
    barcodes = ss.archer_barcodes
    # Find the Archer VCF we are interested in.
    archer = ArcherOps(token=archer_token, url=args.archer_url)
    archer_vcf, archer_name = archer.get_targ_vcf()
    archer_vrnts = ArcherVrnts(archer_vcf).vrnts
    # Get reported variants from CGD
    cgd = CgdOps(cgd_token)
    reported = collect_reported(barcodes, args.runid, cgd)
    for entry in reported:
        if entry not in archer_vrnts:
            # We have a new variant for the file, will need to send the new target file to Archer.
            send = True
            archer_vcf += write_string(entry)
    # If we have new variants, prepare the new file.
    if send:
        new_targets = archer_version_incr(archer_name)
        targets_name_reformat = archer_name_reformat(new_targets)
        outfile = outfile_loc(args.outloc, targets_name_reformat)
        writer = ArcherWriteVrnts(outfile)
        writer.write_archer(archer_vcf)
        resp = archer.put_dep_file(outfile)
        print(resp)

if __name__ == "__main__":
    main()
