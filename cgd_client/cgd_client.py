#!/usr/bin/env python

# Galaxy wrapper for cgd_client.jar.
# JAVA8_PATH and CGD_CLIENT_CONFIG must be defined in the Galaxy contrib/ohsu_exacloud_env.sh file.
# USAGE: send_to_cgd.py -h
# CODED BY: John Letaw

from __future__ import print_function
from snp_profile import SnpProfile, CompareProfiles

import argparse
import json
import os
import sys
import shutil

# Including this because is it quite useful...
# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

VERSION = '1.2.6.0'


def supply_args():
    """                                                                                                                            Populate args,
    https://docs.python.org/2.7/library/argparse.html                                                                              """
    parser = argparse.ArgumentParser(description='Galaxy wrapper for cgd_client.jar.')

    parser.add_argument('stdout_log', help='')
    parser.add_argument('endpoint', help='')
    parser.add_argument('java8_path', help='Java 8 PATH as defined in JAVA8_PATH.')
    parser.add_argument('--report_vcf', help='Output VCF if utilizing '
                                             'reportvariants endpoint.')
    parser.add_argument('--report_bed', help='Output BED if utilizing '
                                             'reportvariants endpoint.')
    parser.add_argument('--json_out', help='JSON will be written to this file.', default='cgd_profile_to_send.json')
    parser.add_argument("--pipeline_out", help='')
    parser.add_argument("--cgd_url", help='')
    parser.add_argument("--runid", help='')
    parser.add_argument("--barcodeid", help='')
    parser.add_argument("--qcversion", help='')
    parser.add_argument("--cnvcalls", help='')
    parser.add_argument("--cnvpdf", help='')
    parser.add_argument("--cgd_client", help="Location of the cgd_client.")
    parser.add_argument("--cgd_config", help="Location of the cgd_client config file.")
    parser.add_argument("--tissue",
                        help="Which tissue is the client receiving data for.  This is most applicable to tumor/normal workflows.")
    parser.add_argument("--servicebase",
                        help="The service host name and port + service base. e.g. kdlwebprod02:8080/cgd")

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    args = parser.parse_args()
    return args


def eprint(*args, **kwargs):
    """
    http://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)


def rename_fastqc_output(runid, barcodeid, endpoint, ext):
    """
    CGD needs the filename to be restructured.
    Applies to FastQC and CNV PDF only.
    """
    ext = '.' + ext

    if endpoint == "uploadqcsheet":
        newfile = "/tmp/" + '_'.join([runid, barcodeid, "R1"]) + ext
    elif endpoint == "uploadqcsheetrtwo":
        newfile = "/tmp/" + '_'.join([runid, barcodeid, "R2"]) + ext
    elif endpoint == "cnvpdf":
        newfile = "/tmp/" + '_'.join([runid, barcodeid]) + ext
    else:
        print("Not sending FastQC.")
        return None

    return newfile


def split_url(url, n):
    """
    Split the manual URL and take n elements.
    Not currently in use, what was this even for?
    """

    return url.split('/')[-n:]


def run_cmd(cmd):
    """
    Run the command via subprocess.
    """
    print("Running the following command: ")
    print('\t'.join(cmd))

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if stderr:
        raise Exception(stderr)
    elif type(json.loads(stdout)) is list:
        pass
    elif 'errors' in json.loads(stdout):
        if json.loads(stdout)['errors']:
            raise Exception(json.loads(stdout)['errors'])
    elif 'message' in json.loads(stdout):
        if json.loads(stdout)['message'] != 'ok':
            raise Exception(json.loads(stdout)['message'])
    return stdout


def build_cmd(args, recvd_prof=False):
    """
    Build the command that will send data to the CGD.
    """
    cmd = [args.java8_path, '-jar', args.cgd_client, "-n", args.endpoint, "-c", args.cgd_config]
    newfile = ""

    if args.servicebase:
        cmd.extend(["-s", args.servicebase])
    if args.endpoint == "uploadqcsheet" or args.endpoint == "uploadqcsheetrtwo":
        # For FastQC files, we will create a new file name.
        newfile = rename_fastqc_output(args.runid, args.barcodeid, args.endpoint, 'html')
        print("Copying to " + newfile)
        shutil.copyfile(args.pipeline_out, newfile)
        cmd.extend(["-f", newfile])
    elif args.endpoint == "cnvpdf":
        newfile = rename_fastqc_output(args.runid, args.barcodeid, args.endpoint, 'pdf')
        print("Copying to " + newfile)
        shutil.copyfile(args.pipeline_out, newfile)
        cmd.extend(["-f", newfile])
    elif args.endpoint == "annotationcomplete" or args.endpoint == "annotate" or args.endpoint == "reportedvariants" \
            or (args.endpoint == "snpProfile" and not recvd_prof):
        # The cmd for this endpoint is already set, don't do anything.
        pass
    elif args.endpoint == "updatesamplerun":
        cmd.extend(["-j", args.pipeline_out])
    elif args.endpoint == "snpProfile" and recvd_prof:
        cmd.extend(["-j", args.json_out])
    elif args.endpoint == "none":
        cmd = ["java", "-jar", args.cgd_client, "-f", args.pipeline_out, "-u", args.cgd_url]
    else:
        cmd.extend(["-f", args.pipeline_out])

    if args.runid:
        cmd.append("-r")
        cmd.append(args.runid)
    if args.barcodeid:
        cmd.append("-b")
        cmd.append(args.barcodeid)
    if args.qcversion:
        cmd.append("-v")
        cmd.append(args.qcversion)

    # cmd.append("-d")

    print("The client is receiving data for a " + args.tissue + " sample.")

    return cmd, newfile


def write_vcf_header(outfile):
    """
    Write the VCF header.
    ##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
    :return:
    """
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def prepare_reported(outfile, regions, stdout):
    """

    :return:
    """
    empty = '.'
    for entry in json.loads(stdout):
        if json.loads(stdout):
            chrom = entry['chromosome'][3:]
            pos = entry['positionStart']
            ref = entry['referenceBase']
            alt = entry['variantBase']
            outfile.write('\t'.join([chrom, str(pos), empty, ref, alt, empty, empty, empty]))
            outfile.write('\n')
            start = pos - 1
            regions.write('\t'.join([chrom, str(start), str(pos)]))
            regions.write('\n')

    if not json.loads(stdout):
        outfile.write('\t'.join(['1', '3', empty, 'T', 'C', empty, empty, empty]))
        outfile.write('\n')
        regions.write('\t'.join(['1', '1', '2']))
        regions.write('\n')

    outfile.close()
    regions.close()


def main():
    args = supply_args()

    # Build the command.
    cmd, newfile = build_cmd(args, False)
    # Run the command.
    stdout = run_cmd(cmd)

    if args.endpoint == 'reportedvariants':
        vcf = open(args.report_vcf, 'w')
        regions = open(args.report_bed, 'w')
        write_vcf_header(vcf)
        prepare_reported(vcf, regions, stdout)

        # Get the current profile ready.
        # We run one command here, then another down below.

    if args.endpoint == 'snpProfile':
        json_to_send = SnpProfile(args.pipeline_out).geno_items
        cgd_json = json.loads(stdout)['items']
        compare_snps = CompareProfiles(cgd_json, json_to_send)
        with open(args.json_out, 'w') as to_cgd:
            json.dump(compare_snps.for_cgd, to_cgd)
            #to_cgd.write(compare_snps.for_cgd)
        cmd, newfile = build_cmd(args, True)
        stdout = run_cmd(cmd)

        # Command has been run to retrieve the profile.  Now, we need to figure out if the overlap of the two profiles match.
        # If the match, no need to run anything else.
        # If they don't match, error.
        # If they match, but there are loci that aren't contained within the CGD profile, send the updated loci to CGD.
        # If there is no profile in the CGD, send the created profile.
        # stdout = {"message": "ok", "items": [{"chromosome": "1", "position": 2488153, "genotype": -1},
        #                                      {"chromosome": "1", "position": 16256007, "genotype": 1},
        #                                      {"chromosome": "1", "position": 36937059, "genotype": 2},
        #                                      {"chromosome": "1", "position": 65321388, "genotype": 2},
        #                                      {"chromosome": "1", "position": 215848587, "genotype": 0},
        #                                      {"chromosome": "1", "position": 215960167, "genotype": 2},
        #                                      {"chromosome": "1", "position": 216219781, "genotype": 2},
        #                                      {"chromosome": "1", "position": 237814783, "genotype": 2},
        #                                      {"chromosome": "11", "position": 1267325, "genotype": 0},
        #                                      {"chromosome": "11", "position": 1267917, "genotype": 1},
        #                                      {"chromosome": "11", "position": 1271321, "genotype": 0},
        #                                      {"chromosome": "11", "position": 69462910, "genotype": 1},
        #                                      {"chromosome": "11", "position": 72946020, "genotype": 1},
        #                                      {"chromosome": "11", "position": 108043988, "genotype": 1},
        #                                      {"chromosome": "12", "position": 49444545, "genotype": 1},
        #                                      {"chromosome": "13", "position": 28624294, "genotype": 2},
        #                                      {"chromosome": "13", "position": 73349359, "genotype": 0},
        #                                      {"chromosome": "16", "position": 81929488, "genotype": 1},
        #                                      {"chromosome": "16", "position": 81971403, "genotype": 0},
        #                                      {"chromosome": "17", "position": 11511457, "genotype": 2},
        #                                      {"chromosome": "17", "position": 29508775, "genotype": 0},
        #                                      {"chromosome": "17", "position": 29553485, "genotype": 2},
        #                                      {"chromosome": "17", "position": 37922259, "genotype": 1},
        #                                      {"chromosome": "17", "position": 41244000, "genotype": 0},
        #                                      {"chromosome": "17", "position": 41245466, "genotype": 0},
        #                                      {"chromosome": "17", "position": 62007498, "genotype": 2},
        #                                      {"chromosome": "17", "position": 73552185, "genotype": 2},
        #                                      {"chromosome": "18", "position": 52895531, "genotype": 1},
        #                                      {"chromosome": "18", "position": 60985879, "genotype": 0},
        #                                      {"chromosome": "19", "position": 10267077, "genotype": 1},
        #                                      {"chromosome": "19", "position": 22154732, "genotype": 2},
        #                                      {"chromosome": "19", "position": 22156992, "genotype": 2},
        #                                      {"chromosome": "19", "position": 22157302, "genotype": 1},
        #                                      {"chromosome": "2", "position": 47739551, "genotype": 2},
        #                                      {"chromosome": "20", "position": 31024274, "genotype": 1},
        #                                      {"chromosome": "20", "position": 31386347, "genotype": 2},
        #                                      {"chromosome": "20", "position": 31386449, "genotype": 2},
        #                                      {"chromosome": "20", "position": 39797465, "genotype": 1},
        #                                      {"chromosome": "20", "position": 40743829, "genotype": 1},
        #                                      {"chromosome": "20", "position": 41306600, "genotype": 1},
        #                                      {"chromosome": "20", "position": 57478807, "genotype": 2},
        #                                      {"chromosome": "22", "position": 23653880, "genotype": 1},
        #                                      {"chromosome": "22", "position": 40058186, "genotype": 2},
        #                                      {"chromosome": "3", "position": 47162661, "genotype": 1},
        #                                      {"chromosome": "4", "position": 89670155, "genotype": 1},
        #                                      {"chromosome": "4", "position": 126397321, "genotype": 1},
        #                                      {"chromosome": "4", "position": 187516880, "genotype": 1},
        #                                      {"chromosome": "4", "position": 187534363, "genotype": 1},
        #                                      {"chromosome": "4", "position": 187534375, "genotype": 1},
        #                                      {"chromosome": "4", "position": 187538330, "genotype": 1},
        #                                      {"chromosome": "4", "position": 187557893, "genotype": 1},
        #                                      {"chromosome": "4", "position": 187629497, "genotype": 2},
        #                                      {"chromosome": "5", "position": 13701525, "genotype": 1},
        #                                      {"chromosome": "5", "position": 13701536, "genotype": 2},
        #                                      {"chromosome": "5", "position": 13845107, "genotype": 1},
        #                                      {"chromosome": "5", "position": 13864728, "genotype": 2},
        #                                      {"chromosome": "5", "position": 13913885, "genotype": 1},
        #                                      {"chromosome": "5", "position": 35861068, "genotype": 1},
        #                                      {"chromosome": "5", "position": 35871190, "genotype": 1},
        #                                      {"chromosome": "5", "position": 112162854, "genotype": 2},
        #                                      {"chromosome": "5", "position": 112164561, "genotype": 2},
        #                                      {"chromosome": "5", "position": 112175770, "genotype": 2},
        #                                      {"chromosome": "5", "position": 112177171, "genotype": 2},
        #                                      {"chromosome": "5", "position": 140168070, "genotype": 1},
        #                                      {"chromosome": "5", "position": 140168151, "genotype": 1},
        #                                      {"chromosome": "5", "position": 149460553, "genotype": 0},
        #                                      {"chromosome": "6", "position": 152464839, "genotype": 1},
        #                                      {"chromosome": "6", "position": 152466674, "genotype": 1},
        #                                      {"chromosome": "6", "position": 152469188, "genotype": 1},
        #                                      {"chromosome": "6", "position": 152665261, "genotype": 2},
        #                                      {"chromosome": "6", "position": 152772264, "genotype": 2},
        #                                      {"chromosome": "6", "position": 157405930, "genotype": 0},
        #                                      {"chromosome": "7", "position": 6026988, "genotype": 0},
        #                                      {"chromosome": "7", "position": 55214348, "genotype": 0},
        #                                      {"chromosome": "7", "position": 55238874, "genotype": 2},
        #                                      {"chromosome": "7", "position": 55249063, "genotype": 1},
        #                                      {"chromosome": "7", "position": 82453708, "genotype": 0},
        #                                      {"chromosome": "7", "position": 82455895, "genotype": 0},
        #                                      {"chromosome": "7", "position": 82581859, "genotype": 2},
        #                                      {"chromosome": "7", "position": 82582846, "genotype": 2},
        #                                      {"chromosome": "7", "position": 101844851, "genotype": 0},
        #                                      {"chromosome": "7", "position": 101892328, "genotype": 1},
        #                                      {"chromosome": "7", "position": 101917521, "genotype": 2},
        #                                      {"chromosome": "7", "position": 124499002, "genotype": 1},
        #                                      {"chromosome": "7", "position": 124499003, "genotype": 1},
        #                                      {"chromosome": "7", "position": 140426257, "genotype": 0},
        #                                      {"chromosome": "8", "position": 90958530, "genotype": 1},
        #                                      {"chromosome": "8", "position": 90967711, "genotype": 1},
        #                                      {"chromosome": "8", "position": 90995019, "genotype": 1},
        #                                      {"chromosome": "9", "position": 5081780, "genotype": 0},
        #                                      {"chromosome": "9", "position": 139391636, "genotype": 2},
        #                                      {"chromosome": "9", "position": 139397707, "genotype": 2},
        #                                      {"chromosome": "9", "position": 139407932, "genotype": 2},
        #                                      {"chromosome": "9", "position": 139410177, "genotype": 2},
        #                                      {"chromosome": "9", "position": 139412197, "genotype": 2},
        #                                      {"chromosome": "9", "position": 139418260, "genotype": 2}]}
        #        to_send = CompareProfiles(json.loads(stdout)["items"], json_to_send).for_cgd


#        to_send = CompareProfiles(stdout["items"], json_to_send).for_cgd

        # cmd, newfile = build_cmd(args)
        # stdout = run_cmd(cmd)
        #
        # profile_out.write(to_send)
        # profile_out.close()

        # Hold off on this until we get cgd_client worked out.
        # cmd = ['java8', '-jar', '/home/exacloud/clinical/installedTest/cgd_client-1.2.4.jar', '-c', '/home/exacloud/clinical/installedTest/cgd_client.properties', '-u', 'https://kdlwebuser01.ohsu.edu/cgd_next/service/run/180408_NS500390_0222_AHMMV5BGX5/barcodeId/A01/snpProfile']

    # Write and output log.  This is necessary so that Galaxy knows the process is over.
    outfile = open(args.stdout_log, 'w')
    outfile.write("The process has run.")
    outfile.close()

    # Clean up temp file.
    if args.endpoint == "uploadqcsheet" or args.endpoint == "uploadqcsheetrtwo" or args.endpoint == "cnvpdf":
        os.remove(newfile)


if __name__ == "__main__":
    main()

