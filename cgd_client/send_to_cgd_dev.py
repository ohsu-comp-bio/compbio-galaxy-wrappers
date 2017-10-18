#!/usr/bin/env python

# Galaxy wrapper for cgd_client.jar.
# JAVA8_PATH and CGD_CLIENT_CONFIG must be defined in the Galaxy contrib/ohsu_exacloud_env.sh file.
# USAGE: send_to_cgd.py -h
# CODED BY: John Letaw

from __future__ import print_function
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

VERSION = '1.2.4.1'


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
    parser.add_argument("--pipeline_out", help='')
    parser.add_argument("--cgd_url", help='')
    parser.add_argument("--runid", help='')
    parser.add_argument("--barcodeid", help='')
    parser.add_argument("--qcversion", help='')
    parser.add_argument("--cnvcalls", help='')
    parser.add_argument("--cnvpdf", help='')
    parser.add_argument("--cgd_client", help="Location of the cgd_client.")
    parser.add_argument("--cgd_config", help="Location of the cgd_client config file.")
    parser.add_argument("--tissue", help="Which tissue is the client receiving data for.  This is most applicable to tumor/normal workflows.")
    parser.add_argument("--servicebase", help="The service host name and port + service base. e.g. kdlwebprod02:8080/cgd")

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
    print(stdout)
    if stderr:
        raise Exception(stderr)
#        eprint(stderr)
    return stdout


def build_cmd(args):
    """
    Build the command that will send data to the CGD.
    """
    cmd = [args.java8_path, '-jar', args.cgd_client, "-n", args.endpoint, "-c", args.cgd_config]
    newfile = ""
        
    if args.servicebase:
        cmd.extend(["-s", args.servicebase])
    if args.endpoint == "uploadqcsheet" or args.endpoint == "uploadqcsheetrtwo":
        newfile = rename_fastqc_output(args.runid, args.barcodeid, args.endpoint, 'html') # For FastQC files, we will create a new file name.
        print("Copying to " + newfile)
        shutil.copyfile(args.pipeline_out, newfile)
        cmd.extend(["-f", newfile])
    elif args.endpoint == "cnvpdf":
        newfile = rename_fastqc_output(args.runid, args.barcodeid, args.endpoint, 'pdf')
        print("Copying to " + newfile)
        shutil.copyfile(args.pipeline_out, newfile)
        cmd.extend(["-f", newfile])
    elif args.endpoint == "annotationcomplete" or args.endpoint == \
            "annotate" or args.endpoint == "reportedvariants":
        pass # The cmd for this endpoint is already set, don't do anything.
    elif args.endpoint == "updatesamplerun":
        cmd.extend(["-j", args.pipeline_out])
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
              
    #cmd.append("-d")

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
    for entry in json.loads(stdout):
        chrom = entry['chromosome'][3:]
        pos = entry['positionStart']
        ref = entry['referenceBase']
        alt = entry['variantBase']
        empty = '.'
        outfile.write('\t'.join([chrom, str(pos), empty, ref, alt, empty, empty, empty]))
        outfile.write('\n')

        start = pos - 1
        regions.write('\t'.join([chrom, str(start), str(pos)]))
        regions.write('\n')

    outfile.close()
    regions.close()


def main():
    
    args = supply_args()
    # Build the command.
    cmd, newfile = build_cmd(args)
    # Run the command.
    stdout = run_cmd(cmd)

    if args.endpoint == 'reportedvariants':
        vcf = open(args.report_vcf, 'w')
        regions = open(args.report_bed, 'w')
        write_vcf_header(vcf)
        print(stdout)
        prepare_reported(vcf, regions, stdout)

    # Write and output log.  This is necessary so that Galaxy knows the process is over.
    outfile = open(args.stdout_log, 'w')
    outfile.write("The process has run.")
    outfile.close()

    # Clean up temp file.
    if args.endpoint == "uploadqcsheet" or args.endpoint == "uploadqcsheetrtwo" or args.endpoint == "cnvpdf":
        os.remove(newfile)

if __name__ == "__main__":
    main()

