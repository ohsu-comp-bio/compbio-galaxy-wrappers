#!/usr/bin/env python

# Galaxy wrapper for cgd_client.jar.
# JAVA8_PATH and CGD_CLIENT_CONFIG must be defined in the Galaxy contrib/ohsu_exacloud_env.sh file.
# USAGE: send_to_cgd.py -h
# CODED BY: John Letaw

# Not providing this code until we have finalized SNP profile reqs for lab.
# from snp_profile import SnpProfile, CompareProfiles
import argparse
import json
import logging
import os
import sys
import shutil

# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

VERSION = '1.2.9.0'


def supply_args():

    parser = argparse.ArgumentParser(description='Galaxy wrapper for cgd_client.jar.')

    # parser.add_argument('stdout_log', help='Output file, mainly so that you can see if process succeeded in Galaxy.')
    parser.add_argument('--endpoint', help='CGD endpoint to send data, required.')
    parser.add_argument('--java8_path', help='Specify java 8 path, in the case you have multiple java installations.')
    parser.add_argument('--report_vcf', help='Output VCF if utilizing '
                                             'reportvariants endpoint.')
    parser.add_argument('--report_bed', help='Output BED if utilizing '
                                             'reportvariants endpoint.')
    parser.add_argument('--json_out', help='JSON will be written to this file.', default='cgd_profile_to_send.json')
    parser.add_argument("--pipeline_out", help='Output to send to CGD.')
    parser.add_argument("--cgd_url", help='CGD URL to send data to.')
    parser.add_argument("--runid", help='Run ID associated with import.')
    parser.add_argument("--barcodeid", help='Barcode ID associated with import')
    parser.add_argument("--qcversion", help='Attached QC version, only useful for SeattleSeq.')
    parser.add_argument("--cnvcalls", help='CNV calls to be sent.')
    parser.add_argument("--cnvpdf", help='CNV PDF to be sent.')
    parser.add_argument("--cgd_client", help="Location of the cgd_client.")
    parser.add_argument("--cgd_config", help="Location of the cgd_client config file.")
    parser.add_argument("--include_chr", action="store_true", help="Include the chr prefix in reported variant output.")
    parser.add_argument("--tissue", help="Which tissue is the client receiving data for. Deprecated.")
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
    elif endpoint == "geneFusionReport":
        newfile = "/tmp/" + '_'.join([runid, barcodeid]) + ext
    else:
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
        if json.loads(stdout)['message'] == 'error_patient_not_found':
            return None
        elif json.loads(stdout)['message'] != 'ok':
            raise Exception(json.loads(stdout)['message'])
    return stdout


def build_cmd(args, recvd_prof=False):
    """
    Build the command that will send data to the CGD.
    """
    if args.java8_path:
        cmd = [args.java8_path, '-jar', args.cgd_client, "-n", args.endpoint, "-c", args.cgd_config]
    else:
        cmd = ['java', '-jar', args.cgd_client, "-n", args.endpoint, "-c", args.cgd_config]
    newfile = ""

    if args.servicebase:
        cmd.extend(["-s", args.servicebase])
    if args.endpoint == "uploadqcsheet" or args.endpoint == "uploadqcsheetrtwo":
        # For FastQC files, we will create a new file name.
        newfile = rename_fastqc_output(args.runid, args.barcodeid, args.endpoint, 'html')
        logging.info("Copying FastQC to " + newfile)
        shutil.copyfile(args.pipeline_out, newfile)
        cmd.extend(["-f", newfile])
    elif args.endpoint == "cnvpdf":
        newfile = rename_fastqc_output(args.runid, args.barcodeid, args.endpoint, 'pdf')
        logging.info("Copying CNV PDF to " + newfile)
        shutil.copyfile(args.pipeline_out, newfile)
        cmd.extend(["-f", newfile])
    elif args.endpoint == "geneFusionReport":
        newfile = rename_fastqc_output(args.runid, args.barcodeid, args.endpoint, 'html')
        logging.info("Copying gene fusion HTML report to " + newfile)
        shutil.copyfile(args.pipeline_out, newfile)
        cmd.extend(["-f", newfile])
    elif args.endpoint == "annotationcomplete" \
            or args.endpoint == "completeRun" \
            or args.endpoint == "completeSampleRun" \
            or args.endpoint == "annotate" \
            or args.endpoint == "annotateRun" \
            or args.endpoint == "annotateSampleRun" \
            or args.endpoint == "reportedvariants" \
            or (args.endpoint == "snpProfile" and not recvd_prof):
        # The cmd for this endpoint is already set, don't do anything.
        pass
    elif args.endpoint == "updatesamplerun" or args.endpoint == "metrics":
        cmd.extend(["-j", args.pipeline_out])
    elif args.endpoint == "snpProfile" and recvd_prof:
        cmd.extend(["-j", args.json_out])
    elif args.endpoint == "none":
        cmd = [args.java8_path, "-jar", args.cgd_client, "-f", args.pipeline_out, "-u", args.cgd_url]
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

    logging.info("The client is receiving data for a " + args.tissue + " sample.")

    return cmd, newfile


def write_vcf_header(outfile):
    """
    Write the VCF header.
    ##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
    :return:
    """
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def prepare_reported(outfile, regions, stdout, inc_chr=False):
    """

    :return:
    """
    empty = '.'
    for entry in json.loads(stdout):
        if json.loads(stdout):
            if inc_chr:
                chrom = entry['chromosome']
            else:
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
        if inc_chr:
            outfile.write('\t'.join(['chr1', '3', empty, 'T', 'C', empty, empty, empty]))
            outfile.write('\n')
            regions.write('\t'.join(['chr1', '1', '2']))
            regions.write('\n')
        else:
            outfile.write('\t'.join(['1', '3', empty, 'T', 'C', empty, empty, empty]))
            outfile.write('\n')
            regions.write('\t'.join(['1', '1', '2']))
            regions.write('\n')

    outfile.close()
    regions.close()


def main():
    args = supply_args()
    # outfile = open(args.stdout_log, 'w')
    # Set up logger.
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    outfile = logging.FileHandler('stdout_log')
    outfile.setLevel(logging.DEBUG)
    logger.addHandler(outfile)
    # logging.basicConfig(filename=args.stdout_log, level=logging.DEBUG)
    # Build the command.
    cmd, newfile = build_cmd(args, False)
    # Run the command and write command to log.
    logger.info("Running the following command:")
    logger.info('\t'.join(cmd))
    stdout = run_cmd(cmd)

    # Write CGD return json to log.
    logger.info("From CGD:")
    logger.info(json.loads(stdout))

    if args.endpoint == 'reportedvariants':
        vcf = open(args.report_vcf, 'w')
        regions = open(args.report_bed, 'w')
        write_vcf_header(vcf)
        prepare_reported(vcf, regions, stdout, args.include_chr)

        # Get the current profile ready.
        # We run one command here, then another down below.

    if args.endpoint == 'snpProfile' and stdout:
        # Command has been run to retrieve the profile.  Now, we need to figure out if the overlap of the
        # two profiles match.
        # If the match, no need to run anything else.
        # If they don't match, error.
        # If they match, but there are loci that aren't contained within the CGD profile, send the updated loci to CGD.
        # If there is no profile in the CGD, send the created profile.
        # stdout = {"message": "ok", "items": [{"chromosome": "1", "position": 2488153, "genotype": -1},
        #                                      {"chromosome": "1", "position": 16256007, "genotype": 1},
        #                                      {"chromosome": "1", "position": 36937059, "genotype": 2},
        #                                      {"chromosome": "1", "position": 65321388, "genotype": 2},
        #                                      {"chromosome": "9", "position": 139418260, "genotype": 2}]}
        # TODO: Deal with (Exception: error_patient_not_found)
        # This is on hold until the lab decides whether they want to use it.  We will be implementing a version of the
        # profiling tool for HOP, which will compare all samples against each other and produce a table.
        pass

        # json_to_send = SnpProfile(args.pipeline_out).geno_items
        # cgd_json = json.loads(stdout)['items']
        # compare_snps = CompareProfiles(cgd_json, json_to_send)
        # with open(args.json_out, 'w') as to_cgd:
        #     json.dump(compare_snps.for_cgd, to_cgd)
        # cmd, newfile = build_cmd(args, True)
        # stdout = run_cmd(cmd)
        # out_metric = {'snp_profile_pvalue': compare_snps.pvalue,
        # 'snp_profile_total': compare_snps.total, 'snp_profile_mismatch': compare_snps.mismatch}
        # json.dump(out_metric, outfile)

    outfile.close()

    # Clean up temp file.
    if (args.endpoint == "uploadqcsheet" or args.endpoint == "uploadqcsheetrtwo"
            or args.endpoint == "cnvpdf" or args.endpoint == "geneFusionReport"):
        os.remove(newfile)


if __name__ == "__main__":
    main()
