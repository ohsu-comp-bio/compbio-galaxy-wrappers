#!/usr/bin/env python

# Galaxy wrapper for cgd_client.jar.
# JAVA8_PATH and CGD_CLIENT_CONFIG must be defined in the Galaxy contrib/ohsu_exacloud_env.sh file.
# 1.2.9.5 - Added support for chimeric junctions endpoint

import argparse
import json
import logging
import os
import shutil
import sys

import requests

from snp_profile import SnpProfile


# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32
if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

VERSION = '1.2.9.5'


def supply_args():

    parser = argparse.ArgumentParser(description='Galaxy wrapper for cgd_client.jar.')

    # parser.add_argument('stdout_log', help='Output file, mainly so that you can see if process succeeded in Galaxy.')
    parser.add_argument('--endpoint', help='CGD endpoint to send data, required.', required=True)
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
    parser.add_argument("--sampleid", help='Sample ID associated with import')
    parser.add_argument("--qcversion", help='Attached QC version, only useful for SeattleSeq.')
    parser.add_argument("--cnvcalls", help='CNV calls to be sent.')
    parser.add_argument("--cnvpdf", help='CNV PDF to be sent.')
    parser.add_argument("--cgd_client", help="Location of the cgd_client.")
    parser.add_argument("--cgd_config", help="Location of the cgd_client config file.")
    parser.add_argument("--include_chr", action="store_true", help="Include the chr prefix in reported variant output.")
    parser.add_argument("--servicebase", help="The service host name and port + service base. e.g. kdlwebprod02:8080/cgd")

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    args = parser.parse_args()
    return args


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

def run_cmd(logger, cmd, rdm):
    """
    Run the command via subprocess.
    
    The response from CGD can be a SimpleResponse with a message and error(s) or a list of json objects.  When we know 
    a json list will be received we don't parse it or log it, we just return it. 
    """
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if stderr:
        raise Exception(stderr)

    result = json.loads(stdout)

    if type(result) is list:
        logger.info(f"List response from CGD: {stdout[:200]}...")
        return result        
    elif 'errors' in result:
        logger.info(f"Error from CGD: {result}")
        if result['errors']:
            if result['errors'][0] == 'Could not find patient to provide previously reported variants' and rdm:
                return result
            else:
                raise Exception(result['errors'])
    elif 'message' in result and result['message'] == 'error_case_not_found':
        logger.warning(f"Case not found: {result}")
        return None
    elif 'message' in result and result['message'] == 'error_patient_not_found':
        # error_patient_not_found was probably changed to error_case_not_found so we should remove this reference
        raise ValueError(f"We didn't think 'error_patient_not_found' was used any more: {cmd} --> {result}")
    elif 'message' in result:
        logger.info(f"Message from CGD: {result}")
        return result
    else:
        raise ValueError(f"Response was not a list, message, or error: {result}")

def build_cmd(args):
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
    elif (args.endpoint == "annotationcomplete" or args.endpoint == "completeRun"
          or args.endpoint == "completeSampleRun" or args.endpoint == "annotate" or args.endpoint == "annotateRun"
          or args.endpoint == "annotateSampleRun" or args.endpoint == "reportedvariants"):
        # The cmd for this endpoint is already set, don't do anything.
        pass
    elif args.endpoint == "updatesamplerun" or args.endpoint == "metrics":
        cmd.extend(["-j", args.pipeline_out])
    elif args.endpoint == "snpProfile":
        cmd.extend(["-j", args.json_out])
    elif args.endpoint == 'requestVariants' or args.endpoint == 'uploadTranscriptEffects':
        cmd.extend(["-j", args.pipeline_out])
    elif args.endpoint == "none":
        cmd = [args.java8_path, "-jar", args.cgd_client, "-f", args.pipeline_out, "-u", args.cgd_url]
    elif not args.pipeline_out:
        raise ValueError(f"No file specified and endpoint parameter is unknown or missing: {args.endpoint}")
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

    return cmd, newfile


def write_vcf_header(outfile):
    """
    Write the VCF header.
    ##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
    :return:
    """
    outfile.write("##fileformat=VCFv4.2\n")
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def prepare_reported(outfile, regions, json_data, inc_chr=False):
    """

    :return:
    """
    empty = '.'
    for entry in json_data:
        if entry != 'message' and entry != 'errors':
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

    if not json_data:
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

    if 'errors' in json_data:
        if json_data['errors'][0] == 'Could not find patient to provide previously reported variants':
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


def check_conn(url, timeout=3):
    """
    Before trying to hit an endpoint, check to see if connection is open.
    :return:
    """
    try:
        request = requests.get(url, timeout=timeout)
        return True
    except (requests.ConnectionError, requests.Timeout) as exception:
        raise ConnectionError(exception)


def check_sample(samp):
    """
    Check to see if the sample is prefixed with RDM.
    :param samp:
    :return:
    """
    return samp.startswith('RDM-')


def main():
    args = supply_args()

    # Set up logger.
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    outfile = logging.FileHandler('stdout_log')
    outfile.setLevel(logging.DEBUG)
    logger.addHandler(outfile)

    # Also capture on stdout.
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # If the sample_id is passed, check to see if it starts with RDM.
    rdm = False
    if args.sampleid:
        rdm = check_sample(args.sampleid)

    # TODO: Consider moving SnpProfile into its own tool. Right now the --json_out parameter is used to give a filename to a json file 
    #       that is created here and then sent to CGD and also passed on to the next tool. This is a little confusing. I think there should
    #       be a parameter for saving the JSON list response from CGD and it would make sense to use "--json_out" for that. 
    if args.endpoint == 'snpProfile':
        json_to_send = SnpProfile(args.pipeline_out).geno_items
        with open(args.json_out, 'w') as to_cgd:
            json.dump(json_to_send, to_cgd)

    # Build the command.
    cmd, newfile = build_cmd(args)
    
    # Run the command and write command to log.
    logger.info("Running the following command:")
    logger.info('\t'.join(cmd))
    
    # TODO: This makes servicebase a required parameter, but cgd client could use its configuration to figure out the host 
    if check_conn(args.servicebase):
        json_response = run_cmd(logger, cmd, rdm)

    if not json_response:
        # There must have been a problem, it will have been logged 
        pass 
    elif args.endpoint == 'reportedvariants':
        vcf = open(args.report_vcf, 'w')
        regions = open(args.report_bed, 'w')
        write_vcf_header(vcf)
        prepare_reported(vcf, regions, json_response, args.include_chr)
    elif args.endpoint == 'requestVariants':
        write_response(json_response, args.json_out)

    outfile.close()

    # Clean up temp file.
    if (args.endpoint == "uploadqcsheet" or args.endpoint == "uploadqcsheetrtwo"
            or args.endpoint == "cnvpdf" or args.endpoint == "geneFusionReport"):
        os.remove(newfile)

def write_response(data, file_name):
    '''
    Write data to file
    '''
    with open(file_name, 'w') as file:
        json.dump(data, file, indent=2)

if __name__ == "__main__":
    main()
