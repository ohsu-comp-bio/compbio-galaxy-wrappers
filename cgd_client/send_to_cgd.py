#!/usr/bin/env python

###  <command interpreter="python">send_to_cgd.py
###    $pipeline_output $endpoint $cgd_url $output $runid $barcodeid $qcversion
###  </command>

### Galaxy wrapper for cgd_client.jar.
### CGD_CLIENT is hard coded, but this is not expected to move.
# VERSION: 1.0.1
# TOOL VERSION: 1.2.3

import argparse
import subprocess
from subprocess import Popen, STDOUT, PIPE
import os
import sys
import shutil

def renameOutput(runid, barcodeid, endpoint):

    """
    CGD needs the filename to be restructured.
    Applies to FastQC only.
    """

    if endpoint == "uploadqcsheet":
        newfile = "/tmp/" + '_'.join([runid, barcodeid, "R1"]) + ".html"
    elif endpoint == "uploadqcsheetrtwo":
        newfile = "/tmp/" + '_'.join([runid, barcodeid, "R2"]) + ".html"
    else:
        print("Not sending FastQC.")
        return None

    return newfile

def splitUrl(url, n):
    
    """
    Split the manual URL and take n elements.
    """

    return url.split('/')[-n:]

def main():
    
#    CGD_CLIENT="/opt/installed/cgd_client-1.0.7.jar"
#    CGD_CLIENT="/home/exacloud/clinical/installedTest/cgd_client-1.0.9.jar"

    # Java 8 required from cgd_client 1.2 and on.
    JAVA = "/opt/installed/java8/bin/java"

    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--pipeline_out", help='')
    parser.add_argument("--cgd_url", help='')
    parser.add_argument(dest='stdout_log', help='')
    parser.add_argument(dest='endpoint', help='')
    parser.add_argument("--runid", help='')
    parser.add_argument("--barcodeid", help='')
    parser.add_argument("--qcversion", help='')
    parser.add_argument("--cgd_client", help="Location of the cgd_client.", default="/home/exacloud/clinical/installedTest/cgd_client-1.2.3.jar")
#    parser.add_argument("--config", help="Location of the cgd_client config file.", default="/home/exacloud/clinical/installedTest/cgd_client.properties")
    parser.add_argument("--tissue", help="Which tissue is the client receiving data for.  This is most applicable to tumor/normal workflows.")
    parser.add_argument("--servicebase", help="The service host name and port + service base. e.g. kdlwebprod02:8080/cgd")
#    parser.add_argument("--is_fastqc", action='store_true', help='If this is a FastQC file, and a URL is manually set, a temp file must be created so that CGD understands the file format.')

    args = parser.parse_args()

    CONFIG = "/home/exacloud/clinical/installedTest/cgd_client.properties"

    cmd = [JAVA, "-jar", args.cgd_client, "-n", args.endpoint, "-c", CONFIG]

    if args.servicebase:
        cmd.extend(["-s", args.servicebase])
    if args.endpoint == "uploadqcsheet" or args.endpoint == "uploadqcsheetrtwo":
        newfile = renameOutput(args.runid, args.barcodeid, args.endpoint) # For FastQC files, we will create a new file name.
        print("Copying to " + newfile)
        shutil.copyfile(args.pipeline_out, newfile)
        cmd.extend(["-f", newfile])
    elif args.endpoint == "annotationcomplete" or args.endpoint == "annotate":
        pass # The cmd for this endpoint is already set, don't do anything.
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
              
    cmd.append("-d")

    print("The client is receiving data for a " + args.tissue + " sample.")

    print("We are running this command:")
    print(' '.join(cmd))

    print(cmd)
    proc = subprocess.call(cmd)

    outfile = open(args.stdout_log, 'w')
    outfile.write("The process has run.")
    outfile.close()

    ## Clean up temp file.
    if args.endpoint == "uploadqcsheet" or args.endpoint == "uploadqcsheetrtwo":
        os.remove(newfile)

if __name__ == "__main__":
    main()

