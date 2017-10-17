#!/usr/bin/env python

### Command to wrap:
### R --vanilla --slave -f /var/lib/cnv_service/scripts/cna-analysis-25Apr14.R --args job86/amplicons.txt job86/samples.txt job86/counts.csv
### Need to do something with the output files so they are seen in Galaxy.  They look like this:
### calls.RDM.15.01189.tumor.txt
### out.RDM.15.01189.tumor.pdf


import sys
import subprocess
import os
import shutil

VERSION = '0.1.0'

def build_cmd(program, *argv):
    """
    Build the command we will run.
    """

    build_cmd = ["R", "--vanilla", "--slave", "-f", program, "--args"]
    for arg in argv:
        build_cmd.append(arg)
#    build_cmd.append('--label-sig-genes-only=T')

    return build_cmd


def main():

    program = "/home/groups/clinical/installedTest/cna-analysis-25Apr14.R"
    amplicons = sys.argv[1]
    samples = sys.argv[2]
    counts = sys.argv[3]
    pdf_outfile = sys.argv[4]
    calls_outfile = sys.argv[5]

    cmd = build_cmd(program, amplicons, samples, counts)
    print("Running the following command: ")
    print(' '.join(cmd))
    proc = subprocess.call(cmd)

    work_dir = os.listdir('.')
    print(work_dir)
    pdf_out = [f for f in work_dir if "tumor.pdf" in f]
    print(pdf_out)
    calls_out = [f for f in work_dir if "tumor.txt" in f]
    print(calls_out)

    if len(pdf_out) == 1:
        shutil.copy(pdf_out[0], pdf_outfile)
    else:
        raise Exception("Too many elements with extension tumor.pdf")
    
    if len(calls_out) == 1:
        shutil.copy(calls_out[0], calls_outfile)
    else:
        raise Exception("Too many elements with extension tumor.txt")


if __name__ == "__main__":
    main()
