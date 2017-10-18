#!/usr/bin/env python

# Remove the UMI from R2 of input paired-end FASTQ files regardless of gzip status.
# Write new files, currently not gzipped.
# Include parameter that allows removal of different size UMIs.
# Include option to either remove or not remove ligating sequence.
# USAGE: python umi_reformat_fastq.py
# CODED BY: John Letaw

import argparse
import subprocess
import os

VERSION = '0.3.0'

def supply_args():
    """                                                                                                                                                                     
    Populate args.                                                                                                                                                          
    https://docs.python.org/2.7/library/argparse.html                                                                                                                       
    """
    parser = argparse.ArgumentParser(description='Include UMIs in the proper position of the FASTQ file.')
    parser.add_argument('input_fastq', help='Input FASTQ R1 File')
    parser.add_argument('input_fastq_2', help='Input FASTQ R2 File')
    parser.add_argument('output_fastq', help='Output FASTQ R1 File')
    parser.add_argument('output_fastq_2', help='Output FASTQ R2 File')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def run_cmd(cmd, outfile):
    """
    Run command.
    """
    print('Running the following command:')
    print('\t'.join(cmd))

    p = subprocess.call(cmd, stdout=outfile)


def file_type(filename):
    """
    Determine the file type based on bytes at start of file.
    """

    magic_dict = {
        "\x1f\x8b\x08": "gz",
        "\x42\x5a\x68": "bz2",
        "\x50\x4b\x03\x04": "zip"
        }

    max_len = max(len(x) for x in magic_dict)

    with open(filename) as f:
        file_start = f.read(max_len)
    for magic, filetype in magic_dict.items():
        if file_start.startswith(magic):
            return filetype


def zcat_cmd(filename):
    """
    Construct command to run a zcat.
    """
    return ['gzip', '-cd', filename]


def create_handle(filetype, temp_out, infq):
    """
    Return a file handle based on whether we are working with a gzipped file or not.
    """
    if filetype == "gz":
        cmd = zcat_cmd(infq)
        with open(temp_out, 'w') as outfile:
            run_cmd(cmd, outfile)
        handle_fq = open(temp_out, 'r')
    else:
        handle_fq = open(infq, 'r')

    return handle_fq


def main():

    args = supply_args()
    temp_out = '/tmp/umi_reformat_fastq_temp'
    temp_out_2 = '/tmp/umi_reformat_fastq_temp2'

    filetype = file_type(args.input_fastq)
    filetype_2 = file_type(args.input_fastq_2)
    handle_fq = create_handle(filetype, temp_out, args.input_fastq)
    handle_fq_2 = create_handle(filetype_2, temp_out_2, args.input_fastq_2)
    handle_out = open(args.output_fastq, 'w')
    handle_out_2 = open(args.output_fastq_2, 'w')

    umis = {}

    with handle_fq_2 as myfq2:
        while myfq2:
            idline = myfq2.readline()
            seq = myfq2.readline()
            spacer = myfq2.readline()
            quals = myfq2.readline()

            umi = seq[:12]
            new_seq = seq[22:]
            new_qual = quals[22:]
            id_begin = idline.split(' ')[0]
            umis[id_begin] = umi

            if idline:
                new_idline = idline.split()[0] + ':' + umi + ' ' + idline.split()[1] + '\n'
            else:
                break

            print(new_idline)
            
            handle_out_2.write(new_idline)
            handle_out_2.write(new_seq)
            handle_out_2.write(spacer)
            handle_out_2.write(new_qual)

    with handle_fq as myfq:
        while myfq:
            idline = myfq.readline()
            seq = myfq.readline()
            spacer = myfq.readline()
            quals = myfq.readline()

            id_begin = idline.split(' ')[0]

            if idline:
                new_idline = idline.split()[0] + ':' + umis[id_begin] + ' ' + idline.split()[1] + '\n'
            else:
                break

            handle_out.write(new_idline)
            handle_out.write(seq)
            handle_out.write(spacer)
            handle_out.write(quals)


    try:
        os.remove(temp_out)
        os.remove(temp_out_2)
    except:
        pass

    handle_out.close()
    handle_out_2.close()
    handle_fq.close()
    handle_fq_2.close()


if __name__ == "__main__":
    main()
