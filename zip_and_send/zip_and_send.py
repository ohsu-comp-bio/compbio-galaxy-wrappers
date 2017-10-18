#!/usr/bin/env python

import argparse
import zipfile
import subprocess
import os

def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument(dest='gene_qc', help='')
    parser.add_argument(dest='int_qc', help='')
    parser.add_argument(dest='fastqc', help='')
    parser.add_argument(dest='vcf', help='')
    parser.add_argument(dest='seattle_seq', help='')
    parser.add_argument('-o', dest='output', required=True, help='')

    args = parser.parse_args()

    new_name = args.output.replace(' ', '_').replace('/', '_')
    file_path = "/tmp/" + new_name + ".zip"
    z = zipfile.ZipFile(file_path, 'w')

    z.write(args.gene_qc, arcname = "gene_qc.txt")
    z.write(args.int_qc, arcname = "int_qc.txt")
    z.write(args.fastqc, arcname = "fastqc.html")
    z.write(args.vcf, arcname = "variants.vcf")
    z.write(args.seattle_seq, arcname = "seattle_seq.tsv")

    z.close()

    cgd_url = "http://10.78.40.45:8080/cgd/chimera/intervalqc/" + new_name
    cmd = "java -jar /mnt/lustre1/CompBio/cgdtools/cgd_client.jar -f %s -u %s" % (file_path, cgd_url)
    process = subprocess.call(args=cmd, shell=True, stdout=subprocess.PIPE)

    os.remove(file_path)
    

if __name__ == "__main__":
    main()

