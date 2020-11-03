#!/usr/bin/env python

## USAGE:
#   transfer_nanostring.py --config nanostring_config.json --batch 20190314_208421591020

## ARGUMENTS:
#
## RETURN:
#
import argparse
import os
import sys
import shutil
import zipfile
import json
import re
from nanostring_client import nCounter

VERSION="1.0.0"


def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='Transfer batch from ftp to mounted X drive')
    parser.add_argument('--config', type=str, help='nanostring_config.json', default='nanostring_config.json')
    parser.add_argument('--batch', type=str, help='<batch id>_RCC.ZIP')
    args = parser.parse_args()
    return args

def prep_dest_dir(batch, dest_data_dir, ss_dir):
    """
    Creates directory if directory does exist.
    Copies samplesheet file if it exists in samplesheet and doesn't already exist
    """

    data_dir = os.path.join(dest_data_dir, batch)
    if os.path.isdir(data_dir):
        print("Directory exists")
    else:
        os.mkdir(data_dir)
    #find samplesheet
    ssregex = re.compile(batch + ".*" + "_samplesheet.txt")
    ss=None
    for file in os.listdir(ss_dir):
        if ssregex.match(file):
            ss = file
            print(ss)
    if ss is not None:
        # if samplesheet doesn't already exist then copy
        if not os.path.isfile(os.path.join(data_dir, ss)):
            shutil.copy(os.path.join(ss_dir, ss), data_dir, follow_symlinks=True)
    else:
        print("No samplesheet exists for " + batch)




def main():
    args = supply_args()
    print(args.batch)
    with open(args.config, 'r') as cfgfile:
        cfg = json.load(cfgfile)

    source_host = cfg['nCounterFTP']['host']
    source_user = cfg['nCounterFTP']['user']
    source_password = cfg['nCounterFTP']['password']
    dest_data_dir = cfg['Xdrive']['automated_datadir']
    ss_dir = cfg['Xdrive']['samplesheet_dir']

    data_dir = os.path.join(dest_data_dir, args.batch)
    zfile = args.batch + "_RCC.ZIP"
    zfile_path = os.path.join(data_dir, zfile)

    with nCounter(source_host, source_user, source_password) as ftp_host:
        print("ftp host established")

        #ftp_host = ftputil.FTPHost("10.125.46.25", "technician", "MAX123")
        # check if zip file and directory structure already exist in automated data
        if os.path.isfile(zfile_path):
            print(os.path.join(args.batch + "already exists and files exist"))
        else:
            # create automated_data directory and copy samplesheet
            prep_dest_dir(args.batch, dest_data_dir, ss_dir)
            ftp_host.download_batch(args.batch, data_dir)

        if os.path.isfile(zfile_path):
            with zipfile.ZipFile(zfile_path,'r') as zobj:
                zobj.extractall(data_dir)
    #move zip file to zipfile directory, jsut remove, RCCData already being backedup
    # shutil.move(zfile_path, os.path.join(args.dest_data_dir, "RCCData", zfile))
    os.remove(zfile_path)


if __name__ == "__main__":
    main()
