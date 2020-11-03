#!/usr/bin/env python

## USAGE:
#   Metadata locates CollectInfo.txt with plate info
## NOTES:
#
## ARGUMENTS:
#   abfile
## RETURN:
#   metadata

import re
import argparse
import os


VERSION = "1.0.0.0"
def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='DSP metadata')
    parser.add_argument('--archive_path', type=str, help='path to HSR archive', default='/Volumes/HSRNanoString')

    args = parser.parse_args()
    return args


def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result

def parse_collectinfo(collectinfo):

    with open(collectinfo, 'r') as cifile:
        next(cifile) #first line is \n
        second_line = cifile.readline()
        print(second_line)
        cimatrix = []
        for line in cifile:
            items = line.strip().split()
            cimatrix.append(items)

    plate = second_line.strip().split("   ")[0].split(":")[1].strip()
    scan = second_line.strip().split("   ")[1].split(":")[1].strip()
    return {plate:[scan,cimatrix]}

def main():
    args = supply_args()
    #this might be better with pop find * -name "CollectInfo.txt", certainly is faster
    collectinfo_list = find_all("CollectInfo.txt", args.archive_path)

    cil={}
    #this one takes a while, searches directories for c
    for collectinfo in collectinfo_list:
        cil.update(parse_collectinfo(collectinfo))


    with open('dsp_metadata.txt', 'w') as fw:
        for key,value in cil.items():
            fw.write("Plate" '\t' + "Scan" + '\n')
            fw.write(key + '\t' + value[0] + '\n')





if __name__ == "__main__":
    main()