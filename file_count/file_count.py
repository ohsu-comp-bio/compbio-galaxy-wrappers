#!/usr/bin/env python

from __future__ import print_function
import argparse
import os

VERSION = '0.3.0'


def supply_args():
    """                                                                                                                            Populate args,
    https://docs.python.org/2.7/library/argparse.html                                                                              """
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('input_cnt', help='Input number of files that should be in directory.')
    parser.add_argument('base_dir', help='Directory to count files within.')
    parser.add_argument('file_list', help='File in .list format to write to.')
    parser.add_argument('--mother', help='KDL accession of mother, implies trio.')
    parser.add_argument('--father', help='KDL accession of father, implies trio.')
    parser.add_argument('--proband', help='KDL accession of proband, implies trio.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    args = parser.parse_args()
    return args

def full_path(base_dir, fname):
    """
    Return absolute PATH based on base_dir and file name.
    :return:
    """
    return '/'.join([base_dir.rstrip('/'), fname])

def trio_determine(mother, father, proband):
    """
    Based on arguments, determine if we are going to be processing a trio.
    :return:
    """
    if mother and father and proband:
        return True
    return False

def check_trio_members(dir_list, *args):
    """
    Check to see if the member accession exists in the dirlist.
    :param member:
    :return:
    """
    new_dir_list = []
    for accession in args:
        for entry in dir_list:
            if accession == entry.rstrip('/').split('/')[-1].split('.')[0]:
                new_dir_list.append(entry)
    return new_dir_list

def main():
    args = supply_args()
    is_trio = trio_determine(args.mother, args.father, args.proband)
    dir_list = [name for name in os.listdir(args.base_dir)]
    if is_trio:
        dir_list = check_trio_members(dir_list, args.mother, args.father, args.proband)
    dir_len = len(dir_list)
    status = (dir_len == int(args.input_cnt))
    if not status:
        raise Exception("Counts are unequal(" + str(dir_len) + ',' + args.input_cnt + "), stopping.")
    else:
        with open(args.file_list, 'w') as to_write:
            for entry in dir_list:
                to_write.write(full_path(args.base_dir, entry))
                to_write.write('\n')

if __name__ == "__main__":
    main()
