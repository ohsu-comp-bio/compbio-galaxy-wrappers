#!/usr/bin/env python

from bioblend import galaxy
import getpass
import argparse
import json
import shutil
import os
import grp
import errno
import pysam
import requests.packages.urllib3
import subprocess
requests.packages.urllib3.disable_warnings()

VERSION='0.2.1'

def mkdir_p(path):
    ### Emulate mkdir -p functionality.

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def get_access_groups():
    """
    Find all access groups associated with the current user.
    """
    my_groups = [(grp.getgrgid(g).gr_name, int(g)) for g in os.getgroups()]
    return my_groups


def set_new_group(my_groups, group):
    """
    Set a new current group for the process to run under."
    """
    for entry in my_groups:
        if entry[0] == group:
            os.setgid(entry[1])


def index_bam(bai_origin, dest_bai, dest_bam):
    ### Index the BAM based on whether it already exists in Galaxy or not.

    if bai_origin != None:
        print("Copying " + bai_origin + " to " + dest_bai + ".")
        shutil.copyfile(bai_origin, dest_bai) 
        os.chmod(dest_bai, 0440)
    else:
        print("BAI not found in Galaxy, reindexing...")
        pysam.index(dest_bam, dest_bai)


# These functions run commands using sg and subprocess.  For users that
# have more than 16 groups.
def run_sg_copy_cmd(origin, dest, group):
    """
    Build the command, using sg, that will run a copy operation as a specific group.
    Don't know of python libs for this, will use subprocess.
    -a = -rlptgoD in rsync
    """

#    cmd = "sg %s \"cp %s %s\"" % (group, origin, dest)
    cmd = "sg %s \"rsync --chmod=u+rw,g+rw,o-rwx %s %s\"" % (group, origin, dest) 
    print(cmd)
    process = subprocess.call(args=cmd, shell=True, stdout=subprocess.PIPE)

#    print(subprocess.check_output(["pidof", cmd]))


def run_sg_index_cmd(filename, group):
    """
    Perform an indexing step under sg utility.
    """

    cmd = "sg %s \"samtools index %s\"" % (group, filename)
    print(cmd)
    process = subprocess.call(args=cmd, shell=True, stdout=subprocess.PIPE)


def main():

    ### Store current user, so we can connect them to a Galaxy API key.
    curr_user = getpass.getuser()

    ### Set argparse options.
    parser = argparse.ArgumentParser(description='Move BAM files from Galaxy to warm storage.')
#    parser.add_argument('--dummy_input', help="Dummy input file.")
    parser.add_argument('--galaxy_url', default='https://exaclinical.ohsu.edu/galaxy', help='URL of the Galaxy instance.')
    ### Temporarily set a default history id to test with.
    parser.add_argument('--history_id', default='8d4d7622a593869c', help='Galaxy history id, defined after a Galaxy history is created.')
    parser.add_argument('--sample_id', default='DNA-15-01448-1', help='Illumina SampleSheet sample id.  This will be used to create the BAM and BAI files.')
    parser.add_argument('--bam_path', help='Path where BAM files will be deposited.')
    parser.add_argument('--run_id', help='A subdirectory with the same name as the run_id will be created in the bam_path directory.')

    parser.add_argument('--input', help="Input file to be moved.")
    parser.add_argument('--output', default='/tmp/default.log', help='Outfile')
    args = parser.parse_args()

    api_file = '/home/users/' + curr_user + '/.galaxy_api_key.json'
    with open(api_file, 'r') as f:
        api_key = json.load(f)

    # Code to find access groups and set a default access group based on where the BAM files are going.
    my_groups = get_access_groups()
    # Choose your path, based NFS groups number limitations.
    if len(my_groups) > 16:
        use_sg = True
    else:
        use_sg = False

#    set_new_group(my_groups, "CorlessLab")

    print(my_groups)
    print(len(my_groups))

    gi = galaxy.GalaxyInstance(url=args.galaxy_url, key=api_key[curr_user])

    this_hist = gi.histories.show_history(args.history_id, contents=True)

    for entry in this_hist:
        # Make the name an argument, or do something else, so we can move other BAM files that Print Reads.
        if "Print Reads" in entry['name'] and "BAM" in entry['name'] and entry['deleted'] != True:

            dataset_id = entry['id']
            bam_origin = gi.datasets.show_dataset(dataset_id)['file_path']
            bai_origin = gi.datasets.show_dataset(dataset_id)['metadata_bam_index']

            ### Change this behavior to automatically create an index with Samtools if there is none.
            if bam_origin == args.input:

                new_path = args.bam_path + args.run_id + '/'
                dest_bam = new_path + args.sample_id + '.bam'
                dest_bai = new_path + args.sample_id + '.bai'

                print("Copying " + bam_origin + " to " + dest_bam + ".")
                if use_sg == False:
                    mkdir_p(new_path)
                    if not os.path.isfile(dest_bam):
                        shutil.copyfile(bam_origin, dest_bam)
                        os.chmod(dest_bam, 0440)
                ### Check to see if the index file was found in Galaxy, if not, make one.
                        index_bam(bai_origin, dest_bai, dest_bam)
                    else:
                        if not os.path.isfile(dest_bai):
                            print("BAM file has been copied, but there is no index.")
                            index_bam(bai_origin, dest_bai, dest_bam)
                else:
                    # Convert to argument.
                    mkdir_p(new_path)
                    run_sg_copy_cmd(bam_origin, dest_bam, "CorlessLab")
                    run_sg_index_cmd(dest_bam, "CorlessLab")

            elif bam_origin == None:
                raise Exception("No BAM filepath found in Galaxy for " + bam_origin)


    handle_out = open(args.output, 'w')
    handle_out.close()

if __name__ == "__main__":
    main()






