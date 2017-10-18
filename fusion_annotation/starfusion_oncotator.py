#!/usr/bin/env python
# USAGE: python starfusion_oncotator.py <full path to starfusionOutput> <oncotator_dbdir>
# EXAMPLE: python ~python/starfusion_oncotator.py
#      star-fusion.fusion_candidates.final.abridged.FFPM
#      AnnotationSources/Oncotator/oncotator_v1_ds_April052016
# Script takes in starfusion results and writes two FICTITIOUS maflite files. Then runs oncotator on mafs
# FICTITIOUS, why? Turns breakpoint to single nt deletions because Oncotator annotates snps/indels only.
# Two files are created one for left & right breakpoints and mafs indicate single nt deletion
# ARGS:
#   starfusionOutput (str): star-fusion.fusion_candidates.final.abridged.FFPM
# RETURN:
#   auto generated oncotator output files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Including this because is it quite useful...
# https://docs.python.org/2/library/subprocess.html
# https://github.com/google/python-subprocess32

import sys
import subprocess
import os
def convert_starfusion_to_bedpe(starfusionOutput):
    #with open(starfusionOutput, 'r') as starfusion_ffpm:
    bedpe =[]
    for line in starfusionOutput:
        if not line.startswith('#'):
            vals = line.strip().split()
            valsL = vals[5].split(':')
            valsR = vals[7].split(':')
            chrL = valsL[0]
            posL = valsL[1]
            geneL = vals[4].split('^')[0]
            strandL = valsL[2]
            chrR = valsR[0]
            posR = valsR[1]
            geneR = vals[6].split('^')[0]
            strandR = valsR[2]
            linebedpe=[chrL, str(int(posL) - 1), posL, chrR, str(int(posR) - 1), posR, '-'.join([geneL, geneR]), '0', strandL, strandR, vals[1], vals[2], vals[3], vals[4].split('^')[1], vals[4].split('^')[1]]
            linebedpe.extend(vals[8:15])
            bedpe.append(linebedpe)
    return bedpe


def convert_starfusion_to_maflites(starfusionOutput):
    """
    Convert starfusion results into maflite format. Simplest input for oncotator. Two lists, left & right breakpoint
    Args:
        starfusionOutput (file): opened for read file type
    :Returns:
        maflite_leftls (list of list)
        maflite_rightls (list of list)
    """
    maflite_leftls = []
    maflite_rightls = []
    for line in starfusionOutput:
            if not line.startswith('#'):
                    vals = line.strip().split()
                    valsL = vals[5].split(':')
                    valsR = vals[7].split(':')
                    chrL = valsL[0]
                    posL = valsL[1]
                    dinucL = vals[9][0]
                    chrR = valsR[0]
                    posR = valsR[1]
                    dinucR = vals[11][0]
                    maflite_leftls.append([chrL, str(int(posL)-1), posL, dinucL, '-'])
                    maflite_rightls.append([chrR, str(int(posR)-1), posR, dinucR, '-'])
    return maflite_leftls, maflite_rightls


def write_mafs(maflite, outfile):
    """
    :param maflite
    :param outfile
    Input maflite is a list of list
    Close the handle after processing.
    """
    outfile.write('\t'.join(["chr", "start", "end", "ref_allele", "alt_allele\n"]))
    for i in maflite:
        # handle_out.write(str(item) for item in maflite[i])
        outfile.writelines('\t'.join(i))
        outfile.write('\n')

def runOncotator(oncotator_dbdir, inputmaf, output_filename):
    """
    Runs bash command to run Oncotator. Originally wanted to cut output, but I can't seem to figure out how to get Oncotator to output into stdout, rather than write to a file.
    So for the sake of not adding a single bash command after python, then to fusionannotation.R I kept the command here.
    Args:
        oncotator_dbdir (str): directory string
        inputmaf: name of maflite file
        output_filename: name of oncotator output file

    Returns: nothing, writes to directory
    Todo:
    Change hardcoded db_dir to argparse
    Figure out how write oncotator to stdout, and Popen to cut -f whatever columns I need from oncotator
    """
    #output_file = anno_dir + output_filename
    output_file = output_filename
    log_file = "oncotator.log"
    cmd = ['Oncotator', '--input_format=MAFLITE', '--db-dir', oncotator_dbdir, '--output_format=SIMPLE_TSV', '--log_name', log_file,
           inputmaf, output_file, 'hg19']
    #cmd = ['Oncotator', '-i', 'MAFLITE', '--db-dir', oncotator_dbdir, '-o', 'SIMPLE_TSV', '--log_name', log_file, inputmaf, output_file, 'hg19']
    subprocess.call(cmd)
    # I don't know how to make this output go to stdout instead write to file.
    # So I'm just going to cut it in R


    # Oncotator -i MAFLITE --db-dir oncotator_v1_ds_April052016 -o SIMPLE_TSV --log_name log.out left_brkpnt.maf oncotated_left_output hg19


def main():
    """
    Commented files are testing
    # starfusionOutput=open("/Users/patterja/Workspace/FusionAnnotation/star-fusion.fusion_candidates.final.abridged.FFPM", 'r')
    # handleL_out=open("maflite_left.txt", 'w')
    :return:
    """

    starfusion_filename = sys.argv[1]
    oncotator_dbdir = sys.argv[2]

    #Version 2 no longer specifies directory structure
    #anno_dir = os.path.dirname(starfusion_filename) + '/anno_files/'
    #try:
    #    os.stat(anno_dir)
    #except:
    #    os.mkdir(anno_dir)
    mafleft_filename = "left_brkpnt.maf"
    mafright_filname = "right_brkpnt.maf"


    with open(starfusion_filename, 'r') as starfusionOutput:
        maflite_leftls, maflite_rightls=convert_starfusion_to_maflites(starfusionOutput)

    # Oncotator Input
    with open(mafleft_filename, 'w') as outfile_mafleft:
        write_mafs(maflite_leftls, outfile_mafleft)

    with open(mafright_filname, 'w') as outfile_mafright:
        write_mafs(maflite_rightls, outfile_mafright)

    # Run Oncotator, output is automatic
    # oncotator_dbdir = '/home/users/patterja/BioCoders/DataResources/AnnotationSources/Oncotator/oncotator_v1_ds_April052016'
    # see fusion_annotation.py

    runOncotator(oncotator_dbdir, mafleft_filename, "oncotated_left_output")
    runOncotator(oncotator_dbdir, mafright_filname, "oncotated_right_output")


if __name__ == "__main__":
    main()
