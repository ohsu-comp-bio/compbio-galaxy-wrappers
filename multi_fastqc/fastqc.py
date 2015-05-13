# wrapper for RNASeq Explorer
import os
from argparse import ArgumentParser
import pandas as pd
import subprocess
import logging
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def explorer_argparse():
    
    parser = ArgumentParser(description='RNASeq Explorer: Gene Expression vs Drug Resistance')
    parser.add_argument("-i", "--input", dest = "input", nargs = "+", type=str, required=True, default=None, help="List of fastqc_data.txt files.")
    parser.add_argument("-o1", "--output1", dest = "output1", required=True, type=str, default=None, help="Label Output File One")
    # parser.add_argument("-o2", "--output2", dest = "output2", required=True, type=str, default=None, help="Label Output File Two")

    args = parser.parse_args()
    
    return args

def get_modules(sing_file):
    file_dict = {}
    with open(sing_file) as f:  
        this_mod = ''
        for line in f.readlines():
            line = line.strip() # trim endline
            if (line[:2] == ">>" and line[:12] != ">>END_MODULE"): # for each module grab summary data
                module = line[2:-5] # grab module name
                status = line[-4:] # and overall status pass/warn/fail
                if module not in file_dict:
                    file_dict[module] = []
                    this_mod = module
            # removing headers
            elif (line[:2] != ">>" and line[:1] != "#"): # grab details under each module
                cols = line.split('\t')
                file_dict[this_mod].append(cols)
                #file_dict[this_mod] = file_dict[this_mod] + line + '\n'
    return file_dict

def plot_per_seq_quality(in_files):
    with PdfPages('plot_per_seq_quality.pdf') as pdf:
        for input_f in in_files:
            f_dict = get_modules(input_f)
            arr_ = np.array(f_dict['Per sequence quality scores'])
            plt.plot(arr_[:,0], arr_[:,1], label=input_f)
        plt.legend(loc='upper left')
        plt.ylabel('Count')
        plt.xlabel('Mean Sequence Quality')
        plt.title('Per Sequence Quality Scores')
        pdf.savefig()
        plt.close()


def cmd_builder(arg_dict, workdir):
    consolidate_exp(arg_dict['exp_files'], workdir)
    consolidate_drug(arg_dict['drug_files'], workdir)
    exp = workdir + "/exp.csv"
    drug = workdir + "/drug.csv"
    return " ".join([drug,
        exp,
        arg_dict['algorithm'],
        arg_dict['sensitivity'],
        arg_dict['output1'],
        arg_dict['output2']])

def cmd_runner(args, workdir):
    cmd = "Rscript " + workdir + "/rnaseq_explorer.R " + args
    logging.info("RUNNING: %s" % (cmd))
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return p.returncode

if __name__ == "__main__":
    args = explorer_argparse()
    wrkdir = os.path.dirname(os.path.realpath(__file__))
    plot_per_seq_quality(vars(args)['input'])
