# wrapper for RNASeq Explorer
import os
from argparse import ArgumentParser
import pandas as pd
import subprocess
import logging

def explorer_argparse():
	
	parser = ArgumentParser(description='RNASeq Explorer: Gene Expression vs Drug Resistance')
	parser.add_argument("-e", "--exp_files", dest = "exp_files", nargs = "+", type=str, required=True, default=None, help="List of gene expression files (raw counts, fpkm)")
	parser.add_argument("-d", "--drug_files", dest = "drug_files", nargs = "+", type=str, required=True, default=None, help="List of drug response files")
	parser.add_argument("-a", "--algorithm", dest = "algorithm", required=True, type=str, default=None, help="Correlation Algorithm")
	parser.add_argument("-s", "--sensitivity", dest = "sensitivity", required=True, type=str, default=None, help="Sensitivity")
	parser.add_argument("-o1", "--output1", dest = "output1", required=True, type=str, default=None, help="Label Output File One")
	parser.add_argument("-o2", "--output2", dest = "output2", required=True, type=str, default=None, help="Label Output File Two")

	args = parser.parse_args()
	
	return args

# read column style expression
def consolidate_exp(file_list, workdir):
	new_ = []
	frame = pd.DataFrame()
	for f in file_list:
		with open(f) as this_file:
			df = pd.read_csv(this_file)
			if len(new_) == 0:
				new_.append(df)
			else:
				del df['gene']
				new_.append(df)
	frame = pd.concat(new_, axis=1)
	frame.to_csv(workdir+"/exp.csv", mode='w', index=False)

# read row style expression
def consolidate_drug(file_list, workdir):
	new_ = []
	new_f = open(workdir+"/drug.csv", "w")
	new_f.write("Patient: Patient ID,Specimen: Lab ID,Heme Malignancy: Diagnosis,Heme Malignancy: Specific Diagnosis,Inhibitor Panel Run: Inhibitor Panel,Inhibitor Panel Run: Run Type,Inhibitor Interpreted Result: Drug,Inhibitor Interpreted Result: Replicant,Inhibitor Interpreted Result: IC10,Inhibitor Interpreted Result: IC25,Inhibitor Interpreted Result: IC50,Inhibitor Interpreted Result: IC75,Inhibitor Interpreted Result: IC90,Inhibitor Interpreted Result: Area under the curve,Inhibitor Interpreted Result: Model Curve")
	for f in file_list:
		with open(f) as this_file:
			for line in this_file:
				if 'Lab ID' not in line:
					new_f.write(line)
	new_f.close()

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

def cmd_runner(args, scriptdir):
	cmd = "Rscript " + scriptdir + "/rnaseq_explorer.R " + args
	# print 'running with this command', cmd
	logging.info("RUNNING: %s" % (cmd))
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	return p.returncode

if __name__ == "__main__":
	args = explorer_argparse()
	one_file = vars(args)['exp_files'][0]
	scriptdir = os.path.dirname(os.path.realpath(__file__))
	wrkdir = os.path.dirname(os.path.realpath(one_file))
	# print 'the path to the file should be here', wrkdir
	cmd_args = cmd_builder(vars(args), wrkdir)
	# print 'running with these arguments', cmd_args
	cmd_runner(cmd_args, scriptdir)