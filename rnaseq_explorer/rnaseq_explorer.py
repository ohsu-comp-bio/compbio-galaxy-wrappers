# wrapper for RNASeq Explorer

from argparse import ArgumentParser
import pandas as pd
import subprocess
import logging

def explorer_argparse():
	
	parser = ArgumentParser(description='RNASeq Explorer: Gene Expression vs Drug Resistance')
	parser.add_argument("-e", "--exp_files", dest = "exp_files", nargs = "+", type=str, required=True, default=None, help="List of gene expression files (raw counts, fpkm)")
	parser.add_argument("-d", "--drug_files", dest = "drug_files", nargs = "+", type=str, required=True, default=None, help="List of drug response files")
	
	args = parser.parse_args()
	
	return args

# read column style expression
def consolidate_exp(file_list):
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
	frame.to_csv("exp.csv", mode='w', index=False)

# read row style expression
def consolidate_drug(file_list):
	new_ = []
	new_f = open("drug.csv", "w")
	new_f.write("Patient: Patient ID,Specimen: Lab ID,Heme Malignancy: Diagnosis,Heme Malignancy: Specific Diagnosis,Inhibitor Panel Run: Inhibitor Panel,Inhibitor Panel Run: Run Type,Inhibitor Interpreted Result: Drug,Inhibitor Interpreted Result: Replicant,Inhibitor Interpreted Result: IC10,Inhibitor Interpreted Result: IC25,Inhibitor Interpreted Result: IC50,Inhibitor Interpreted Result: IC75,Inhibitor Interpreted Result: IC90,Inhibitor Interpreted Result: Area under the curve,Inhibitor Interpreted Result: Model Curve")
	for f in file_list:
		with open(f) as this_file:
			for line in this_file:
				if 'Lab ID' not in line:
					new_f.write(line)
	new_f.close()

def prepare_data():
	cmd = 'R -e "shiny::runApp(port=8002)"'
	logging.info("RUNNING: %s" % (cmd))
	print 'running', cmd
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = p.communicate()
	print stdout
	if len(stderr):
		print stderr
	return p.returncode

if __name__ == "__main__":
	args = explorer_argparse()
	consolidate_exp(args.exp_files)
	consolidate_drug(args.drug_files)
	prepare_data()