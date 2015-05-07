# wrapper for RNASeq Explorer

from argparse import ArgumentParser
import pandas as pd
import subprocess
import logging
import rpy2.robjects as robjects

# import pandas.rpy.common as com

def explorer_argparse():
	
	parser = ArgumentParser(description='RNASeq Explorer: Gene Expression vs Drug Resistance')
	parser.add_argument("-e", "--exp_files", dest = "exp_files", nargs = "+", type=str, required=True, default=None, help="List of gene expression files (raw counts, fpkm)")
	parser.add_argument("-d", "--drug_files", dest = "drug_files", nargs = "+", type=str, required=True, default=None, help="List of drug response files")
	
	args = parser.parse_args()
	
	return args

def read_exp(file_list):
	master = robjects.DataFrame({})
	for f in file_list:
		if len(master) == 0:
			master = robjects.DataFrame.from_csvfile(f, sep=",", header = True)
		else:
			df = robjects.DataFrame.from_csvfile(f, sep=",", header = True)
			sub = df.rx(robjects.IntVector(range(2,len(df)+1)))
			master = robjects.DataFrame.cbind(master, sub)
	return master

def read_drug(file_list):
	master = robjects.DataFrame({})
	for f in file_list:
		if len(master) == 0:
			master = robjects.DataFrame.from_csvfile(f, sep=",", header = True)
		else:
			df = robjects.DataFrame.from_csvfile(f, sep = ",", header = True)
			master = robjects.DataFrame.rbind(master, df)
	return master

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
	RNASeq = read_exp(args.exp_files)
	drug_screen = read_drug(args.drug_files)