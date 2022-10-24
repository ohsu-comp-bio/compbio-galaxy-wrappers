'''
Created on Sep. 27, 2022
This script reads a bunch of genotypes from one file, looks those genotypes up in an existing annovar input file, 
and writes each of those annovar lines to a new file.
@author: pleyte
'''

import argparse
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

stream_handler = logging.StreamHandler()
logging_format = '%(levelname)s: [%(filename)s:%(lineno)s - %(funcName)s()]: %(message)s'

stream_format = logging.Formatter(logging_format)
stream_handler.setFormatter(stream_format)
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)

class Genotype(object):
    def __init__(self, genotype:str):
        parts = genotype.split('-', 4)
        self.chromosome = parts[0]
        self.position = parts[1]
        self.reference = parts[2]
        self.alt = parts[3]        
        
    def __str__(self, *args, **kwargs):
        return self.chromosome + '-' + self.position + '-' + self.reference + '-' + self.alt
    
class FindGenotypeInAvinput(object):
    '''
    classdocs
    '''
    def __init__(self, genotype_file, av_input_file, av_output_file):
        '''
        Constructor
        '''
        self.genotype_file = genotype_file
        self.av_input_file = av_input_file
        self.av_output_file = av_output_file
    
    def _get_genotypes(self):
        '''
        Return a list of genotype object read from the input file
        '''
        logger.info(f"Reading genotypes from {self.genotype_file}")
        genotypes = list()
        with open(self.genotype_file,"r") as file:
            for line in file:                
                genotypes.append(Genotype(line.strip()))
        
        return genotypes
    
    def _is_line_matches_genotype(self, tab_delimited_line, genotype: Genotype):        
        '''
        return true if the tab delimiited line defines the genotype
        '''        
        # Split the tab delimited line into multiple parts 
        line_parts = tab_delimited_line.split("\t")
        if(line_parts[0] == 'chr'+genotype.chromosome and 
           line_parts[1] == genotype.position and 
           line_parts[3] == genotype.reference and 
           line_parts[4] == genotype.alt):
            return True
        else:
            return False
                
    def _get_matching_avinput(self, genotypes):
        '''
        '''
        logger.info(f"Attempting to find {len(genotypes)} genotypes in {self.av_input_file}")
        annovar_lines = list()
        
        logger.info(f"There are {len(genotypes)} genotypes")
        
        with open(self.av_input_file, "r") as file:
            for line in file:
                for genotype in genotypes:
                    if(self._is_line_matches_genotype(line, genotype)):
                        logger.debug(f"Found match for {genotype}")
                        annovar_lines.append(line)
                        genotypes.remove(genotype)
                        break
                    
        logger.info(f"Matched {len(annovar_lines)} annovar genotypes. Did not match {len(genotypes)} genotypes.")                
                    
        return annovar_lines
    
    def _write_output_file(self, avinput_lines):
        '''
        '''
        logger.info(f"Writing {len(avinput_lines)} records to {self.av_output_file}")
        with open(self.av_output_file, "w") as file:
            for line in avinput_lines:
                file.write(line)
                
    
    def run(self):
        '''
        '''
        genotypes = self._get_genotypes()
        avinput_lines = self._get_matching_avinput(genotypes)
        self._write_output_file(avinput_lines)
        
def _parse_args():
    '''
    '''
    parser = argparse.ArgumentParser(description='')
    
    parser.add_argument('-g', '--genotype_file',  
                        help='Input file with one genotype per line (of the form 1-123-A-C)',
                        type=argparse.FileType('r'),
                        required=True)
    
    parser.add_argument('-i', '--av_input_file',  
                        help='Annovar input file (eg annovar.avinput)',
                        type=argparse.FileType('r'),
                        required=True)
    
    parser.add_argument('-o', '--av_output_file', 
                        help='Output file name', 
                        type=argparse.FileType('w'), 
                        required=True)
    
    args = parser.parse_args()
    
    return args

def _main():
    args = _parse_args()
    finder = FindGenotypeInAvinput(args.genotype_file.name, args.av_input_file.name, args.av_output_file.name)
    finder.run()
    
if __name__ == '__main__':
    _main()