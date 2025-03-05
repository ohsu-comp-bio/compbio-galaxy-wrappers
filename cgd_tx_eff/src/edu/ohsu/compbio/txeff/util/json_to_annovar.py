'''
Created on Mar 4, 2025

@author: pleyte
'''
import argparse
import csv
import json
import logging.config
from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig

VERSION = '0.0.1'

class JsonToAnnovar(object):
    '''
    Convert a list of variants in json format to an Annovar avinput file 
    '''
    def read(self, json_file):
        """
        Read json and return list of variants 
        """
        with open(json_file, 'r') as file:
            return json.load(file)

    def write(self, variants: list, output_filename):
        """
        Write list of variants to tab delimited avinput file 
        """
        with open(output_filename, "w") as file:
            tsv_writer = csv.writer(file, delimiter='\t')
            for x in variants:
                positionEnd =  x['positionStart'] + len(x['referenceBase']) - 1
                tsv_writer.writerow([self._get_chromosome(x['chromosome']), x['positionStart'], positionEnd, x['referenceBase'], x['variantBase']])
    
    def _get_chromosome(self, chromosome):
        return chromosome.replace('chr', '')
    
def _parse_args():
    parser = argparse.ArgumentParser(description='Read variants from a json formatted file and write out an Annovar input file.')

    parser.add_argument('-i', '--in_json', 
                help='json list of variants', 
                type=argparse.FileType('r'), 
                required=True)
    
    parser.add_argument('-o', '--out_avinput', 
                help = 'avinput file', 
                type = argparse.FileType('w'), 
                required = False,
                default = "out.avinput")

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    
    return parser.parse_args()

def _main():
    logging.config.dictConfig(TfxLogConfig().utility_config)
    logger = logging.getLogger("edu.ohsu.compbio.txeff.util.json_to_annovar")

    args = _parse_args()

    j2a = JsonToAnnovar()
    variants = j2a.read(args.in_json.name)
    
    j2a.write(variants, args.out_avinput.name)
    
    logger.info(f"Wrote {len(variants)} variants to {args.out_avinput.name}")
    
if __name__ == '__main__':
    _main()    