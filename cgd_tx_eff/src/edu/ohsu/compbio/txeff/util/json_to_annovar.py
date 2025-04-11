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
    def read(self, in_stream):
        """
        Read json and return list of variants 
        """
        return json.load(in_stream)

    def write(self, variants: list, out_stream):
        """
        Write list of variants to tab delimited avinput file 
        """
        tsv_writer = csv.writer(out_stream, delimiter='\t')
        for x in variants:
            position_end =  x['positionStart'] + len(x['referenceBase']) - 1
            id = f"id:{x['genomicVariantId']}"
            tsv_writer.writerow([self._get_chromosome(x['chromosome']), 
                                 x['positionStart'], 
                                 position_end, 
                                 x['referenceBase'], 
                                 x['variantBase'], 
                                 id])
    
    def _get_chromosome(self, chromosome):
        return chromosome.replace('chr', '')
    
def _parse_args():
    parser = argparse.ArgumentParser(description='Read variants from a json formatted file and write out an Annovar input file.')

    parser.add_argument('-i', '--in', 
                dest="input",
                help='json list of variants exported from CGD', 
                type=argparse.FileType('r'), 
                required=True)
    
    parser.add_argument('-o', '--out', 
                dest="output",
                help='avinput file', 
                type=argparse.FileType('w'), 
                required=False,
                default="out.avinput")

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    
    return parser.parse_args()

def _main():
    logging.config.dictConfig(TfxLogConfig().utility_config)
    logger = logging.getLogger("edu.ohsu.compbio.txeff.util.json_to_annovar")

    args = _parse_args()

    j2a = JsonToAnnovar()
    variants = j2a.read(args.input)
    
    j2a.write(variants, args.output)
    
    logger.info(f"Wrote {len(variants)} variants to {args.output.name}")
    
if __name__ == '__main__':
    _main()    