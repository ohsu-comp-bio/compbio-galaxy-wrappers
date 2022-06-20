'''
Created on May 18, 2022

@author: pleyte
'''

import argparse
import logging
import sys
import os
from edu.ohsu.compbio.txeff import tx_eff_annovar, tx_eff_hgvs, tx_eff_vcf

VERSION = '0.0.1'

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

stream_handler = logging.StreamHandler()
logging_format = '%(levelname)s: [%(filename)s:%(lineno)s - %(funcName)s()]: %(message)s'

stream_format = logging.Formatter(logging_format)
stream_handler.setFormatter(stream_format)
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)

def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Read variants from a VCF, add Annovar variant effects, correct nomenclature using HGVS, and write out a new VCF.')
    
    parser.add_argument('-i', '--in_vcf', 
                help='Input VCF', 
                type=argparse.FileType('r'), 
                required=True)

    parser.add_argument('-o', '--out_vcf', 
                    help='Output VCF', 
                    type=argparse.FileType('w'), 
                    required=True)
    
    parser.add_argument('annovar_file', 
                    nargs=argparse.REMAINDER, 
                    help='One or more variant_function and exonic_variant_function Annovar files.',
                    type=argparse.FileType('r'))

    parser.add_argument('-r', '--require_match',
                        help='Only add transcript effects to variants that are found by both Annovar and HGVS/UTA.',
                        action='store_true')

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    
    args = parser.parse_args()
    
    if len(args.annovar_file) == 0:
        print("At least one Annovar input file is necessary")
        sys.exit(os.EX_CONFIG)

    return args

def _main():
    '''
    main function
    '''            
    args = _parse_args()

    # Use tx_eff_annovar to read annovar records 
    annovar_file_names = [x.name for x in args.annovar_file]
    annovar_records = tx_eff_annovar.get_annovar_records(annovar_file_names)
    
    # Use tx_eff_hgvs to fix the nomenclature
    tx_eff_hgvs.identify_hgvs_datasources() 
    merged_transcripts = tx_eff_hgvs.get_updated_hgvs_transcritpts(annovar_records, require_match = args.require_match)
    
    # Use tx_eff_vcf to write the transcript effects to a VCF
    tx_eff_vcf.create_vcf_with_transcript_effects(args.in_vcf.name, args.out_vcf.name, merged_transcripts)

if __name__ == '__main__':
    _main()
    
    
    
    