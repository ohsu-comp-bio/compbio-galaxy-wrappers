'''
Created on May 18, 2022

@author: pleyte
'''

import argparse
import logging
from edu.ohsu.compbio.txeff import tx_eff_annovar, tx_eff_hgvs, tx_eff_vcf
from edu.ohsu.compbio.txeff.tx_eff_ccds import TxEffCcds

VERSION = '0.2.5'

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
    
    parser.add_argument('-r', '--require_match',
                        help='Only add transcript effects to variants that are found by both Annovar and HGVS/UTA.',
                        action='store_true')

    parser.add_argument('-i', '--in_vcf', 
                help='Input VCF', 
                type=argparse.FileType('r'), 
                required=True)

    parser.add_argument('-c', '--ccds_map', 
                help='Input CSV or GFF with Annovar-to-CCDS mappings. Use refseq_to_ccds.py to convert GFF to CSV.', 
                type=argparse.FileType('r'), 
                required=True)
    
    parser.add_argument('--annovar_variant_function', 
                        help='Annovar variant_function file',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('--annovar_exonic_variant_function', 
                        help='Annovar exonic_variant_function file',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('-o', '--out_vcf', 
                    help='Output VCF', 
                    type=argparse.FileType('w'), 
                    required=True)

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    
    args = parser.parse_args()
    
    return args

def _main():
    '''
    main function
    '''            
    args = _parse_args()

    # Use tx_eff_annovar to read annovar records 
    annovar_records = tx_eff_annovar.get_annovar_records(args.annovar_variant_function.name, args.annovar_exonic_variant_function.name)

    # Use tx_eff_hgvs to fix the nomenclature
    tx_eff_hgvs.identify_hgvs_datasources() 
    merged_transcripts = tx_eff_hgvs.get_updated_hgvs_transcritpts(annovar_records, require_match = args.require_match)

    # Use tx_eff_ccds to add CCDS transcripts
    tx_eff_ccds = TxEffCcds(args.ccds_map.name)
    tx_eff_ccds.add_ccds_transcripts(merged_transcripts)
     
    # Use tx_eff_vcf to write the transcript effects to a VCF
    tx_eff_vcf.create_vcf_with_transcript_effects(args.in_vcf.name, args.out_vcf.name, merged_transcripts)

if __name__ == '__main__':
    _main()
    
    
    
    