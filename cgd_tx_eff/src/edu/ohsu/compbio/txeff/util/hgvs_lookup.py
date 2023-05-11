'''
Use the HGVS python library to lookup a single variant. This script allows us to test and debug 
HGVS processing of a variant without having to run a VCF through tx_eff_control.py.  

Created on Mar. 21, 2023

@author: pleyte
'''
import logging
import os
from argparse import ArgumentParser
import hgvs.assemblymapper
import hgvs.dataproviders.uta
from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig
from edu.ohsu.compbio.txeff import tx_eff_hgvs
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript

class HgvsLookup(object):
    def __init__(self):
        '''
        Constructor
        '''
        self.logger = logging.getLogger(__name__)
        
        if os.environ.get('HGVS_SEQREPO_DIR') is None:
            self.logger.warning("The HGVS_SEQREPO_DIR environment variable is not defined. The remote seqrepo database will be used.")
        else:
            logging.info(f"Using SeqRepo {os.environ.get('HGVS_SEQREPO_DIR')}")
                
        if os.environ.get('UTA_DB_URL') is None:
            logging.warning("The UTA_DB_URL environment variable is not defined. The remote UTA database will be used.")
        else:
            logging.info(f"Using UTA Database at {os.environ.get('UTA_DB_URL')}")

    def find_genotype(self, genotype: str):
        '''
        Search UTA for transcripts associated with the genotype (e.g. 7-12345-C-G)
        '''
        chromosome, position, ref, alt = genotype.split('-')
        variant = VariantTranscript(chromosome, position, ref, alt)
        
        transcripts = tx_eff_hgvs._lookup_hgvs_transcripts([variant])
        logging.info(f"Found {len(transcripts)} transcripts associated with {genotype}")
        
        self.print_result_transcripts(transcripts)

    def find_gene(self, gene:str):
        hdp = hgvs.dataproviders.uta.connect()
        gene = hdp.get_gene_info(gene)
        print(f'Gene: {gene}')
            
    def print_result_transcripts(self, transcripts: list):
        for transcript in transcripts:
            print(f'Variant: {transcript.chromosome}-{transcript.position}-{transcript.reference}-{transcript.alt}')
            print(f'  Gene: {transcript.hgnc_gene}')
            print(f'  g.: {transcript.sequence_variant}')
            print(f'  c.: {transcript.hgvs_c_dot}')
            print(f'  p.: {transcript.hgvs_p_dot_three}')
            print(f'  Transcript: {transcript.refseq_transcript}')
            print(f'  Protein Transcript: {transcript.protein_transcript}')
            print(f'  Splicing: {transcript.splicing if transcript.splicing else "No"}')
            print('--------')
        
if __name__ == '__main__':
    logging.config.dictConfig(TfxLogConfig().utility_config)
    parser = ArgumentParser(description='Use the HGVS python library to lookup a single variant')
    parser.add_argument("-v", "--variant_genotype", help="Variant genotype (e.g. 7-123456-C-G)", type=str, required=False)
    parser.add_argument("-g", "--gene", help="Gene name", type=str, required=False)

    args = parser.parse_args()
    
    hgvs_lookup = HgvsLookup()
    
    if args.variant_genotype is not None:
        hgvs_lookup.find_genotype(args.variant_genotype)
    elif args.gene is not None:
        hgvs_lookup.find_gene(args.gene)
    else:
        print("Please specify genotype or gene")
