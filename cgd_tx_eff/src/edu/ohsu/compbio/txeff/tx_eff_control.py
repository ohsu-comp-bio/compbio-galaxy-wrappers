'''
Created on May 18, 2022

@author: pleyte
'''

import argparse
import logging.config
import time
from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig
from edu.ohsu.compbio.txeff.util.tx_eff_pysam import PysamTxEff
from edu.ohsu.compbio.txeff import tx_eff_annovar, tx_eff_hgvs, tx_eff_vcf
from edu.ohsu.compbio.txeff.tx_eff_ccds import TxEffCcds

VERSION = '0.5.5'

def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Read variants from a VCF, add Annovar variant effects, correct nomenclature using HGVS, and write out a new VCF.')

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

    parser.add_argument('--reference_fasta',
                        help='Reference genome, in FASTA format.  The associated index is expected to be in the same directory.',
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
    print(f"tfx_cgd {VERSION} is starting...")
    
    logging.config.dictConfig(TfxLogConfig().log_config)
    print(f"Log level={logging.root.getEffectiveLevel()}, file={logging.root.handlers[0].baseFilename}")

    start_time = time.time()
    
    args = _parse_args()

    # Use tx_eff_annovar to read annovar records 
    annovar_records = tx_eff_annovar.get_annovar_records(args.annovar_variant_function.name, args.annovar_exonic_variant_function.name)

    # Load the reference genome with pysam.
    pysam_file = PysamTxEff(args.reference_fasta)

    # Use tx_eff_hgvs to fix the nomenclature
    tx_eff_hgvs.identify_hgvs_datasources()

    merged_transcripts = tx_eff_hgvs.get_updated_hgvs_transcripts(annovar_records, pysam_file)
    # Close the reference FASTA
    pysam_file.my_fasta.close()

    # Use tx_eff_ccds to add CCDS transcripts
    tx_eff_ccds = TxEffCcds(args.ccds_map.name)
    tx_eff_ccds.add_ccds_transcripts(merged_transcripts)
     
    # Use tx_eff_vcf to write the transcript effects to a VCF
    tx_eff_vcf.create_vcf_with_transcript_effects(args.in_vcf.name, args.out_vcf.name, merged_transcripts)
    
    print(f"Wrote {len(merged_transcripts)} transcripts to {args.out_vcf.name}")

    end_time = time.time()
    total_time = end_time - start_time
    time_str = time.strftime("%Mm:%Ss", time.gmtime(total_time))
    print(f"Completed in {time_str}")

if __name__ == '__main__':
    _main()