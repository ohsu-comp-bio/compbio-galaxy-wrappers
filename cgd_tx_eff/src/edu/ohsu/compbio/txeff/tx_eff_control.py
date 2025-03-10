'''
Created on May 18, 2022

@author: pleyte
'''
import argparse
from collections import Counter
import datetime
import json
import logging.config
import sys
import time

import biocommons.seqrepo
import hgvs

from edu.ohsu.compbio.txeff.tx_eff_annotate import TxEffAnnotate
from edu.ohsu.compbio.txeff.tx_eff_annovar import TxEffAnnovar
from edu.ohsu.compbio.txeff.tx_eff_ccds import TxEffCcds
from edu.ohsu.compbio.txeff.tx_eff_hgvs import TxEffHgvs
from edu.ohsu.compbio.txeff.tx_eff_writer import TxEffWriter
from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig
from edu.ohsu.compbio.txeff.util.tx_eff_pysam import PysamTxEff
from edu.ohsu.compbio.txeff.variant import Variant


VERSION = '0.7.9'

def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Read variants from a VCF, add Annovar variant effects, correct nomenclature using HGVS, and write out a new VCF.')
    
    parser.add_argument('--annovar_variant_function', 
                        help='Annovar variant_function file',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('--annovar_exonic_variant_function', 
                        help='Annovar exonic_variant_function file',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('--ccds_map', 
                help='Input CSV or GFF with Annovar-to-CCDS mappings.', 
                type=argparse.FileType('r'), 
                required=True)

    parser.add_argument('--reference_fasta',
                        help='Reference genome, in FASTA format.  The associated index is expected to be in the same directory.',
                        required=True)

    parser.add_argument('--out',
                    help='Output file  with variants and their transcript effects (json)', 
                    type=argparse.FileType('w'), 
                    required=True)

    parser.add_argument('--variants',
                help='list of variants received from cgd', 
                type=argparse.FileType('r'), 
                required=False)
    
    parser.add_argument('--benchmark',
                    help='Write benchmarking information to file', 
                    action='store_true')
    
    parser.add_argument('--sequence_source',
                    help='Method for looking up reference sequences (can be a url, path, or "ncbi")',
                    required=False)
    
    parser.add_argument('--threads', default = 1, type=int,
                    help="Number of threads to use", )
    
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    
    args = parser.parse_args()
    
    if args.sequence_source \
       and not (args.sequence_source.lower().startswith('http') or args.sequence_source.startswith("/") or args.sequence_source == 'ncbi'):
        raise ValueError(f'--sequence_source is expected to be a url, directory, or "ncbi": "{args.sequence_source}"')

    if args.threads > 1 and args.benchmark:
        raise ValueError('Benchmarking is not possible when multiple threads are used')
    
    return args

def _get_time(real_seconds, user_seconds):
    """
    Take the real time and user time in seconds and return a string that breaks the time into hours, minutes, and seconds. 
    """
    
    dt_real = datetime.timedelta(seconds = real_seconds)
    dt_user = datetime.timedelta(seconds = user_seconds)
    
    return str(dt_real), str(dt_user)

def _is_same_variant(v1, v2):
    """
    Return true if the variants have the same genotype. v1 is a dict and v2 is a Variant object. 
    """
    return (v1['chromosome'].replace('chr', '') == v2.chromosome and 
            v1['positionStart'] == v2.position and
            v1['referenceBase'] == v2.reference and
            v1['variantBase'] == v2.alt)

    
def _check_output_for_all_variants(logger, input_variants, transcripts):
    """
    Compare the list of variants exported from CGD with the transcript effects being output, and identify
    any variants that are in the input but not the output. 

    """
    output_variants = {x._id: Variant(x.chromosome, x.position, x.reference, x.alt, x._id) for x in transcripts}

    if len(input_variants) != len(output_variants):
        logger.info(f"The number of input variants is not equal to the number of output variants: {len(input_variants)} != {len(output_variants)}")

    counter = Counter()    
    for x in input_variants:
        counter['Total'] += 1
        output_variant = output_variants.get(str(x['genomicVariantId']))
        if not output_variant:
            logger.info(f"Genomic variant id {x['genomicVariantId']} not found in output variant transcripts")
            counter['Input variants not found in output'] += 1
        elif not _is_same_variant(x, output_variant):
            counter['Variants with same id but different genotype'] += 1
        else:
            counter['Variants found and have matching genotype'] += 1
            
    logger.info(f"Input/output variant comparison: {counter}")

def _main():
    '''
    main function
    '''
    print(f"tfx_cgd {VERSION} is starting...")
    print("Python version: " + sys.version.replace('\n', ''))
    print(f"biocommons.seqrepo version: {biocommons.seqrepo.__version__}")
    print(f"hgvs version: {hgvs.__version__}")
    
    logging.config.dictConfig(TfxLogConfig().log_config)
    logger = logging.getLogger("edu.ohsu.compbio.txeff.tx_eff_control")

    if logging.root.handlers[0].stream: 
        output = str(logging.root.handlers[0].stream.name)
    else:
        output = logging.root.handlers[0].baseFilename

    print(f"Log level={logging.root.getEffectiveLevel()}, output={output}")

    start_time_real = time.perf_counter()
    start_time_user = time.process_time()

    args = _parse_args()

    print(f"Thread count: {args.threads}")

    # Use tx_eff_annovar to read annovar records
    annovar_records = TxEffAnnovar().get_annovar_records(args.annovar_variant_function.name, args.annovar_exonic_variant_function.name)

    # Load the reference genome with pysam.
    pysam = PysamTxEff(args.reference_fasta)

    # Load the refseq to ccds mappings
    tx_eff_ccds = TxEffCcds(args.ccds_map.name)

    # Look for additional transcripts in the HGVA/UTA database and merge them with the annovar records.
    with TxEffHgvs(pysam = pysam, 
                   refseq_ccds_map = tx_eff_ccds.refseq_to_ccds_map,
                   sequence_source = args.sequence_source,
                   threads = args.threads,
                   benchmark = args.benchmark) as tx_eff_hgvs:
        merged_transcripts = tx_eff_hgvs.get_updated_hgvs_transcripts(annovar_records)

    # Close the reference FASTA
    pysam.my_fasta.close()

    # Add additional annotations to each variant
    tx_eff_annotate = TxEffAnnotate(args.sequence_source)
    tx_eff_annotate.annotate(merged_transcripts)
    
    # If the cgd export is provided we can check if all the variants requested have transcripts in the output
    if args.variants:
        _check_output_for_all_variants(logger, json.load(args.variants), merged_transcripts)
    
    # Write the transcript effects to file
    TxEffWriter(args.out).write(merged_transcripts)
    
    print(f"Wrote {len(merged_transcripts)} transcripts to {args.out.name}")

    stop_time_real = time.perf_counter()
    stop_time_user = time.process_time()
    real_time, user_time = _get_time(stop_time_real - start_time_real, stop_time_user - start_time_user) 
    
    print(f"Time: real {real_time}, user {user_time}")
    logger.info(f"Time: real {real_time}, user {user_time}")

if __name__ == '__main__':
    _main()