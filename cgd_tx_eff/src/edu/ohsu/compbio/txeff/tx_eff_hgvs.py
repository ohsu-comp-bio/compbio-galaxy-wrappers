'''
Created on Apr 20, 2022

@author: pleyte
'''

import argparse
import csv
from collections import defaultdict
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import logging
import os
from hgvs.dataproviders.uta import UTABase
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript
from hgvs.exceptions import HGVSInvalidVariantError, HGVSUsageError

VERSION = '0.0.1'

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

stream_handler = logging.StreamHandler()
logging_format = '%(levelname)s: [%(filename)s:%(lineno)s - %(funcName)s()]: %(message)s'

stream_format = logging.Formatter(logging_format)
stream_handler.setFormatter(stream_format)
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)

CHROM_MAP = {'1': 'NC_000001.10', '2': 'NC_000002.11', '3': 'NC_000003.11', '4': 'NC_000004.11', '5': 'NC_000005.9',
             '6': 'NC_000006.11', '7': 'NC_000007.14', '8': 'NC_000008.10', '9': 'NC_000009.11', '10': 'NC_000010.10',
             '11': 'NC_000011.9', '12': 'NC_000012.11', '13': 'NC_000013.10', '14': 'NC_000014.8', '15': 'NC_000015.9',
             '16': 'NC_000016.9', '17': 'NC_000017.10', '18': 'NC_000018.9', '19': 'NC_000019.9', '20': 'NC_000020.10',
             '21': 'NC_000021.8', '22': 'NC_000022.10', 'X': 'NC_000023.10', 'Y': 'NC_000024.9', 'MT': 'NC_012920.1'}

class Variant(object):
    def __init__(self, chromosome, position, ref, alt):
        self.chromosome = chromosome
        self.position = position
        self.reference = ref
        self.alt = alt

    def __str__(self):
        '''
        String representation of the AnnovarVariantFunction object
        '''
        return f'[Variant: {self.chromosome}-{self.position}-{self.reference}-{self.alt}]'
    
    def __eq__(self, other):
        return self.chromosome == other.chromosome  and self.position == other.position and self.reference == other.reference and self.alt == other.alt
    
    def __hash__(self):  
        return hash((self.chromosome, self.position, self.reference, self.alt))

def _noneIfEmpty(value: str):
    '''
    Return None if the string is an empty string.
    ''' 
    if value == '':
        return None
    return value


def _read_annovar_transcripts(input_filename):
    '''
    '''
    transcripts = list()

    with open(input_filename) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            transcript = VariantTranscript(row['chromosome'], row['position'], row['reference'], row['alt'])
            transcript.variant_effect = row['variant_effect']
            transcript.variant_type = row['variant_type']
            transcript.hgvs_amino_acid_position = _noneIfEmpty(row['hgvs_amino_acid_position'])
            transcript.hgvs_base_position = _noneIfEmpty(row['hgvs_base_position'])
            transcript.exon = row['exon']
            transcript.hgnc_gene = row['hgnc_gene']
            transcript.hgvs_c_dot = row['hgvs_c_dot']
            transcript.hgvs_p_dot_one = row['hgvs_p_dot_one']
            transcript.hgvs_p_dot_three = row['hgvs_p_dot_three']
            transcript.splicing = row['splicing']
            transcript.refseq_transcript = row['refseq_transcript']
            transcripts.append(transcript)
    
    return transcripts

def _correct_indel_coords(pos, ref, alt):
    """
    Using a VCF position, create coords that are compatible with HGVS nomenclature.
    Since we are already determining at this stage whether the event is an ins or del, also
    include the ins or del strings in the result.
    substitution event -> ac:g.[pos][ref]>[alt]
    :return:
    """
    lref = len(ref)
    lalt = len(alt)
    if lref == 1 and lalt == 1:
        # Substitution case
        change = '>'.join([ref, alt])
        new_pos = str(pos) + change
        return new_pos
    elif lref == lalt:
        # Multi-nucleotide substitution case
        # NG_012232.1: g.12_13delinsTG
        new_start = str(pos)
        new_end = str(int(pos) + lref - 1)
        new_pos = '_'.join([new_start, new_end]) + 'delins' + alt
        return new_pos
    elif lref > lalt:
        # Deletion case
        shift = lref - lalt
        if shift == 1:
            new_pos = str(int(pos) + 1) + 'del'
            return new_pos
        else:
            new_start = str(int(pos) + 1)
            new_end = str(int(pos) + shift)
            new_pos = '_'.join([new_start, new_end]) + 'del'
            return new_pos
    elif lalt > lref:
        # Insertion case
        new_start = str(pos)
        new_end = str(int(pos) + 1)
        new_pos = '_'.join([new_start, new_end]) + 'ins' + alt[1:]
        return new_pos
    else:
        # OTHER case
        raise Exception(f"Change type not supported: {pos}:{ref}>{alt}")
        
def _lookup_hgvs_transcripts(variants: list):
    '''
    Return the HGVS transcripts associated with a list of variants  
    '''
    # Initialize the HGVS connection
    hdp = hgvs.dataproviders.uta.connect()
    am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method='splign')
    hgvs_parser = hgvs.parser.Parser()
    
    transcripts = list()
    for variant in variants:
        variant_transcripts = __lookup_hgvs_transcripts(hgvs_parser, hdp, am, variant)
        logger.info(f"HGVS found {len(variant_transcripts)} transcripts for {variant}")
        
        # jDebug: this constraint may be removed. I expect that HGVS will find a transcript whenever Annovar did, but if that is not
        # true then I want to know about it. 
        if len(variant_transcripts) == 0:
            logger.warning(f'HGVS could not find any transcripts for variant {variant} which has transcripts known to Annovar.')
        else:
            transcripts.extend(variant_transcripts)
    
    return transcripts
        
def __lookup_hgvs_transcripts(hgvs_parser: hgvs.parser.Parser, hdp: UTABase, am: hgvs.assemblymapper.AssemblyMapper, variant: Variant):
    '''
    Use HGVS to find the transcripts for a variant
    '''
    hgvs_chrom = CHROM_MAP.get(variant.chromosome)
    
    if hgvs_chrom == None:
        raise Exception("Unknown chromosome: {chrom}-{pos}-{ref}-{alt}")
    
    # Look up the variant using HGVS            
    pos_part = _correct_indel_coords(variant.position, variant.reference, variant.alt)
    new_hgvs = hgvs_chrom + ':g.' + pos_part

    # jDebug: consider catching HGVSParseError
    var_g = hgvs_parser.parse_hgvs_variant(new_hgvs)
    
    tx_list = hdp.get_tx_for_region(str(var_g.ac), 'splign', str(var_g.posedit.pos.start), str(var_g.posedit.pos.end))
    
    hgvs_transcripts = list()
    
    for hgvs_transcript in tx_list:
        try:
            # HGVS's dicts aren't labeled so you have to access everything by index.   
            var_c = am.g_to_c(var_g, str(hgvs_transcript[0]))
            var_p = am.c_to_p(var_c)
            
            # Convert the three letter amino acid seq to a one letter and remove the transcript: prefix. 
            # Remove 'transcript:' to be left with only p.
            var_p1 = var_p.format(conf={"p_3_letter": False}).replace(var_p.ac+':','')
            var_p3 = var_p.format(conf={"p_3_letter": True}).replace(var_p.ac+':','')
            c_dot = var_c.type +'.' + str(var_c.posedit)

            variant_transcript = VariantTranscript(variant.chromosome, variant.position, variant.reference, variant.alt)
            
            if var_p3 == 'p.?':
                logger.debug(f"HGVS variant has ambiguous 'p.?' for {variant.chromosome}-{variant.position}-{variant.reference}-{variant.alt}, {var_c.ac}, {var_p.ac}")
                variant_type = None
            elif isinstance(var_p.posedit, hgvs.edit.AARefAlt):
                # jDebug: need to come back here and decide what to do with these.  All seem to be 'delins'
                # This type is not tested by John's VCF. Need to use a different VCF to get here. 
                variant_type = var_p.posedit.type
                logger.debug(f"HGVS variant does not have a position: ref={var_p.posedit.ref}, alt={var_p.posedit.alt}, type={var_p.posedit.type}, str={str(var_p.posedit)}. Keeping.")
            else:
                variant_transcript.hgvs_amino_acid_position = var_p.posedit.pos.start.pos
                variant_type = None
            
            variant_transcript.hgvs_base_position = var_c.posedit.pos.start.base
            variant_transcript.hgvs_c_dot = c_dot
            variant_transcript.hgvs_p_dot_one = var_p1
            variant_transcript.hgvs_p_dot_three = var_p3
            variant_transcript.refseq_transcript = var_c.ac
            variant_transcript.variant_type = variant_type
            variant_transcript.protein_transcript = var_p.ac
            
            hgvs_transcripts.append(variant_transcript)
        
        except HGVSUsageError as e:
            if("non-coding transcript" in e.args[0]):
                logger.info("Transcript probably non-coding, ignore: %s", str(e))
            else:
                raise(e)
        except HGVSInvalidVariantError as e:            
            logger.warn("Invalid variant: %s", str(e))
            raise(e)
    
    return hgvs_transcripts

def __get_unmatched_annovar_transcripts(annovar_dict: defaultdict(list), hgvs_dict: defaultdict(list)):
    '''
    After HGVS and Annovar transcripts have been merged, this function can be called to find Annovar transcripts that were not paired 
    with an HGVS trancsript. 
    '''
    transcripts = list()
    
    for (transcript_key, annovar_transcript_list) in annovar_dict.items():
        assert len(annovar_transcript_list) == 1, f'Annovar transcript map has more than one transcript matching {transcript_key}'
        
        # Check the HGVS dictionary for a key matching the annovar key. If there is a match, then the annovar transcript has already  
        # been processed. If HGVS does not have the key then the annovar transcript has not been looked at.
        if hgvs_dict.get(transcript_key) == None:
            logger.debug(f"Adding unmatched Annovar transcript {transcript_key}")
            transcripts.extend(annovar_transcript_list)
    
    return transcripts

    
def _merge_annovar_with_hgvs(annovar_transcripts: list, hgvs_transcripts: list, require_match=False):
    '''
    Given a list of transcripts from Annovar and HGVS, find those with the same genotype and transcript, and merge them into a single record.
    
    The optional ``require_match`` parameter can be used indicate that only transcripts which are found in both the Annovar 
    and HGVS lists, are returned.       
    '''
    transcripts = list()
    
    # Collect annovar records into a map keyed by genotype and transcript 
    annovar_dict = defaultdict(list)    
    key_maker = lambda x: "-".join([x.chromosome, x.position, x.reference, x.alt, x.refseq_transcript])
    for annovar_rec in annovar_transcripts:
        if annovar_rec.splicing == 'splicing':
            # jDebug: i think i have this wrong. I think HGVS can be merged with the splice variant. It is the intronic variant that HGVS doesn't have. 
            # Place splice variants directly in the final result because they won't be merged with an HGVS transcript since HGVS doesn't 
            # have splicing information.   
            transcripts.append(annovar_rec)
        else:
            annovar_dict[key_maker(annovar_rec)].append(annovar_rec)
    
    # Collect hgvs records into a map keyed by genotype and transcript    
    hgvs_dict = defaultdict(list)
    for hgvs_rec in hgvs_transcripts:
        hgvs_dict[key_maker(hgvs_rec)].append(hgvs_rec)

    # Iterate over every HGVS variant-transcript and see if there is a matching Annovar transcript
    for (transcript_key, hgvs_transcript_list) in hgvs_dict.items():
        # The dictionary stores the HGVS transcript in a list, but there will only be one HGVS transcript in the list. 
        assert len(hgvs_transcript_list) == 1
        hgvs_transcript = hgvs_transcript_list[0]
        
        annovar_matches = annovar_dict.get(transcript_key)
        if not annovar_matches:
            logger.debug(f"HGVS {transcript_key} does not match any Annovar transcripts")
            if require_match == False:
                transcripts.append(hgvs_transcript)
        elif len(annovar_matches) == 1:
            logger.debug(f"Merging matched HGVS and Annovar transcripts having key {transcript_key}")
            transcripts.append(_merge(hgvs_transcript, annovar_matches[0]))
        else:
            # jDebug: i think i have this wrong:
            # Annovar should only have one instance of each variant-transcript tuple. The only time you might have two Annovar transcripts is if 
            # one was, say, intronic, and the other is splicing. But splicing variants go straight into the ``transcripts`` list rather than
            # being added to the ``annovar_dict``.   
            raise Exception(f'Multiple Annovar transcripts found matching key {transcript_key}')

    # Add unmatched Annovar transcripts to the final list
    if require_match == False:
        transcripts.extend(__get_unmatched_annovar_transcripts(annovar_dict, hgvs_dict))

    return transcripts


def _merge(hgvs_transcript: VariantTranscript, annovar_transcript: VariantTranscript):
    '''
    Combine Annovar and HGVS data relating to the same transcript into a single record.  
    '''
    
    new_transcript = VariantTranscript(hgvs_transcript.chromosome, hgvs_transcript.position, hgvs_transcript.reference, hgvs_transcript.alt)
    
    assert hgvs_transcript.chromosome == annovar_transcript.chromosome
    assert hgvs_transcript.position == annovar_transcript.position
    assert hgvs_transcript.reference == annovar_transcript.reference
    assert hgvs_transcript.alt == annovar_transcript.alt
    
    if annovar_transcript.splicing == 'splicing':
        logger.error("jDebug: merging annovar splice variant")
        
    # Amino acid position can be empty, and will match between annovar and HGVS
    # jDebug: assert str(hgvs_transcript.hgvs_amino_acid_position) == str(annovar_transcript.hgvs_amino_acid_position), f"{hgvs_transcript.hgvs_amino_acid_position} != {annovar_transcript.hgvs_amino_acid_position}" 
    new_transcript.hgvs_amino_acid_position = hgvs_transcript.hgvs_amino_acid_position
    if str(hgvs_transcript.hgvs_amino_acid_position) != str(annovar_transcript.hgvs_amino_acid_position):
        logger.warning(f"HGVS and Annovar do not agree on amino acid position: {hgvs_transcript.hgvs_amino_acid_position} != {annovar_transcript.hgvs_amino_acid_position}")
    
    # base position may not match what HGVS says can be empty and will match between annovar and hgvs    
    new_transcript.hgvs_base_position = hgvs_transcript.hgvs_base_position
    if str(hgvs_transcript.hgvs_base_position) != str(annovar_transcript.hgvs_base_position):
        logger.debug(f"Using HGVS base_position rather than Annovar:  {hgvs_transcript.hgvs_base_position} != {annovar_transcript.hgvs_base_position}")
    
    # Only Annovar provides exon 
    assert hgvs_transcript.exon == None
    new_transcript.exon = annovar_transcript.exon
    
    # Only Annovar provides gene
    assert annovar_transcript.hgnc_gene != None
    assert hgvs_transcript.hgnc_gene == None
    new_transcript.hgnc_gene = annovar_transcript.hgnc_gene
    
    # Use HGVS's c. because Annovar's is not always correct
    assert hgvs_transcript.hgvs_c_dot != None
    new_transcript.hgvs_c_dot = hgvs_transcript.hgvs_c_dot
    
    if hgvs_transcript.hgvs_c_dot != annovar_transcript.hgvs_c_dot:
        logger.debug(f"Using HGVS c_dot rather than Annovar: {hgvs_transcript.hgvs_c_dot} != {annovar_transcript.hgvs_c_dot} ")
    
    # Use HGVS's p. because Annovar's is not always correct 
    assert hgvs_transcript.hgvs_p_dot_one != None
    new_transcript.hgvs_p_dot_one = hgvs_transcript.hgvs_p_dot_one
    
    if hgvs_transcript.hgvs_p_dot_one != annovar_transcript.hgvs_p_dot_one:
        logger.debug(f"Using HGVS p_dot (one letter) rather than Annovar: {hgvs_transcript.hgvs_p_dot_one} != {annovar_transcript.hgvs_p_dot_one}")

    assert hgvs_transcript.hgvs_p_dot_three != None
    new_transcript.hgvs_p_dot_three = hgvs_transcript.hgvs_p_dot_three
    
    if hgvs_transcript.hgvs_p_dot_three != annovar_transcript.hgvs_p_dot_three:
        logger.debug(f"Using HGVS p_dot (three letter) rather than Annovar: {hgvs_transcript.hgvs_p_dot_three} != {annovar_transcript.hgvs_p_dot_three}")
        
    # jDebug: i don't think HGVS ever tells us about splice variants
    # jDebug: special splice logic probalby needs to be moved to the top of this function since above validation will fail since splcie will have empty values
    assert hgvs_transcript.splicing == None
    new_transcript.splicing = annovar_transcript.splicing
    
    # jDebug: this will fail if we are accepting transcripts that are not matched in Annovar and HGVS. It is common for HGVS to provide multiple versions of the same
    #         transcript, while Annovar only provides the older version. 
    assert hgvs_transcript.refseq_transcript == annovar_transcript.refseq_transcript
    new_transcript.refseq_transcript = annovar_transcript.refseq_transcript
        
    # Variant effect is only provided by Annovar
    assert hgvs_transcript.variant_effect == None
    assert annovar_transcript.variant_effect != None
    new_transcript.variant_effect = annovar_transcript.variant_effect
    
    # Variant type is only provided by Annovar
    assert hgvs_transcript.variant_type == None
    assert annovar_transcript.variant_type != None
    new_transcript.variant_type = annovar_transcript.variant_type 
    
    # Only HGVS provides protein transcript
    # jDebug: I think this is only going to be true when it is not a splice variant 
    assert hgvs_transcript.protein_transcript != None
    new_transcript.protein_transcript = hgvs_transcript.protein_transcript
    
    return new_transcript


def get_summary(args, annovar_transcripts, annovar_variants, hgvs_transcripts, merged_transcripts):
    '''
    '''
    results = dict()
    results['arg_require_match'] = args.require_match
    results['annovar_transcript_count'] = len(annovar_transcripts)
    results['annovar_distinct_variant_count'] = len(annovar_variants)
    results['hgvs_transcript_count'] = len(hgvs_transcripts)
    
    results['hgvs_distinct_variant_count'] = len(set(map(lambda x: Variant(x.chromosome, x.position, x.reference, x.alt), hgvs_transcripts)))
    
    # Collect annovar records into a map keyed by genotype and transcript
    key_maker = lambda x: "-".join([x.chromosome, x.position, x.reference, x.alt, x.refseq_transcript])
    annovar_dict = defaultdict(list)
     
    annovar_splice_variant_transcript_count = 0
    
    for annovar_rec in annovar_transcripts:
        if annovar_rec.splicing == 'splicing':
            annovar_splice_variant_transcript_count += 1
        else:
            annovar_dict[key_maker(annovar_rec)].append(annovar_rec)
    
    results['annovar_splice_variant_transcript_count'] = annovar_splice_variant_transcript_count
    
    # Collect hgvs records into a map keyed by genotype and transcript    
    hgvs_dict = defaultdict(list)
    for hgvs_rec in hgvs_transcripts:
        hgvs_dict[key_maker(hgvs_rec)].append(hgvs_rec)
        
    matched_annovar_and_hgvs_transcript_count = 0
    unmatched_annovar_transcript_count = 0
    unmatched_hgvs_transcript_count = 0
    
    for transcript in merged_transcripts:
        transcript_key = key_maker(transcript)
        if annovar_dict.get(transcript_key) and hgvs_dict.get(transcript_key):
            matched_annovar_and_hgvs_transcript_count += 1
        elif annovar_dict.get(transcript_key) and not hgvs_dict.get(transcript_key):
            unmatched_annovar_transcript_count += 1
        elif hgvs_dict.get(transcript_key) and not annovar_dict.get(transcript_key):
            unmatched_hgvs_transcript_count += 1
        else:
            raise Exception("Unexpected summary condition")
        
    results['matched_annovar_and_hgvs_transcript_count'] = matched_annovar_and_hgvs_transcript_count
    results['unmatched_annovar_transcript_count'] = unmatched_annovar_transcript_count
    results['unmatched_hgvs_transcript_count'] = unmatched_hgvs_transcript_count
    
    results['merged_transcript_count'] = len(merged_transcripts)
    
    merged_distinct_variant_count = len(set(map(lambda x: Variant(x.chromosome, x.position, x.reference, x.alt), merged_transcripts)))
    results['merged_distinct_variant_count'] = merged_distinct_variant_count 
    
    # Essential sanity check - make sure counts add up
    sanity_check_total_transcripts = results['merged_transcript_count'] == results['matched_annovar_and_hgvs_transcript_count'] + results['unmatched_annovar_transcript_count'] + results['unmatched_hgvs_transcript_count']
    assert sanity_check_total_transcripts, 'Sanity check failed: Total number of transcripts does not equal sum of matched, and unmatched.'
         
    # Non-critical sanity check: all variants from annovar have at least one transcript in the final output
    sanity_check_variant_coverage = merged_distinct_variant_count == results['annovar_distinct_variant_count']
    if not sanity_check_variant_coverage:
        logger.warning(f"Not all of the Annovar variants made it into the list of variant-transcripts ({merged_distinct_variant_count}/{results['annovar_distinct_variant_count']})")

    sanity_check = sanity_check_total_transcripts and sanity_check_variant_coverage
    results['sanity_check'] = sanity_check

    return results


def _log_summary(results: dict):
    '''
    ''' 
    
    logger.info(f"Argument require_match: {results['arg_require_match']}")
    logger.info(f"Number of annovar transcripts: {results['annovar_transcript_count']}")
    logger.info(f"Number of Annovar transcripts that are splice variants: {results['annovar_splice_variant_transcript_count']}")
    logger.info(f"Number of distinct variants from Annovar: {results['annovar_distinct_variant_count']}")    
    logger.info(f"Number of HGVS transcripts: {results['hgvs_transcript_count']}")
    logger.info(f"Number of distinct variants from HGVS: {results['hgvs_distinct_variant_count']}")
    logger.info(f"Number of transcripts matched in Annovar and HGVS: {results['matched_annovar_and_hgvs_transcript_count']}")
    logger.info(f"Number of Annovar transcripts not matched with HGVS transcripts: {results['unmatched_annovar_transcript_count']}")
    logger.info(f"Number of HGVS transcripts not matched with Annovar transcripts: {results['unmatched_hgvs_transcript_count']}")
    logger.info(f"Number of distinct variants in merged transcript list: {results['merged_distinct_variant_count']}")
    logger.info(f"Total number of transcripts in final list: {results['merged_transcript_count']}")
    logger.info(f"Sanity check: {'Passed' if results['sanity_check'] else 'Failed' }")
    

def _write_csv(file_name, records: list):
    '''
    Write the transcript records to csv file
    ''' 
    fields = ['chromosome', 'position', 'reference', 'alt', 
              'variant_effect', 'variant_type', 'hgvs_amino_acid_position', 'hgvs_base_position', 
              'exon', 'hgnc_gene', 'hgvs_c_dot', 'hgvs_p_dot_one', 'hgvs_p_dot_three', 
              'splicing', 'refseq_transcript', 'protein_transcript']
    
    with open(file_name, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(fields)
        
        for rec in records:
            csv_writer.writerow([rec.chromosome, rec.position, rec.reference, rec.alt, 
                                rec.variant_effect, rec.variant_type, rec.hgvs_amino_acid_position, rec.hgvs_base_position,
                                rec.exon, rec.hgnc_gene, rec.hgvs_c_dot, rec.hgvs_p_dot_one, rec.hgvs_p_dot_three,
                                rec.splicing, rec.refseq_transcript, rec.protein_transcript])


def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Read files generated by Annovar and write out a CSV file with variant transcript effects.')

    parser.add_argument('-i', '--in_file',  
                        help='Input CSV (generated by tx_eff_annovar.py)',
                        type=argparse.FileType('r'),
                        required=True)
        
    parser.add_argument('-o', '--out_file', 
                        help='Output CSV', 
                        type=argparse.FileType('w'), 
                        required=True)
    
    parser.add_argument('-r', '--require_match',
                        help='Keep only transcripts that are found in both HGVS and Annovar',
                        action='store_true')
    
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    
    args = parser.parse_args()
        
    if os.environ.get('HGVS_SEQREPO_DIR') == None:
        logger.warning("The HGVS_SEQREPO_DIR environment variable is not defined. The remote seqrepo database will be used.")
    
    if os.environ.get('UTA_DB_URL') == None:
        logger.warning("The UTA_DB_URL environment variable is not defined. The remote UTA database will be used.")
        

    return args


def _main():
    '''
    main function
    '''
    args = _parse_args()
    
    annovar_transcripts = _read_annovar_transcripts(args.in_file.name)
    logger.debug(f'Read {len(annovar_transcripts)} Annovar transcripts from {args.in_file.name}')
    
    variants = set(map(lambda x: Variant(x.chromosome, x.position, x.reference, x.alt), annovar_transcripts))
    logger.debug(f'{len(variants)} distinct variants')
    
    hgvs_transcripts = _lookup_hgvs_transcripts(variants)
    logger.debug(f'{len(hgvs_transcripts)} transcripts from HGVS')
    
    merged_transcripts = _merge_annovar_with_hgvs(annovar_transcripts, hgvs_transcripts, require_match = args.require_match)
    
    _log_summary(get_summary(args, annovar_transcripts, variants, hgvs_transcripts, merged_transcripts))
    
    logger.info(f"Writing {args.out_file.name}")
    _write_csv(args.out_file.name, merged_transcripts)


if __name__ == '__main__':
    _main()