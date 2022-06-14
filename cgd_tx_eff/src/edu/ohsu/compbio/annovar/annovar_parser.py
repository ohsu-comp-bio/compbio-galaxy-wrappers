'''
Created on Apr 14, 2022

@author: pleyte
'''

import csv
from enum import Enum
import hgvs.parser
import logging
from collections import defaultdict

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

stream_handler = logging.StreamHandler()
logging_format = '%(levelname)s: [%(filename)s:%(lineno)s - %(funcName)s()]: %(message)s'

stream_format = logging.Formatter(logging_format)
stream_handler.setFormatter(stream_format)
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)

class AnnovarFileType(Enum):
    '''
    Annovar file types supported by this parser
    '''
    ExonicVariantFunction = 1
    VariantFunction = 2
    
class AnnovarVariantFunction(object):
    '''
    Simple Annovar variant function object
    '''
    def __init__(self, chromosome, position, ref, alt):
        '''
        Constructor
        '''
        
        assert type(position) == int or position.isnumeric(), 'Position must be an integer' 
        assert alt.find(',') == -1, 'ALT may only have one value'
        
        self.chromosome = chromosome
        self.position = int(position)
        self.reference = ref
        self.alt = alt
        
        self.variant_effect = None
        self.variant_type = None

        self.hgvs_amino_acid_position = None 
        self.hgvs_base_position = None 
        self.exon = None 
        self.hgnc_gene = None 
        self.hgvs_c_dot = None 
        self.hgvs_p_dot_one = None 
        self.hgvs_p_dot_three = None 
        self.splicing = None 
        self.refseq_transcript = None

    def __str__(self):
        '''
        String representation of the AnnovarVariantFunction object
        '''
        return f'[AnnovarVariantFunction: genotype={self.chromosome}-{self.position}-{self.reference}-{self.alt}, transcript={self.refseq_transcript}, variant_effect={self.variant_effect}, variant_type={self.variant_type}, aap={self.hgvs_amino_acid_position}, bpos={self.hgvs_base_position}, exon={self.exon}, gene={self.hgnc_gene}, c.={self.hgvs_c_dot}, p1.={self.hgvs_p_dot_one}, p3.={self.hgvs_p_dot_three}, splicing={self.splicing}]'
            
    def __eq__(self, obj):
        return self.chromosome == obj.chromosome \
            and self.position == obj.position \
            and self.reference == obj.reference \
            and self.alt == obj.alt \
            and self.variant_effect == obj.variant_effect \
            and self.variant_type == obj.variant_type \
            and self.hgvs_amino_acid_position == obj.hgvs_amino_acid_position \
            and self.hgvs_base_position == obj.hgvs_base_position \
            and self.exon == obj.exon \
            and self.hgnc_gene == obj.hgnc_gene \
            and self.hgvs_c_dot == obj.hgvs_c_dot \
            and self.hgvs_p_dot_one == obj.hgvs_p_dot_one \
            and self.hgvs_p_dot_three == obj.hgvs_p_dot_three \
            and self.splicing == obj.splicing \
            and self.refseq_transcript == obj.refseq_transcript
            
    
    
class AnnovarParser(object):    
    '''
    Parses Annovar output
    ''' 
    def __init__(self):
        self.hgvs_parser = hgvs.parser.Parser()
    
    def _get_file_type(self, file_name:str):
        '''
        Return the Annovar output type as determined by file extension
        '''
        if file_name.endswith('.exonic_variant_function'):
            return AnnovarFileType.ExonicVariantFunction
        elif file_name.endswith('.variant_function'):
            return AnnovarFileType.VariantFunction
        else:
            raise Exception(f"Unrecognized Annovar file type: {file_name}")
        
    def _unpack_evf_transcript_tuple(self, delimited_transcript: str):
        '''
        Takes the transcript tuple (index 2 in the tsv) from an exonic_variant_function file and parses out the values.  
        '''
        transcript_parts = delimited_transcript.split(':')

        hgnc_gene = transcript_parts[0]
        refseq_transcript = transcript_parts[1]
        raw_exon = transcript_parts[2]

        if raw_exon != 'wholegene':
            # Remove the prefix 'exon' prefix from the raw exon string like "exon4"
            exon = raw_exon.replace('exon','')
            
            hgvs_c_dot = transcript_parts[3]
            full_c = ':'.join([refseq_transcript, hgvs_c_dot])
            
            # This may thrown an exception of type hgvs.exceptions.HGVSParseError but we don't catch it because 
            # it isn't possible to recover. 
            hgvs_basep = self.hgvs_parser.parse(full_c).posedit.pos.start
            
            # p. single amino acid
            hgvs_p = transcript_parts[4]
            full_p = ':'.join([refseq_transcript, hgvs_p])
            
            try:
                parsed_p_dot = self.hgvs_parser.parse(full_p)
                hgvs_three = 'p.' + str(parsed_p_dot.posedit)
                amino_acid_position = parsed_p_dot.posedit.pos.start.pos                
            except hgvs.exceptions.HGVSParseError as e:
                logger.warning(f"Unable to parse {full_p}: {e}")
                hgvs_three = None
                amino_acid_position = None
            

        return amino_acid_position, hgvs_basep, exon, hgnc_gene, hgvs_c_dot, hgvs_p, hgvs_three, refseq_transcript

    def _unpack_vf_transcript_tuple(self, delimited_transcript: str):
        '''
        Takes the transcript tuple (index 1 in the tsv) from a variant_function file and parses out the values.
        '''
        exon = None
        hgvs_c_dot = None

        # There are three different types of information that may be in annovar_row[1]. The types are determined by the number of ':' separated values. 
        tuple_type = len(delimited_transcript.split(':'))
        
        if tuple_type == 1:
            # Don't use transcripts that are labeled with "dist=nnn" 
            if "(dist=" in delimited_transcript:
                refseq_transcript = None
            else:
                refseq_transcript = delimited_transcript
        elif tuple_type == 2:
            # Looks like "NM_000791.4(NM_000791.4:c.-417_-416insGCGCTGCGG)"
            refseq_transcript = delimited_transcript.split('(')[0]
            hgvs_c_dot = delimited_transcript.split(':')[1].rstrip(')')
        elif tuple_type == 3:
            # Looks like "NM_001206844.2(NM_001206844.2:exon5:c.1135-4T>C)"
            refseq_transcript = delimited_transcript.split('(')[0]
            exon = delimited_transcript.split(':')[1][4:]
            hgvs_c_dot = delimited_transcript.split(':')[2].rstrip(')')
        else:
            raise Exception("Tuple format not recognized: " + str(delimited_transcript))
        
        if hgvs_c_dot:
            full_c = ':'.join([refseq_transcript, hgvs_c_dot])
            hgvs_basep = self.hgvs_parser.parse(full_c).posedit.pos.start
        else:
            hgvs_basep = None

        return hgvs_basep, exon, hgvs_c_dot, refseq_transcript

    def _parse_exonic_variant_function_row(self, annovar_row: list):
        '''
        Takes a single row from an Annovar exonic_variant_function file and returns one or more AnnovarVariantFunction objects
        ''' 
        annovar_recs = list()
        
        # genotype fields are in positions 3,4,6,7 and 8,9,11,12 and we want the second group because the first group
        # may have been altered by the convert2annovar.pl script.   
        chrom = str(annovar_row[8]).replace('chr','')
        pos = annovar_row[9]
        ref = annovar_row[11]
        alt = annovar_row[12]
        
        logger.debug(f"Parsing variant {chrom}-{pos}-{ref}-{alt}")
        
        if annovar_row[2] == 'UNKNOWN':
            # Don't bother with unknown types
            logger.debug("Variant is UNKNOWN type.")
            return []

        # The exonic_variant_function file provides variant effect ("VFX in the old code)
        variant_effect = annovar_row[1]
        variant_type = None
        
        transcript_tuples = annovar_row[2].rstrip(',').split(',')
        logger.debug(f"row has {len(transcript_tuples)} transcript tuples")

        for transcript_tuple in transcript_tuples:
            avf = AnnovarVariantFunction(chrom, pos, ref, alt)

            avf.variant_effect = variant_effect
            avf.variant_type = variant_type

            avf.hgvs_amino_acid_position, \
            avf.hgvs_base_position,       \
            avf.exon,                     \
            avf.hgnc_gene,                \
            avf.hgvs_c_dot,               \
            avf.hgvs_p_dot_one,           \
            avf.hgvs_p_dot_three,         \
            avf.refseq_transcript = self._unpack_evf_transcript_tuple(transcript_tuple)

            logger.debug(f"Parsed exonic_variant_function with transcript: {avf}")
            annovar_recs.append(avf)

        return annovar_recs
        

    def _parse_variant_function_row(self, annovar_row: list):
        '''
        Takes a single row from an Annovar variant_function file and returns one or more AnnovarVariantFunction objects.
        '''  
        annovar_recs = list()
        
        # genotype fields are in positions 2,3,5,6 and 7,8,10,11 and we want the second group because the first group
        # may have been altered by the convert2annovar.pl script.
        chrom = str(annovar_row[7]).replace('chr','')
        pos = annovar_row[8]
        ref = annovar_row[10]
        alt = annovar_row[11]
        
        logger.debug(f"Parsing variant {chrom}-{pos}-{ref}-{alt}")
        
        # Handle special cases of splicing variants
        if annovar_row[0] == 'splicing' or annovar_row[0] == 'ncRNA_splicing':            
            variant_type = None
            variant_splicing = annovar_row[0]
        else:
            variant_type = annovar_row[0]
            variant_splicing = None
           
        transcript_tuples = annovar_row[1].split(',')
        
        if(len(transcript_tuples) == 0):
            raise Exception(f"Variant does not have any transcript tuples: {chrom}-{pos}-{ref}-{alt}")
            
        logger.debug(f"row has {len(transcript_tuples)} transcript tuples")
        
        for transcript_tuple in transcript_tuples:
            avf = AnnovarVariantFunction(chrom, pos, ref, alt)
            
            # Variant effect is not present in variant_function files 
            avf.variant_effect = None
            
            avf.variant_type = variant_type
            avf.splicing = variant_splicing
            
            avf.hgvs_base_position, \
            avf.exon,               \
            avf.hgvs_c_dot,         \
            avf.refseq_transcript = self._unpack_vf_transcript_tuple(transcript_tuple)
            
            logger.debug(f"Parsed variant_function record with transcript: {avf}")
            annovar_recs.append(avf)
        
        return annovar_recs


    def parse_file(self, file_name:str, delimiter = '\t'):
        '''
        Read an Annovar file and return a list of variants 
        '''

        annovar_file_type = self._get_file_type(file_name)

        annovar_recs = list()

        with open(file_name, 'r') as annovar_file:
            reader = csv.reader(annovar_file, delimiter = delimiter)
            for row in reader:
                if annovar_file_type == AnnovarFileType.ExonicVariantFunction:
                    assert len(row) == 19, "Exonic variant function file must have 19 columns. Make sure you run annovar with the parameter to include other information from the original VCF"
                    assert row[0].startswith('line')

                    # Parse a row in an exonic_variant_function file
                    annovar_recs.extend(self._parse_exonic_variant_function_row(row))
                elif annovar_file_type == AnnovarFileType.VariantFunction:
                    assert len(row) == 18, "Variant function file must have 18 columns. Make sure you run annovar with the parameter to include other information from the original VCF"
                    
                    # Parse a row in a variant_function file
                    annovar_recs.extend(self._parse_variant_function_row(row))
                else: 
                    raise Exception("Unknown Annovar file type")

        return annovar_recs
    
    def merge(self, annovar_records: list):
        '''
        Take a list of AnnovarVariantFunction beans and combine the ones with the same genotype and transcript into a single record.
        '''
        
        # Collect annovar records into a map keyed by genotype and transcript 
        annovar_dict = defaultdict(list)
        
        key_maker = lambda x: "-".join([x.chromosome, str(x.position), x.reference, x.alt, x.refseq_transcript, (x.splicing or '')])
        for annovar_rec in annovar_records:
            annovar_dict[key_maker(annovar_rec)].append(annovar_rec)

        # Merge the records that come from different files but refer to the same transcript.
        annovar_records = list()
        
        # Each value in the dictionary is list of Annovar records with the same genotype and transcript
        for matching_genotypes in annovar_dict.values():
            left_rec = matching_genotypes[0]
            
            if left_rec.splicing == 'splicing':
                # Splicing variants don't get merged 
                pass
            elif len(matching_genotypes) > 2:
                # Our current workflow only involves two annovar input files so there will only be a maximum of two matching transcripts; this exception ensures
                # that is the case. Once we are through development and testing, this condition can be removed.
                  
                # This function can only merge two records. If you need to merge more than two than two then the code needs to be updated.                  
                raise Exception(f"This transcript has more than 2 instances, which is not supported; see code comment. rec={left_rec}")
            elif len(matching_genotypes) == 2:
                right_rec = matching_genotypes[1]
                logger.debug(f"Merging left={left_rec} and right={right_rec}")
                self._merge(left_rec, right_rec)
            
            annovar_records.append(left_rec)
        
        return annovar_records
            
    def _merge(self, left_rec: AnnovarVariantFunction, right_rec: AnnovarVariantFunction):
        '''
        Merge the right AnnovarVariantFunction record into the left one in order to fill in missing values 
        '''
        left_rec.variant_effect = (left_rec.variant_effect or right_rec.variant_effect)
        left_rec.variant_type = (left_rec.variant_type or right_rec.variant_type)
        left_rec.hgvs_amino_acid_position = (left_rec.hgvs_amino_acid_position or right_rec.hgvs_amino_acid_position) 
        left_rec.hgvs_base_position = (left_rec.hgvs_base_position or right_rec.hgvs_base_position) 
        left_rec.exon = (left_rec.exon or right_rec.exon) 
        left_rec.hgnc_gene = (left_rec.hgnc_gene or right_rec.hgnc_gene) 
        left_rec.hgvs_c_dot = (left_rec.hgvs_c_dot or right_rec.hgvs_c_dot) 
        left_rec.hgvs_p_dot_one = (left_rec.hgvs_p_dot_one or right_rec.hgvs_p_dot_one) 
        left_rec.hgvs_p_dot_three = (left_rec.hgvs_p_dot_three or right_rec.hgvs_p_dot_three) 
        left_rec.splicing = (left_rec.splicing or right_rec.splicing) 
        left_rec.refseq_transcript = (left_rec.refseq_transcript or right_rec.refseq_transcript)        
