'''
Created on Apr 14, 2022

@author: pleyte
'''

from collections import defaultdict
import csv
from enum import Enum
import logging
import os
import re

from hgvs.location import BaseOffsetPosition
import hgvs.parser


def is_annovar_splicing_type(value: str):
    '''
    Return true if the string value matches one of variant types that Annovar uses for splicing variants.  
    '''
    return value in ['splicing', 'ncRNA_splicing']


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

        self._id = None

        self.chromosome = str(chromosome)
        self.position = int(position)
        self.reference = ref
        self.alt = alt
        
        self.variant_effect = None
        self.variant_type = None

        self.amino_acid_position = None 
        self.base_position = None 
        self.exon = None 
        self.gene = None 
        self.c_dot = None 
        self.p_dot1 = None 
        self.p_dot3 = None 
        self.splicing = None 
        self.cdna_transcript = None

    def __str__(self):
        '''
        String representation of the AnnovarVariantFunction object
        '''
        return f'[AnnovarVariantFunction: genotype={self.chromosome}-{self.position}-{self.reference}-{self.alt}, transcript={self.cdna_transcript}, variant_effect={self.variant_effect}, variant_type={self.variant_type}, aap={self.amino_acid_position}, bpos={self.base_position}, exon={self.exon}, gene={self.gene}, c.={self.c_dot}, p1.={self.p_dot1}, p3.={self.p_dot3}, splicing={self.splicing}, _id={self._id}]'
            
    def __eq__(self, obj):
        if obj is None:
            return False
        
        return self.chromosome == obj.chromosome \
            and self.position == obj.position \
            and self.reference == obj.reference \
            and self.alt == obj.alt \
            and self.variant_effect == obj.variant_effect \
            and self.variant_type == obj.variant_type \
            and self.amino_acid_position == obj.amino_acid_position \
            and self.base_position == obj.base_position \
            and self.exon == obj.exon \
            and self.gene == obj.gene \
            and self.c_dot == obj.c_dot \
            and self.p_dot1 == obj.p_dot1 \
            and self.p_dot3 == obj.p_dot3 \
            and self.splicing == obj.splicing \
            and self.cdna_transcript == obj.cdna_transcript
    
    def get_label(self):
        """
        Return this variant transcript as a string.  
        """
        return "-".join([self.chromosome, str(self.position), self.reference, self.alt, self.cdna_transcript])
    
    
class AnnovarParser(object):    
    '''
    Parses Annovar output
    ''' 
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.hgvs_parser = hgvs.parser.Parser()            

    def _unpack_evf_transcript_tuple(self, delimited_transcript: str):
        '''
        Takes the transcript tuple (index 2 in the tsv) from an exonic_variant_function file and parses out the values.  
        '''
        transcript_parts = delimited_transcript.split(':')

        gene = transcript_parts[0]
        cdna_transcript = transcript_parts[1]
        raw_exon = transcript_parts[2]

        if raw_exon == 'wholegene':
            self.logger.debug(f"This transcript's raw_exon value is 'wholegene' which means it is a Start loss: {transcript_parts}")
            amino_acid_position = None 
            hgvs_basep = None
            exon = raw_exon
            c_dot = None
            hgvs_p = None
            hgvs_three = None
        else:
            # Remove the prefix 'exon' prefix from the raw exon string like "exon4"
            exon = int(raw_exon.replace('exon',''))
            
            c_dot = transcript_parts[3]
            full_c = ':'.join([cdna_transcript, c_dot])
            
            # This may thrown an exception of type hgvs.exceptions.HGVSParseError but we don't catch it because 
            # it isn't possible to recover. 
            hgvs_basep = self.hgvs_parser.parse(full_c).posedit.pos.start
            
            # Sometimes basepair position is an integer, but it can also be an offset (eg c.371-3C>T)
            if hgvs_basep and type(hgvs_basep) is BaseOffsetPosition:
                self.logger.debug(f"Converting base pair position to string because it is an offset: {full_c}")
                hgvs_basep = str(hgvs_basep)
            
            # the exonic file may have protein info but we ignore it because 1) we don't trust the format of the p.; 2) it never
            # includes a protein transcript.
            amino_acid_position = None
            hgvs_p = None
            hgvs_three = None

        return amino_acid_position, hgvs_basep, exon, gene, c_dot, hgvs_p, hgvs_three, cdna_transcript

    def _unpack_vf_transcript_tuple(self, delimited_transcript: str):
        '''
        Takes the transcript tuple (index 1 in the tsv) from a variant_function file and parses out the values.
        '''
        if not delimited_transcript.strip():
            raise ValueError("Transcript definition is empty")
        
        exon = None
        c_dot = None

        # There are four different ways that information can be packaged that may in and annovar_row[. 
        # The format can be are determined by the number of ':' separated values. 
        tuple_type = len(delimited_transcript.split(':'))
        
        if tuple_type == 1:
            # type I is just a transcript: "NM_000791.4"
            cdna_transcript = self._parse_transcript_tuple_type1(delimited_transcript)
        elif tuple_type == 2:
            # type II looks like (NM_000791.4:c.-417_-416insGCGCTGCGG)"
            cdna_transcript, c_dot = self._parse_transcript_tuple_type2(delimited_transcript)
        elif tuple_type == 3:
            # type III looks like "NM_001206844.2(NM_001206844.2:exon5:c.1135-4T>C)"
            cdna_transcript, exon, c_dot  = self._parse_transcript_tuple_type3(delimited_transcript)
        elif self._is_multi_exon_transcript_tuple(delimited_transcript):
            # type V defines the transcript at two or more different exons and 
            # looks like "NM_001077690.1(NM_001077690.1:exon1:c.60+1C>-,NM_001077690.1:exon2:c.61-1C>-)"
            cdna_transcript, exon, c_dot = self._parse_transcript_tuple_type5(delimited_transcript)
        else:
            self.logger.warning("Annovar tuple format not recognized: " + str(delimited_transcript))
            cdna_transcript, exon, c_dot = None, None, None

        # Strip the 'exon' prefix from 'exon5'
        if(exon):
            exon = int(exon.replace('exon',''))
        
        # Handle the special case of a UTR that doesn't actually have a c. (eg "NM_001324237.2(NM_001324237.2:exon2:UTR5)" or "...:r.spl)"
        if c_dot and (c_dot.endswith('UTR3') or c_dot.endswith('UTR5') or c_dot.endswith('r.spl')):            
            c_dot = None    
                
        # Extract c. and cdna base position 
        if c_dot:
            try:
                full_c = ':'.join([cdna_transcript, c_dot])
                hgvs_basep = self.hgvs_parser.parse(full_c).posedit.pos.start
        
                # Sometimes basepair position is an integer, but it can also be an offset (eg c.371-3C>T)
                if hgvs_basep and type(hgvs_basep) is BaseOffsetPosition:
                    self.logger.debug(f"Converting base pair position to string because it is an offset: {full_c}")
                    # base position is a string because offsets are like "1135-4"
                    hgvs_basep = str(hgvs_basep)
        
            except hgvs.exceptions.HGVSParseError as e:
                # This doesn't matter because we use the c. from HGVS/UTA not this one from Annovar
                hgvs_basep = None
                c_dot_alt = re.search(r"c\..*\>([\w]+)", c_dot)            
                if c_dot_alt and len(c_dot_alt.group(1)) > 1:
                    self.logger.debug(f"HGVS parser failed to determine base position because the parser expects the c. alt to be just one base, but is {len(c_dot_alt.group(1))}: {c_dot}")
                else:
                    self.logger.debug(f"HGVS parser failed to determine base position from {full_c}: {e}")
        else:
            hgvs_basep = None

        return cdna_transcript, exon, c_dot, hgvs_basep 

    def _is_multi_exon_transcript_tuple(self, unparsed):
        '''
        Return true if the transcript tuple is the fifth type (see _parse_transcript_tuple_type5 and _get_multi_exon_transcript_tuples).
        '''
        return len(self._get_multi_exon_transcript_tuples(unparsed)) >=2
        
    def _get_multi_exon_transcript_tuples(self, unparsed):
        '''
        Return all the transcript definitions from a Type V tuple. Type V transcript tuples define multiple exons for the same transcript. 
        This type looks like NM_x(a,b,...) where "a,b,..." are two or more three part colon delimited tuples.
        Example:
            - NM_001243965.1(NM_001243965.1:exon2:c.278+4C>A,NM_001243965.1:exon3:c.281+1C>A,NM_001243965.1:exon4:c.282-1C>A)
            - NM_001352417.1(NM_001352417.1:exon1:c.60+1C>-,NM_001352417.1:exon2:c.61-1C>-)
            - NM_001306173.2(NM_001306173.2:exon12:c.1365+1G>T,NM_001306173.2:exon12:UTR3)
        This function returns a list of colon separated tuples like ['NM_x:exonN:c.','NM_y:exonM:c.'] 
        '''
        p = re.compile('NM_[0-9]+\.[0-9]:exon[0-9]+:(?:r\.spl|c\.[-+0-9ACTG>]+|UTR3|UTR5)')
        return p.findall(unparsed)
    
    def _parse_exonic_variant_function_row(self, annovar_row: list):
        '''
        Takes a single row from an Annovar exonic_variant_function file and returns one or more AnnovarVariantFunction objects
        ''' 
        annovar_recs = []
        
        # genotype fields are in positions 3,4,6,7. Index 5 is end position.
        chrom = str(annovar_row[3]).replace('chr','')
        pos = annovar_row[4]
        ref = annovar_row[6]
        alt = annovar_row[7]
        
        self.logger.debug(f"Parsing variant {chrom}-{pos}-{ref}-{alt}")
        
        if annovar_row[2] == 'UNKNOWN':
            # Don't bother with unknown types
            self.logger.debug("Variant is UNKNOWN type.")
            return []

        # The exonic_variant_function file provides variant effect ("VFX in the old code)
        variant_effect = annovar_row[1]
        variant_type = None
        
        external_id = self._get_external_id(annovar_row)
        
        transcript_tuples = annovar_row[2].rstrip(',').split(',')
        self.logger.debug(f"row has {len(transcript_tuples)} transcript tuples")

        for transcript_tuple in transcript_tuples:
            avf = AnnovarVariantFunction(chrom, pos, ref, alt)

            avf._id = external_id
            
            avf.variant_effect = variant_effect
            avf.variant_type = variant_type

            avf.amino_acid_position, \
            avf.base_position,       \
            avf.exon,                     \
            avf.gene,                \
            avf.c_dot,               \
            avf.p_dot1,           \
            avf.p_dot3,         \
            avf.cdna_transcript = self._unpack_evf_transcript_tuple(transcript_tuple)

            # Annovar sets the exon to "wholegene" to indicate a start loss. The 
            # annotation "startloss" isn't a term Annovar uses, it is our own custom 
            # type and it replaces whatever the current value is (eg 'frameshift deletion').  
            if avf.exon == "wholegene":
                self.logger.info(f"Start loss (aka wholegene) detected, substituting variant_effect overwriting '{avf.variant_effect}' with 'startloss': {avf}")
                avf.exon = None                
                avf.variant_effect = 'startloss'
                
            self.logger.debug(f"Parsed exonic_variant_function with transcript: {avf}")
            annovar_recs.append(avf)

        return annovar_recs
    
    def _parse_variant_function_row(self, annovar_row: list):
        '''
        Takes a single row from an Annovar variant_function file and returns one or more AnnovarVariantFunction objects.
        '''  
        annovar_recs = []
        
        # genotype fields are in positions 2,3,5,6. Index 4 is end position.  
        chrom = str(annovar_row[2]).replace('chr','')
        pos = annovar_row[3]
        ref = annovar_row[5]
        alt = annovar_row[6]
        
        self.logger.debug(f"Parsing variant {chrom}-{pos}-{ref}-{alt}")
        
        # Handle special cases of splicing variants
        
        if is_annovar_splicing_type(annovar_row[0]):            
            variant_type = None
            variant_splicing = annovar_row[0]
        else:
            variant_type = annovar_row[0]
            variant_splicing = None
           
        # The second column has a list of transcripts and related values
        transcript_tuples = self._split_variant_function_transcript_tuples(annovar_row[1])
        
        if(len(transcript_tuples) == 0):
            raise Exception(f"Variant does not have any transcript tuples: {chrom}-{pos}-{ref}-{alt}")
            
        self.logger.debug(f"row has {len(transcript_tuples)} transcript tuples")
        
        external_id = self._get_external_id(annovar_row)
        
        for transcript_tuple in transcript_tuples:
            avf = AnnovarVariantFunction(chrom, pos, ref, alt)
            
            avf._id = external_id
            
            # Variant effect is not present in variant_function files 
            avf.variant_effect = None
            
            avf.variant_type = variant_type
            avf.splicing = variant_splicing
            
            avf.cdna_transcript,  \
            avf.exon,               \
            avf.c_dot,         \
            avf.base_position = self._unpack_vf_transcript_tuple(transcript_tuple)
            
            if avf.cdna_transcript == None:
                self.logger.debug(f"Ignoring transcript {transcript_tuple}. Probably because it was intergenic and didn't have any useful information")
                continue
            elif avf.cdna_transcript.startswith("NR"):
                self.logger.debug("Ignoring non-coding transcript " + avf.cdna_transcript)
                continue

            self.logger.debug(f"Parsed variant_function record with transcript: {avf}")
            annovar_recs.append(avf)
        
        return annovar_recs

    def _split_variant_function_transcript_tuples(self, unparsed):
        '''
        Take an annovar variant_function field that has delimited information about transcripts 
        and split it into chunks for each transcript. The resulting string need to be further 
        parsed by the _unpack_vf_transcript_tuple function.
        Examples: 
            NM_001005484.1(dist=28)
            NM_001005484.1
            NM_001005277.1(dist=168389),NR_125957.1(dist=25774) 
        '''
        p = re.compile('(?<=,)?([^\\(,]+(\\([^\\)]+\\))?)')
        matches = p.findall(unparsed)
        
        # The first item in the tuple is the broadest match
        return [i[0] for i in matches]

    def _parse_transcript_tuple_type5(self, unparsed):
        '''
        Return the transcript from a definition that has two or more exons. 
        Examples:
            - NM_001077690.1(NM_001077690.1:exon1:c.60+1C>-,NM_001077690.1:exon2:c.61-1C>-)
            - NM_001243965.1(NM_001243965.1:exon2:c.278+4C>A,NM_001243965.1:exon3:c.281+1C>A,NM_001243965.1:exon4:c.282-1C>A) 
        
        Users want to use the lowest exon number, which we assume will always be the one on the left. 
        ''' 
        transcript_tuples = self._get_multi_exon_transcript_tuples(unparsed)
        
        assert len(transcript_tuples) >= 2, "Unrecognised variant function delimited transcript definition"
        
        # Even though there are multiple definitions of this transcript, we can only accept one. So we take the first.
        cdna_transcript, exon, c_dot  = transcript_tuples[0].split(':')

        return cdna_transcript, exon, c_dot
    
    def _parse_transcript_tuple_type3(self, unparsed):
        '''
        Return the values of the transcript definition in a string like "NM_001206844.2(NM_001206844.2:exon5:c.1135-4T>C)"
        '''
        p = re.compile('(?<=\().*:.*:.*(?=\)$)')
        matches = p.findall(unparsed)
        assert len(matches) == 1, "Unrecognised variant function delimited transcript definition"
        
        cdna_transcript, exon, c_dot = matches[0].split(':')
        return cdna_transcript, exon, c_dot
    
    def _parse_transcript_tuple_type2(self, unparsed):
        '''
        Return the values of a transcript definition like "NM_000791.4(NM_000791.4:c.-417_-416insGCGCTGCGG)"
        ''' 
        p = re.compile('(?<=\().*:.*(?=\)$)')
        matches = p.findall(unparsed)
        
        assert len(matches) == 1, "Unrecognised variant function delimited transcript definition"
        
        cdna_transcript, c_dot = matches[0].split(':')
        return cdna_transcript, c_dot
    
    def _parse_transcript_tuple_type1(self, transcript):
        '''
        Return the transcript in a definition like "NM_000791.4"
        '''
        # Don't use the intergenic transcripts that are labeled with "dist=nnn"
        if transcript.find('(dist=') > 0:
            self.logger.debug(f'Found intergenic transcript that will be thrown out: {transcript}')
            return None
        else:
            return transcript
                
    def parse_file(self, annovar_file_type: AnnovarFileType, file_name:str, delimiter = '\t'):
        '''
        Read an Annovar file and return a list of variants 
        '''
        annovar_recs = []

        if (os.path.getsize(file_name) == 0):
            # Sometimes none of the variants are exonic so the exonic file is empty. 
            self.logger.info(f"Annovar file {annovar_file_type} is empty.")
            return annovar_recs
            
        with open(file_name, 'r') as annovar_file:
            reader = csv.reader(annovar_file, delimiter = delimiter)
            for row in reader:
                if annovar_file_type == AnnovarFileType.ExonicVariantFunction:
                    assert row[0].startswith('line')
                    if len(row) != 9:
                        self.logger.info(f"Exonic variant function file is expected have 9 columns. Found {len(row)}.")
                    
                    # Parse a row in an exonic_variant_function file
                    annovar_recs.extend(self._parse_exonic_variant_function_row(row))
                elif annovar_file_type == AnnovarFileType.VariantFunction:
                    if len(row) != 8:
                        self.logger.info(f"Variant function file is expected to have 17 columns. Found {len(row)}. Make sure you run the annotate_variation.pl with the '--otherinfo' parameter to include other information from the original VCF")
                    
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
        
        key_maker = lambda x: "-".join([x.chromosome, str(x.position), x.reference, x.alt, x.cdna_transcript])
        for annovar_rec in annovar_records:
            if annovar_rec.chromosome == None or annovar_rec.position == None or annovar_rec.reference == None or annovar_rec.alt == None or annovar_rec.cdna_transcript == None:
                self.logger.error(f"There is going to be a problem with {annovar_rec}")
                exit
            annovar_dict[key_maker(annovar_rec)].append(annovar_rec)

        # Merge the records that come from different files but refer to the same transcript.
        annovar_records = []
        
        # Each value in the dictionary is list of Annovar records with the same genotype and transcript
        for matching_genotypes in annovar_dict.values():
            left_rec = matching_genotypes[0]
            
            if len(matching_genotypes) == 1:
                # Nothing to merge with
                pass
            elif len(matching_genotypes) == 2:
                # There will be two records for a transcript:
                #  - The transcript is exonic so we get one record in annovar.variant_function and one in annovar.exonic_variant_function 
                #  - The transcript is intronic, and it is splicing so we get two records from annovar.variant_function
                right_rec = matching_genotypes[1]
                
                # Merged transcripts are from the same variants so the external id is expected to be the same as well. 
                if left_rec._id != right_rec._id:
                    raise ValueError(f"Transcripts on the same variant are expected to have the same id: left={left_rec}, right={right_rec}")
                
                self.logger.debug(f"Merging left={left_rec} and right={right_rec}")
                self._merge(left_rec, right_rec)
            elif len(matching_genotypes) == 3:
                # There will be three matches when the transcript is exonic and splicing: One record from annovar.exonic_variant_function 
                # and two records from annovar.variant_function where one is 'exonic' and the other is 'splicing'.
                right1_rec = matching_genotypes[1]
                right2_rec = matching_genotypes[2]  

                # Merged transcripts are from the same variants so the external id is expected to be the same as well. 
                if left_rec._id != right1_rec._id:
                    raise ValueError(f"Transcripts on the same variant are expected to have the same id: left={left_rec}, right1_rec={right1_rec}")
                elif left_rec._id != right2_rec._id:
                    raise ValueError(f"Transcripts on the same variant are expected to have the same id: left={left_rec}, right2_rec={right2_rec}")

                # The two records from annovar.variant_function will label one of the records as exonic and one as splicing
                if not (left_rec.variant_type == 'exonic' or right1_rec.variant_type == 'exonic' or right2_rec.variant_type == 'exonic'):
                    raise ValueError("Transcript has three records and none are exonic: " + left_rec.get_label())
                
                # If the annovar.exonic_variant_function file says the variant is exonic, but the annovar.variant_function file says the same
                # transcript is intronic the information we have is contradictory and we skip it.  Example: Annovar says 1-146466030-C-G 
                # is exonic and intronic. 
                if left_rec.variant_type == 'intronic' or right1_rec.variant_type == 'intronic' or right2_rec.variant_type == 'intronic':
                    self.logger.info(f"Ignoring because Annovar has conflicting variant types for {left_rec.get_label()}: {left_rec.variant_type}, {right1_rec.variant_type}, {right2_rec.variant_type}")
                    continue

                # The two records from annovar.variant_function will label one of the records as exonic and one as splicing
                if not (left_rec.splicing == 'splicing' or right1_rec.splicing == 'splicing' or right2_rec.splicing == 'splicing'):
                    raise ValueError("Transcript has three records and none are are splicing: " + left_rec.get_label())
                
                self.logger.debug(f"Merging three records left={left_rec} and right1={right1_rec}, right2={right2_rec}")
                self._merge(left_rec, right1_rec)
                self._merge(left_rec, right2_rec)
            else:
                # Our current workflow only involves two annovar input files and there will only be a maximum of three matching transcripts.
                # This exception ensures our expectation is true. 
                raise ValueError(f"This transcript has more than 3 annovar records, which is not supported: {left_rec}")
            
            annovar_records.append(left_rec)
        
        return annovar_records
            
    def _merge(self, left_rec: AnnovarVariantFunction, right_rec: AnnovarVariantFunction):
        '''
        Merge the right AnnovarVariantFunction record into the left one in order to fill in missing values 
        '''
        left_rec.variant_effect = (left_rec.variant_effect or right_rec.variant_effect)
        left_rec.variant_type = (left_rec.variant_type or right_rec.variant_type)
        left_rec.amino_acid_position = (left_rec.amino_acid_position or right_rec.amino_acid_position) 
        left_rec.base_position = (left_rec.base_position or right_rec.base_position) 
        left_rec.exon = (left_rec.exon or right_rec.exon) 
        left_rec.gene = (left_rec.gene or right_rec.gene) 
        left_rec.c_dot = (left_rec.c_dot or right_rec.c_dot) 
        left_rec.p_dot1 = (left_rec.p_dot1 or right_rec.p_dot1) 
        left_rec.p_dot3 = (left_rec.p_dot3 or right_rec.p_dot3) 
        left_rec.splicing = (left_rec.splicing or right_rec.splicing) 
        left_rec.cdna_transcript = (left_rec.cdna_transcript or right_rec.cdna_transcript)        

    def _get_external_id(self, annovar_line):
        """
        Checks the last field an annovar vf or annovar evf line for an "id:nnn" attribute and return the id if it is present or None if it is not.
        """
        match = re.search('id:(\d+)', annovar_line[-1])
        if match:
            return match.group(1)
        else:
            return None
        
        