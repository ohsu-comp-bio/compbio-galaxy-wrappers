'''
Created on Mar 4, 2025

@author: pleyte
'''
from json import JSONEncoder
import json

from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript

class VariantTranscriptEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, VariantTranscript):
            return {
                "id": obj._id,
                "chromosome": obj.chromosome, 
                "position": obj.position, 
                "reference": obj.reference, 
                "alt": obj.alt,
                "variantEffect": obj.variant_effect,
                "variantType": obj.variant_type,
                "aminoAcidPosition": obj.amino_acid_position,
                "basePosition": obj.base_position,
                "exon": obj.exon,
                "gene": obj.gene,
                "cDot": obj.c_dot,
                "pDot1": obj.p_dot1,
                "pDot3": obj.p_dot3,
                "splicing": obj.splicing,
                "cdnaTranscript": obj.cdna_transcript,
                "proteinTranscript": obj.protein_transcript,
                "sequenceVariant": obj.sequence_variant,
                "referenceContext": obj.reference_context
            }
        else:
            return super().default(obj)
        
class TxEffWriter(object):
    '''
    Writes variant transcripts to file 
    '''
    def __init__(self, out_file):
        '''
        Constructor
        '''
        self._out_file = out_file
    
    def write(self, variant_transcripts: list):
        '''
        Write transcripts to output file 
        '''
        json.dump(sorted(variant_transcripts), self._out_file, cls=VariantTranscriptEncoder, indent=2)
