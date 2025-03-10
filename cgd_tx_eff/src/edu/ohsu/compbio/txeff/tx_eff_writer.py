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
                "aminoAcidPosition": obj.hgvs_amino_acid_position,
                "basePosition": obj.hgvs_base_position,
                "exon": obj.exon,
                "gene": obj.hgnc_gene,
                "cDot": obj.hgvs_c_dot,
                "pDot1": obj.hgvs_p_dot_one,
                "pDot3": obj.hgvs_p_dot_three,
                "splicing": obj.splicing,
                "cdnaTranscript": obj.refseq_transcript,
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
        # with open(self._out_file, 'w') as f:
            # json.dump(info, f, cls=VariantTranscriptEncoder, indent=2)
        json.dump(variant_transcripts, self._out_file, cls=VariantTranscriptEncoder, indent=2)
