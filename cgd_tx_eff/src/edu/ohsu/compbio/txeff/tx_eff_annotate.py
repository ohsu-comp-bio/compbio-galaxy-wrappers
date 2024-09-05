'''
Created on Aug. 12, 2024

@author: pleyte
'''

import logging
import os

from biocommons.seqrepo import SeqRepo

from edu.ohsu.compbio.txeff.util import chromosome_map
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript


class TxEffAnnotate(object):
    '''
    Add variant annotations to each variant transcript.  
    '''
    def __init__(self):
        '''
        Constructor
        '''
        # Set up a logger for this class
        self.logger = logging.getLogger(__name__)
        
        # Setup the SeqRepo utility that will be used for the Reference Context annotation
        if os.environ.get('HGVS_SEQREPO_DIR') is None:
            raise ValueError("The HGVS_SEQREPO_DIR environment variable is not defined")
        else:
            self.logger.info(f"Using SeqRepo {os.environ.get('HGVS_SEQREPO_DIR')}")
            self.sr = SeqRepo(os.environ.get('HGVS_SEQREPO_DIR'))
        
        self.REFERENCE_CONTEXT_LENGTH_PER_SIDE = 10
        
    def annotate(self, transcripts: list):
        '''
        Add annotations to each variant transcript
        '''
        self.logger.info("Applying Reference context annotations")
        for transcript in transcripts:
            transcript.reference_context = self._get_reference_context(transcript)

    def _get_reference_context(self, transcript: VariantTranscript):
        '''
        Look up the reference sequence that surrounds a variant. 
        '''
        refseq_chromosome = chromosome_map.get_refseq(transcript.chromosome)

        context_start = transcript.position - self.REFERENCE_CONTEXT_LENGTH_PER_SIDE - 1
        if context_start < 0:
            context_start = 0
            
        context_stop = transcript.position + len(transcript.reference) + self.REFERENCE_CONTEXT_LENGTH_PER_SIDE - 1
        
        return self.sr[refseq_chromosome][context_start:context_stop]
