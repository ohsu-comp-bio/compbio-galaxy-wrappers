'''
Created on Aug. 12, 2024

@author: pleyte
'''

import logging
import os
from pathlib import Path

from biocommons.seqrepo.dataproxy import SeqRepoRESTDataProxy
from biocommons.seqrepo.seqrepo import SeqRepo

from edu.ohsu.compbio.txeff.util import chromosome_map
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript


class TxEffAnnotate(object):
    '''
    Add variant annotations to each variant transcript.  
    '''
    def __init__(self, seqrepo = None):
        '''
        Constructor
        '''
        # Set up a logger for this class
        self.logger = logging.getLogger(__name__)
        self.REFERENCE_CONTEXT_LENGTH_PER_SIDE = 10
        
        self._configure_seqrepo(seqrepo)

    def _configure_seqrepo(self, seqrepo = None):
        """
        Configure the SeqRepo class to find sequences using the file repository or the SeqRepo REST service.
        """
        self.logger.debug("Determining SeqRepo repository location")

        env_seqrepo_url = os.environ.get('HGVS_SEQREPO_URL')
        env_seqrepo_directory = os.environ.get('HGVS_SEQREPO_DIR')
        
        if seqrepo and seqrepo.lower().startswith("http"):
            self._use_seqrepo_service(seqrepo)
        elif seqrepo and seqrepo.startswith("/"):
            self._use_seqrepo_file_repository(seqrepo)
        elif seqrepo:
            raise ValueError("The seqrepo parameter must be a url or a path not: " + seqrepo)
        elif env_seqrepo_url:
            if not env_seqrepo_url.lower().startswith('http'):
                self.logger.warning("The SeqRepo service is requested but the URL appears invalid: " + env_seqrepo_url)
            self._use_seqrepo_service(env_seqrepo_url)
        elif env_seqrepo_directory:
            if not Path(env_seqrepo_directory).is_dir():
                self.logger.warning("The SeqRepo file repository is requested but the path appears invalid: " + env_seqrepo_directory)
            self._use_seqrepo_file_repository(env_seqrepo_directory)
        else:
            raise ValueError("Unable to locate a SeqRepo instance. A url or path must be specified in the environment or as a parameter") 

    def _use_seqrepo_file_repository(self, path):
        """
        """
        self.logger.info(f"Using SeqRepo file repository at " + path)
        self.sr = SeqRepo(path)
    
    def _use_seqrepo_service(self, url):
        """
        """
        self.logger.info(f"Using SeqRepo service at " + url)
        self.sr = SeqRepoRESTDataProxy(url)
            
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
        
        # Use service or file repository depending on configuration 
        if type(self.sr) == SeqRepoRESTDataProxy:
            return self.sr.get_sequence(refseq_chromosome, context_start, context_stop)
        else:
            return self.sr[refseq_chromosome][context_start:context_stop]
