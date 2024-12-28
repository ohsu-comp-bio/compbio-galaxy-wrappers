from edu.ohsu.compbio.annovar.annovar_parser import AnnovarVariantFunction

class VariantTranscript(AnnovarVariantFunction):
    '''
    '''
    def __init__(self, chromosome, position, ref, alt):
        '''
        Create new VariantTranscript using chromosome, position, ref, and alt.
        '''
        super().__init__(chromosome, position, ref, alt)
            
        self.protein_transcript = None
        
        # An hgvs.sequencevariant.SequenceVariant with g. notation 
        self.sequence_variant = None
        
        # Reference sequence of 10 bases before and after the variant 
        self.reference_context = None

    def __eq__(self, obj):
        if obj is None:
            return False
        elif super.__eq__(self,obj):
            return self.protein_transcript == obj.protein_transcript
        else:           
            return False
    
    def __str__(self):
        '''
        String representation of the VariantTranscript object
        '''
        return f'[AnnovarVariantFunction: genotype={self.chromosome}-{self.position}-{self.reference}-{self.alt}-transcript={self.refseq_transcript}-{self.protein_transcript}' + ('-splicing' if self.splicing else '') 

    def get_copy(self):
        '''
        Return a copy of this trancript 
        ''' 
        transcript = VariantTranscript(self.chromosome, self.position, self.reference, self.alt)
        transcript.variant_effect = self._noneIfEmpty(self.variant_effect)
        transcript.variant_type = self._noneIfEmpty(self.variant_type)
        transcript.hgvs_amino_acid_position = self._noneIfEmpty(self.hgvs_amino_acid_position)
        transcript.hgvs_base_position = self._noneIfEmpty(self.hgvs_base_position)
        transcript.exon = self._noneIfEmpty(self.exon)
        transcript.hgnc_gene = self._noneIfEmpty(self.hgnc_gene)
        transcript.sequence_variant = self._noneIfEmpty(self.sequence_variant)
        transcript.hgvs_c_dot = self._noneIfEmpty(self.hgvs_c_dot)
        transcript.hgvs_p_dot_one = self._noneIfEmpty(self.hgvs_p_dot_one)
        transcript.hgvs_p_dot_three = self._noneIfEmpty(self.hgvs_p_dot_three)
        transcript.splicing = self._noneIfEmpty(self.splicing)
        transcript.refseq_transcript = self._noneIfEmpty(self.refseq_transcript)
        transcript.protein_transcript = self._noneIfEmpty(self.protein_transcript)
        transcript.reference_context = self._noneIfEmpty(self.reference_context)

        return transcript

    def _noneIfEmpty(self, value: str):
        '''
        Return None if the string is an empty string.
        ''' 
        if value == '':
            return None
        return value
    
    def __lt__(self, other):
        '''
        The overloaded ``__lt__`` function is used to compare variant transcripts by the number of fields that are non-null. 
        
        Only objects that have the same variant genotype and transcript may be compared (transcript versions may differ). 
        Primary comparison is determined by score (see _get_self_score()). When scores are the same then transcript version 
        is compared.   
        '''
        self_genotype, self_transcript_version = self.get_label().split('.')
        other_genotype, other_transcript_version = other.get_label().split('.')
        
        # Comparison can only be between same variants with same transcript (minus the transcript version)
        assert self_genotype == other_genotype, f"attempt to compare different variants: {self_genotype} and {other_genotype}"
                
        self_score = self._get_self_score()
        other_score = other._get_self_score()
        
        if self_score != other_score:
            return self_score < other_score
        
        # If scores are the same then the comparison is made between transcript versions
        return int(self_transcript_version) < int(other_transcript_version)
        
    def _get_self_score(self):
        '''
        Return a score based on how many the fields are filled in. This method can be improved in the following ways
        - Give p_dot and c_dot a higher score when the values come from HGVS/UTA rather than Annovar
        - Scoring some fields (e.g. g-dot, protein protein script, exon) may not be relevent for non coding variants. 
        '''
        score = 0
        
        # cDNA Fields
        if self.refseq_transcript:
            score += 1

        if self.hgvs_c_dot:
            score += 1
            
        if self.variant_type:
            score += 1
            
        if self.hgvs_base_position:
            score += 1

        if self.exon:
            score += 1
        
        if self.hgnc_gene:
            score += 1
                
        # Protein fields
        if self.protein_transcript:
            score += 1
            
        if self.hgvs_amino_acid_position:
            score += 1

        if self.hgvs_p_dot_one:
            assert self.hgvs_p_dot_three, "Failure to set values for p1 and p3 suggests a bug somewhere"
            assert self.protein_transcript, "Protein genotype is of no use without protein transcript"
            score += 1
                
        if self.hgvs_p_dot_three:
            assert self.hgvs_p_dot_one, "Failure to set values for p3 and p1 suggests a bug somewhere"
            score += 1
            
        if self.variant_effect:            
            score += 1

        return score