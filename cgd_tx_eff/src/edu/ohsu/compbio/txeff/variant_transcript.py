from edu.ohsu.compbio.annovar.annovar_parser import AnnovarVariantFunction

class VariantTranscript(AnnovarVariantFunction):
    
    def __init__(self, chromosome, position, ref, alt):
        '''
        Create new VariantTranscript using chromosome, position, ref, and alt.
        '''
        super().__init__(chromosome, position, ref, alt)
            
        self.protein_transcript = None

    def __eq__(self, obj):
        if super.__eq__(self,obj):
            return self.protein_transcript == obj.protein_transcript
        else:           
            return False
    
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
        transcript.hgvs_c_dot = self._noneIfEmpty(self.hgvs_c_dot)
        transcript.hgvs_p_dot_one = self._noneIfEmpty(self.hgvs_p_dot_one)
        transcript.hgvs_p_dot_three = self._noneIfEmpty(self.hgvs_p_dot_three)
        transcript.splicing = self._noneIfEmpty(self.splicing)
        transcript.refseq_transcript = self._noneIfEmpty(self.refseq_transcript)
        transcript.protein_transcript = self._noneIfEmpty(self.protein_transcript)
        return transcript

    def _noneIfEmpty(self, value: str):
        '''
        Return None if the string is an empty string.
        ''' 
        if value == '':
            return None
        return value
