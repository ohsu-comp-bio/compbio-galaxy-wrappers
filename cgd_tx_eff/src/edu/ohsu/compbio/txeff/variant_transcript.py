from edu.ohsu.compbio.annovar.annovar_parser import AnnovarVariantFunction

class VariantTranscript(AnnovarVariantFunction):
    
    def __init__(self, chromosome, position, ref, alt):
        super().__init__(chromosome, position, ref, alt)
    
        self.protein_transcript = None
    