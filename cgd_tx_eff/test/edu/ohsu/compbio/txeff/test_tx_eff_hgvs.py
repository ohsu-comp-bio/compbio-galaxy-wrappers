import unittest
import edu.ohsu.compbio.txeff.tx_eff_hgvs as tx_eff_hgvs 
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript

class TxEffHgvsTest(unittest.TestCase):
    
    def test_merge_annovar_with_hgvs_dont_require_match(self):
        '''
        Merge two lists of transcripts where one matches and the other doesn't; keep all three. 
        '''
        annovar_a = self._create_variant('1', 1, 'A', 'T', 'aaa.1', gene='gene_a', variant_effect='stoploss', variant_type='exonic')
        annovar_b = self._create_variant('2', 2, 'A', 'T', 'bbb.1', gene='gene_b', variant_effect='nonsynonymous SNV', variant_type='intronic')
        
        hgvs_b = self._create_variant('2', 2, 'A', 'T', 'bbb.1', c_dot='c.1A>T', p_dot_one='p.L1P', p_dot_three='p.Leu1Pro', protein_transcript='NP_001')
        hgvs_c = self._create_variant('3', 3, 'A', 'T', 'ccc.1', c_dot='c.3A>T', p_dot_one='p.L3P', p_dot_three='p.Leu3Pro', protein_transcript='NP_003')
    
        merged_transcripts = tx_eff_hgvs._merge_annovar_with_hgvs([annovar_a, annovar_b], [hgvs_b,hgvs_c], False)
        self.assertEqual(len(merged_transcripts), 3, 'Merged list should have one merged transcript and two unmerged transcripts')
        
        
    def test_merge_annovar_with_hgvs_require_match(self):
        '''
        Merge two lists of transcripts and get only those that are common to both lists
        '''
        annovar_a = self._create_variant('1', 1, 'A', 'T', 'aaa.1', gene='gene_a', variant_effect='stoploss', variant_type='exonic')
        annovar_b = self._create_variant('2', 2, 'A', 'T', 'bbb.1', gene='gene_b', variant_effect='nonsynonymous SNV', variant_type='intronic')

        hgvs_b = self._create_variant('2', 2, 'A', 'T', 'bbb.1', c_dot='c.1A>T', p_dot_one='p.L1P', p_dot_three='p.Leu1Pro', protein_transcript='NP_001')
        hgvs_c = self._create_variant('3', 3, 'A', 'T', 'ccc.1', c_dot='c.3A>T', p_dot_one='p.L3P', p_dot_three='p.Leu3Pro', protein_transcript='NP_003')

        merged_transcripts = tx_eff_hgvs._merge_annovar_with_hgvs([annovar_a, annovar_b], [hgvs_b,hgvs_c], True)
        
        self.assertEqual(len(merged_transcripts), 1, 'Only one transcript is common to both lists')
        
    
    def test_get_unmatched_annovar_transcripts(self):
        '''
        Given two dictionary objects, create a list of values based on the the symmetric difference of the keys  
        '''
        a = {'A':'a', 'B':'b', 'C':'c', 'D':'d'}
        b = {'C':'c', 'D':'d', 'E':'e', 'F':'f'}
        c = tx_eff_hgvs._get_unmatched_annovar_transcripts(a,b)
        
        self.assertEqual(c, ['a', 'b', 'e', 'f'], 'Symmetric difference')
        
        
    def _create_variant(self, chromosome, pos, ref, alt, transcript, gene=None, c_dot=None, p_dot_one=None, p_dot_three=None, variant_effect=None, variant_type=None, protein_transcript=None):
        '''
        Helper method for creating variant-transcript objects
        ''' 
        variant_transcript = VariantTranscript(chromosome, pos, ref, alt)
        variant_transcript.refseq_transcript = transcript
        variant_transcript.hgnc_gene = gene
        variant_transcript.hgvs_c_dot = c_dot
        variant_transcript.hgvs_p_dot_one = p_dot_one
        variant_transcript.hgvs_p_dot_three = p_dot_three
        variant_transcript.variant_effect = variant_effect
        variant_transcript.variant_type = variant_type
        variant_transcript.protein_transcript = protein_transcript
        
        return variant_transcript
        