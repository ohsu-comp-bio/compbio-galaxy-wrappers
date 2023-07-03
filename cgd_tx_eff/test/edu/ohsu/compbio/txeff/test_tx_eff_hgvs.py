import unittest
import edu.ohsu.compbio.txeff.tx_eff_hgvs as tx_eff_hgvs 
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript
from edu.ohsu.compbio.txeff.util.tx_eff_pysam import PysamTxEff

FILE_HG37_REFERENCE_FASTA = '/opt/bioinformatics/Broad/Homo_sapiens_assembly19.fasta'

class TxEffHgvsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._pysam_file = PysamTxEff(FILE_HG37_REFERENCE_FASTA)
        
    def test__merge_annovar_with_hgvs(self):
        '''
        Merge two lists of transcripts and get only those that are common to both lists
        '''
        annovar_a = self._create_variant('1', 1, 'A', 'T', 'aaa.1', gene='gene_a', variant_effect='stoploss', variant_type='exonic')
        annovar_b = self._create_variant('2', 2, 'A', 'T', 'bbb.1', gene='gene_b', variant_effect='nonsynonymous SNV', variant_type='intronic')

        hgvs_b = self._create_variant('2', 2, 'A', 'T', 'bbb.1', c_dot='c.1A>T', p_dot_one='p.L1P', p_dot_three='p.Leu1Pro', protein_transcript='NP_001')
        hgvs_c = self._create_variant('3', 3, 'A', 'T', 'ccc.1', c_dot='c.3A>T', p_dot_one='p.L3P', p_dot_three='p.Leu3Pro', protein_transcript='NP_003')

        merged_transcripts, unmerged_transcripts = tx_eff_hgvs._merge_annovar_with_hgvs([annovar_a, annovar_b], [hgvs_b,hgvs_c])
        
        self.assertEqual(len(merged_transcripts), 1, 'Only one transcript is common to both lists')
        self.assertEqual(len(unmerged_transcripts), 2, 'Two transcripts are common to both lists')
        
    
    def test__get_unmatched_annovar_transcripts(self):
        '''
        Given two dictionary objects, create a list of values based on the the symmetric difference of the keys  
        '''        
        t1 = self._create_variant('1', 1, 'A', 'T', 'NM_1.1')
        t2 = self._create_variant('2', 2, 'A', 'T', 'NM_2.2')
        t3 = self._create_variant('3', 3, 'A', 'T', 'NM_3.3')
        t4 = self._create_variant('4', 4, 'A', 'T', 'NM_4.4')
        t5 = self._create_variant('5', 5, 'A', 'T', 'NM_5.5')
        t6 = self._create_variant('6', 6, 'A', 'T', 'NM_6.6')
        
        # Dict a contains 1,2,3,4    
        a = {t1.get_label():t1, t2.get_label():t2, t3.get_label():t3, t4.get_label():t4}
        
        # Dict b contains 3,4,5,6
        b = {t3.get_label():t3, t4.get_label():t4, t5.get_label():t5, t6.get_label():t6}
        
        # Unmatched will be 1,2,5,6
        c = tx_eff_hgvs._get_unmatched_annovar_transcripts(a,b)
        
        self.assertEqual(len(c), 4, 'Symmetric difference') 
        
        # Use genotype labels rather than try to compare objects
        labels = [x.get_label() for x in c]
        self.assertTrue(t1.get_label() in labels, 't1 is unmatched')
        self.assertTrue(t2.get_label() in labels, 't2 is unmatched')

        self.assertTrue(t3.get_label() not in labels, 't3 not in list')
        self.assertTrue(t4.get_label() not in labels, 't4 not in list')
        
        self.assertTrue(t5.get_label() in labels, 't5 is unmatched')
        self.assertTrue(t6.get_label() in labels, 't6 is unmatched')
        
    def test__get_the_best_transcripts(self):
        '''
        Given three versions of a transcript return then one with the most values and the latest version
        '''
        merged_transcripts = [self._create_variant('1', 1, 'A', 'T', 'NM_000.1', gene='AAA', p_dot_one='p.(L1P)', variant_effect='X'),
                              self._create_variant('1', 1, 'A', 'T', 'NM_000.3', gene='AAA'),
                              self._create_variant('1', 1, 'A', 'T', 'NM_000.2', gene='AAA', p_dot_one='p.(L1P)', variant_effect='Y')]

        best = tx_eff_hgvs._get_the_best_transcripts(merged_transcripts)
        
        self.assertEqual(len(best), 1, 'Only one transcript expected')
        self.assertEqual(best[0].get_label(), '1-1-A-T-NM_000.2', 'Most complete, lastest version')
    
    def test__correct_indel_coords(self):
        '''
        Given a genotype the _correct_indel_coords function returns the nomenclature describing the variant change (the g., minus the g. prefix) 
        ''' 
        # Substitution 
        genotype = '1-123-G-A'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = tx_eff_hgvs._correct_indel_coords(chromosome, int(position), ref, alt, self._pysam_file)
        self.assertEqual(pos_part, '123G>A', "g. is incorrect")

        # Insertion 
        genotype = '1-123-G-AC'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = tx_eff_hgvs._correct_indel_coords(chromosome, int(position), ref, alt, self._pysam_file)
        self.assertEqual(pos_part, '123_124insC', "g. is incorrect")

        # Deletion
        genotype = '1-123-AC-G'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = tx_eff_hgvs._correct_indel_coords(chromosome, int(position), ref, alt, self._pysam_file)
        self.assertEqual(pos_part, '124del', "g. is incorrect")
        
        # Indel
        genotype = '5-112175461-CAGTTCACTTGA-CGTC'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = tx_eff_hgvs._correct_indel_coords(chromosome, int(position), ref, alt, self._pysam_file)
        self.assertEqual(pos_part, '112175462_112175472delinsGTC', "g. is incorrect")
        
        # Repeat
        genotype = '1-153924029-AG-A'
        chromosome, position, ref, alt = genotype.split('-')
        pos_part = tx_eff_hgvs._correct_indel_coords(chromosome, int(position), ref, alt, self._pysam_file)
        self.assertEqual(pos_part, '153924033del', "g. is incorrect")
        
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
