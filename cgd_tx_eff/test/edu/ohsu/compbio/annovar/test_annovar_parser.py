'''
https://stackoverflow.com/questions/1896918/running-unittest-with-typical-test-directory-structure

rename "test" directory to "tests":
* It is recommended to use "tests" instead of "test", because "test" is a Python build in module. 
* https://docs.python.org/2/library/test.html


'''
import unittest
from edu.ohsu.compbio.annovar.annovar_parser import AnnovarParser
from edu.ohsu.compbio.annovar.annovar_parser import AnnovarFileType

class AnnovarParserTest(unittest.TestCase):

    def test_get_file_type(self):
        annovar_parser = AnnovarParser();
        file_type = annovar_parser._get_file_type(".variant_function")
        self.assertEqual(file_type, AnnovarFileType.VariantFunction)
        
    def test_unpack_vf_transcript_tuple(self):
        '''
        Test each of the three tuple types that a *.variant_function file uses for transcripts 
        '''
        annovar_parser = AnnovarParser();
    
        # Type I
        hgvs_basep, exon, hgvs_c_dot, refseq_transcript = annovar_parser._unpack_vf_transcript_tuple("NM_001039724.3")
        self.assertIsNone(hgvs_basep)
        self.assertIsNone(exon)
        self.assertIsNone(hgvs_c_dot)
        self.assertEqual(refseq_transcript, "NM_001039724.3")
        
        # Type II
        hgvs_basep, exon, hgvs_c_dot, refseq_transcript = annovar_parser._unpack_vf_transcript_tuple("NM_001290354.2(NM_001290354.2:c.-5432_-5431insGCGCTGCGG)")
        self.assertEqual(hgvs_basep.base, -5432)
        self.assertIsNone(exon)
        self.assertEqual(hgvs_c_dot, "c.-5432_-5431insGCGCTGCGG")
        self.assertEqual(refseq_transcript, "NM_001290354.2")
        
        # Type III
        hgvs_basep, exon, hgvs_c_dot, refseq_transcript = annovar_parser._unpack_vf_transcript_tuple("NM_001206844.2(NM_001206844.2:exon5:c.1135-4T>C)")
        self.assertEqual(hgvs_basep.base, 1135)
        self.assertEqual(exon, '5')
        self.assertEqual(hgvs_c_dot, 'c.1135-4T>C')
        self.assertEqual(refseq_transcript, 'NM_001206844.2')
        
        # Unrecognized
        self.assertRaises(Exception, annovar_parser._unpack_vf_transcript_tuple, "a:b:c:d")
        
