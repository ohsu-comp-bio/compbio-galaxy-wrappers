'''
https://stackoverflow.com/questions/1896918/running-unittest-with-typical-test-directory-structure

rename "test" directory to "tests":
* It is recommended to use "tests" instead of "test", because "test" is a Python build in module. 
* https://docs.python.org/2/library/test.html


'''
import unittest
from edu.ohsu.compbio.annovar.annovar_parser import AnnovarParser

class AnnovarParserTest(unittest.TestCase):
        
    def test_unpack_vf_transcript_tuple(self):
        '''
        Test each of the three tuple types that a *.variant_function file uses for transcripts 
        '''
        annovar_parser = AnnovarParser();
    
        # Type I
        refseq_transcript, exon, hgvs_c_dot, hgvs_basep = annovar_parser._unpack_vf_transcript_tuple("NM_001039724.3")
        self.assertEqual(refseq_transcript, "NM_001039724.3")
        self.assertIsNone(exon)
        self.assertIsNone(hgvs_c_dot)        
        self.assertIsNone(hgvs_basep)
        
        # Type II
        refseq_transcript, exon, hgvs_c_dot, hgvs_basep = annovar_parser._unpack_vf_transcript_tuple("NM_001290354.2(NM_001290354.2:c.-5432_-5431insGCGCTGCGG)")
        self.assertEqual(refseq_transcript, "NM_001290354.2")
        self.assertIsNone(exon)
        self.assertEqual(hgvs_c_dot, "c.-5432_-5431insGCGCTGCGG")        
        self.assertEqual(hgvs_basep, '-5432')
        
        # Type III
        refseq_transcript, exon, hgvs_c_dot, hgvs_basep = annovar_parser._unpack_vf_transcript_tuple("NM_001206844.2(NM_001206844.2:exon5:c.1135-4T>C)")
        self.assertEqual(refseq_transcript, 'NM_001206844.2')
        self.assertEqual(exon, 5)
        self.assertEqual(hgvs_c_dot, 'c.1135-4T>C')        
        self.assertEqual(hgvs_basep, '1135-4')
        
        # Type V
        refseq_transcript, exon, hgvs_c_dot, hgvs_basep = annovar_parser._unpack_vf_transcript_tuple("NM_001242366.3(NM_001242366.3:exon3:c.1142+1T>-,NM_001242366.3:exon4:c.1143-1T>-)")
        self.assertEqual(refseq_transcript, 'NM_001242366.3')
        self.assertEqual(exon, 3)
        self.assertEqual(hgvs_c_dot, 'c.1142+1T>-')        
        # For some reason hgvs.parser.Parser.parse can't figure out the base pair in that c.
        self.assertEqual(hgvs_basep, None)
        
        # Invalid format: empty 
        self.assertRaises(ValueError, annovar_parser._unpack_vf_transcript_tuple, "")
        
        # Invalid format: delimited but still not a proper Type V  
        refseq_transcript, exon, hgvs_c_dot, hgvs_basep = annovar_parser._unpack_vf_transcript_tuple("a:b:c:d:e:f")
        self.assertIs(refseq_transcript, None)
        self.assertIs(exon, None)
        self.assertIs(hgvs_c_dot, None)
        self.assertIs(hgvs_basep, None)
         
        refseq_transcript, exon, hgvs_c_dot, hgvs_basep = annovar_parser._unpack_vf_transcript_tuple("a:b:c:d:e")
        self.assertIs(refseq_transcript, None)
        self.assertIs(exon, None)
        self.assertIs(hgvs_c_dot, None)
        self.assertIs(hgvs_basep, None)
                
    def test__parse_variant_function_row(self):
        '''
        Parse the list transcripts information that is found in an annovar variant_function file
        '''
        data = ("NM_001190458.2,"
                "NM_001077690.1(NM_001077690.1:exon1:c.60+1C>-,NM_001077690.1:exon2:c.61-1C>-),"
                "NM_001289861.1(NM_001289861.1:exon20:c.3002-2A>G)")
        annovar_parser = AnnovarParser();
        transcript_tuples = annovar_parser._split_variant_function_transcript_tuples(data)        
        self.assertEqual(len(transcript_tuples), 3)        
    
    def test__parse_transcript_tuple_type5_two_exons(self):
        '''
        Parse the transcript definition that is defined at two different exons, and the left one is chosen.  
        '''        
        data = "NM_001077690.1(NM_001077690.1:exon1:c.60+1C>-,NM_001077690.1:exon2:c.61-1C>-)"
        
        annovar_parser = AnnovarParser();
        refseq_transcript, exon, hgvs_c_dot = annovar_parser._parse_transcript_tuple_type5(data)
        
        self.assertEqual(refseq_transcript, 'NM_001077690.1', "Incorrect refseq transcript")
        self.assertEqual(exon, 'exon1', "Incorrect exon")
        self.assertEqual(hgvs_c_dot, 'c.60+1C>-', "Incorrect c.")
    
    def test__parse_transcript_tuple_type5_invalid(self):
        '''
        The "Type 5" tuple type defines a transcript that can be at two or more exons. This test function makes sure
        an exception would be thrown if a Type-V format specified only one exon.         
        '''
        # This definition is invalid because it only defines a single exon 
        data = "NM_001077690.1(NM_001077690.1:exon1:c.60+1C>-)"
        annovar_parser = AnnovarParser();
        self.assertRaises(AssertionError, annovar_parser._parse_transcript_tuple_type5, data)
    
    def test__parse_transcript_tuple_type5_six_exons(self):    
        '''
        The "Type 5" tuple type defines a transcript that can be at two or more exons. This test function makes sure
        that the correct exon is returned even when there are six exons.     
        ''' 
        data = "NM_1.2(NM_1.2:exon1:c.1+1C>-,NM_1.2:exon2:c.1-1C>-,NM_1.2:exon3:c.1-1C>-,NM_1.2:exon4:c.1-1C>-,NM_1.2:exon5:c.1-1C>-,NM_1.2:exon6:c.1-1C>-)"        
        annovar_parser = AnnovarParser();
        refseq_transcript, exon, hgvs_c_dot = annovar_parser._parse_transcript_tuple_type5(data)
        self.assertEqual(exon, 'exon1', "Incorrect exon")

    def test__parse_transcript_tuple_type3(self):
        data = "NM_001289861.1(NM_001289861.1:exon20:c.3002-2A>G)"
        
        annovar_parser = AnnovarParser();
        refseq_transcript, exon, hgvs_c_dot = annovar_parser._parse_transcript_tuple_type3(data)
        
        self.assertEqual(refseq_transcript, 'NM_001289861.1', "Incorrect refseq transcript")
        self.assertEqual(exon, 'exon20', "Incorrect exon")
        self.assertEqual(hgvs_c_dot, 'c.3002-2A>G', "Incorrect c.")
    
    def test__parse_transcript_tuple_type2(self):
        data = "NM_000791.4(NM_000791.4:c.-417_-416insGCGCTGCGG)"
        annovar_parser = AnnovarParser();
        refseq_transcript, hgvs_c_dot = annovar_parser._parse_transcript_tuple_type2(data)
        
        self.assertEqual(refseq_transcript, 'NM_000791.4', "Incorrect refseq transcript")        
        self.assertEqual(hgvs_c_dot, 'c.-417_-416insGCGCTGCGG', "Incorrect c.")

    def test__parse_transcript_tuple_type1(self):
        annovar_parser = AnnovarParser()

        data = "NM_000791.4"        
        refseq_transcript = annovar_parser._parse_transcript_tuple_type1(data)        
        self.assertEqual(refseq_transcript, 'NM_000791.4', "Incorrect refseq transcript")        

        data = "NM_001126117.1(dist=778)"
        refseq_transcript = annovar_parser._parse_transcript_tuple_type1(data)        
        self.assertEqual(refseq_transcript, None, "Intergenic transcript should be ignored")        
