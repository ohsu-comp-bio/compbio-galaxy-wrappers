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
        
    
    