'''
Created on Apr 9, 2025

@author: pleyte
'''
import unittest
from unittest.mock import patch

import cgd_client
import pathlib as pl


class CgdClientTest(unittest.TestCase):

    @patch('sys.argv', ['cgd_client.py', 
                        '--endpoint', 'requestVariants', 
                        '--json_out', '/tmp/variants.json', 
                        '--java8_path', '/Library/Java/JavaVirtualMachines/openjdk-8.jdk/Contents/Home/bin/java', 
                        '--cgd_client', '/Users/pleyte/.m2/repository/edu/ohsu/kdl/cgd_client/1.2.10-SNAPSHOT/cgd_client-1.2.10-SNAPSHOT.jar', 
                        '--cgd_config', '/Users/pleyte/Documents/office/ohsu/cgd_client/cgd_client.properties', 
                        '--servicebase', 'https://kdlcgdwebdev2.ohsu.edu/tfx_cgd/', 
                        '--pipeline_out', 'test-data/variant_request.json'])
    def test_main_endpoint_variantRequest(self):
        cgd_client.main()
        
        self.assertTrue(pl.Path("/tmp/variants.json").is_file(), "Variant request response file")

    @patch('sys.argv', ['cgd_client.py', 
                        '--endpoint', 'reportedvariants', 
                        '--java8_path', '/Library/Java/JavaVirtualMachines/openjdk-8.jdk/Contents/Home/bin/java', 
                        '--cgd_client', '/Users/pleyte/.m2/repository/edu/ohsu/kdl/cgd_client/1.2.10-SNAPSHOT/cgd_client-1.2.10-SNAPSHOT.jar', 
                        '--cgd_config', '/Users/pleyte/Documents/office/ohsu/cgd_client/cgd_client.properties', 
                        '--servicebase', 'https://kdlcgdwebdev2.ohsu.edu/tfx_cgd/',
                        '--runid', '231121_NS500390_0123_AHG7WKBGXV',
                        '--barcodeid', 'E02', 
                        '--report_vcf', '/tmp/reported_variants.vcf',
                        '--report_bed', '/tmp/reported_variants.bed'])
    def test_main_endpoint_reportedvariants(self):
        cgd_client.main()
        self.assertTrue(pl.Path("/tmp/reported_variants.vcf").is_file(), "reported variants vcf")
        self.assertTrue(pl.Path("/tmp/reported_variants.bed").is_file(), "reported variants bed")
    
    def test_main_endpoint_snpProfile(self):
        cgd_client.main()
            
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()