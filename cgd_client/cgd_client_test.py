'''
Created on Apr 9, 2025

@author: pleyte
'''
import unittest
from unittest.mock import patch

import cgd_client
import pathlib as pl

JAVA_BIN = '/usr/bin/java'
CGD_CLIENT_JAR = 'cgd_client-1.2.10.jar'
CGD_CLIENT_PROPERTIES = 'cgd_client.properties'
SERVICE_BASE = 'http://localhost:8080/cgd'

# Run id and barcode id must be fore a run that has a patient 
RUN_ID = ''
BARCODE_ID = ''

class CgdClientTest(unittest.TestCase):

    @patch('sys.argv', ['cgd_client.py', 
                        '--endpoint', 'requestVariants', 
                        '--json_out', '/tmp/variants.json', 
                        '--java8_path', JAVA_BIN, 
                        '--cgd_client', CGD_CLIENT_JAR, 
                        '--cgd_config', CGD_CLIENT_PROPERTIES, 
                        '--servicebase', SERVICE_BASE, 
                        '--pipeline_out', 'test-data/variant_request.json'])
    def test_main_endpoint_variantRequest(self):
        cgd_client.main()
        
        self.assertTrue(pl.Path("/tmp/variants.json").is_file(), "Variant request response file")

    @patch('sys.argv', ['cgd_client.py', 
                        '--endpoint', 'reportedvariants', 
                        '--java8_path', JAVA_BIN, 
                        '--cgd_client', CGD_CLIENT_JAR, 
                        '--cgd_config', CGD_CLIENT_PROPERTIES, 
                        '--servicebase', SERVICE_BASE,
                        '--runid', RUN_ID,
                        '--barcodeid', BARCODE_ID, 
                        '--report_vcf', '/tmp/reported_variants.vcf',
                        '--report_bed', '/tmp/reported_variants.bed'])
    def test_main_endpoint_reportedvariants(self):
        cgd_client.main()
        self.assertTrue(pl.Path("/tmp/reported_variants.vcf").is_file(), "reported variants vcf")
        self.assertTrue(pl.Path("/tmp/reported_variants.bed").is_file(), "reported variants bed")
            
if __name__ == "__main__":
    unittest.main()