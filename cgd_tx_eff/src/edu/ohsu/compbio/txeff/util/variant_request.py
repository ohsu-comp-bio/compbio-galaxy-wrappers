'''
Created on Mar 14, 2025

@author: pleyte
'''
import argparse
import json

from edu.ohsu.compbio.txeff import tx_eff_control

VERSION = '0.0.1'

class VariantRequest(object):
    '''
    Create a json message used for requesting variants that need to be updated with transcript effects  
    '''
    def get_tfx_version(self):
        '''
        Return the version of the Transcript Effects control module 
        '''
        return tx_eff_control.VERSION
    
    def get_request(self, run_id, dna_barcode_id, tfx_version, batch_size = None) -> dict:
        '''
        Place the parameters in a dict object
        '''
        return {'runId': run_id, 'dnaBarcodeId': dna_barcode_id, 'tfxVersion': tfx_version, 'batchSize': batch_size} 
        
    def write(self, out_stream, request: dict):
        '''
        '''
        json.dump(request, out_stream, indent=2)
    
def _parse_args():
    parser = argparse.ArgumentParser(description='Create the json file used for request a list of variants in CGD that need transcript effects')

    parser.add_argument('--out', 
                dest="output",
                help='json file to place request parameters', 
                type=argparse.FileType('w'), 
                required=True)
    
    parser.add_argument("--run_id", help="Sample run id", type=str, required=True)
    parser.add_argument("--dna_barcode_id", help="Sample run DNA barcode id", type=str, required=True)
    parser.add_argument("--batch_size", help="Maximum number of variants to receive", type=str, default=None)

    return parser.parse_args()

if __name__ == '__main__':
    args = _parse_args()
    
    vr = VariantRequest()
    tfx_version = vr.get_tfx_version()
    request = vr.get_request(args.run_id, args.dna_barcode_id, tfx_version, args.batch_size)
    vr.write(args.output, request)
    print(f"Request parameters written to {args.output.name}")