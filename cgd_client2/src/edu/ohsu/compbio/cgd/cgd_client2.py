'''
Created on Mar 19, 2025
@author: pleyte
'''
import argparse
import json
import logging.config
import os
from urllib.parse import urljoin

import requests
from requests_oauthlib import OAuth1


VERSION = '0.0.1'
class CgdClient2(object):
    '''
    classdocs
    '''
    def __init__(self, uri, secret):
        '''
        Constructor
        '''
        self.logger = logging.getLogger(__name__)        
        self._uri = uri
        self._secret = secret
        self._result = None

    def _get_auth(self):
        return OAuth1(client_key = 'cgd-key', client_secret = self._secret, signature_method="HMAC-SHA1")

    def _get_uri(self, path):
        """
        """
        return urljoin(self._uri, path)

    def get(self, path):
        """
        Send a GET request to CGD
        """
        uri = self._get_uri(path)        
        self._result = requests.get(uri, auth=self._get_auth())
        self.logger.info(f"GET {uri} result: {self._result}")

    def post(self, path, file):
        '''
        '''
        uri = self._get_uri(path)       
        data = json.loads(file.read())        
        self._result = requests.post(uri, json = data, auth = self._get_auth()) 
        self.logger.info(f"GET {uri} result: {self._result}")

    def log_simple_response(self):
        '''
        '''
        j = json.loads(self._result.content)
        if j['errors']:
            self.logger.warning(f"SimpleMessage: {j}")
        else:
            self.logger.info(f"SimpleMessage: {j}")

    def write(self, output):
        """
        Write the response to file 
        """
        if not self._result:
            raise ValueError("write called before GET or POST request")

        print(f"jDebug: result content={self._result.content}")

        output.write(self._result.content)
        self.logger.info(f"Wrote response to {output.name}")

def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Communicate with CGD')

    parser.add_argument('--in', 
                        help='Send file to CGD',
                        dest="input_file",
                        type=argparse.FileType('r'),
                        required=False)

    parser.add_argument('--base_uri', 
                        help='Base URI (eg http://localhost:8080/cgd/)',
                        required=True)

    parser.add_argument('--path', 
                        help='API path. The first character should not be a slash',
                        required=True)

    parser.add_argument('--simple',
                    help='Indicate that CGD''s response will be a SimpleResponse', 
                    action='store_true')

    parser.add_argument('--out', 
                        dest="output_file",
                        help='Write CGD''s response to file',
                        type=argparse.FileType('wb'),
                        required=False)

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    args = parser.parse_args()

    if args.path.startswith('/'):
        logger.warning("The path begins with a forward slash which will cause any path in the base_uri to be removed")

    if not args.simple and not args.output_file:
        logger.warning("Output type not specified")

    return args

stdout_log_config = { 
            'version': 1,
            'disable_existing_loggers': False,
            'formatters': {
                'standard': { 
                    'format': '%(levelname)s: %(name)s::%(module)s:%(lineno)s: %(message)s'
                },
            },
            'handlers': {
                'default': {                     
                    'formatter': 'standard',
                    'class': 'logging.StreamHandler',
                    'stream': 'ext://sys.stdout'
                },
            },
            'loggers': { 
                '': {  # root logger
                    'level': 'INFO',
                    'handlers': ['default'],
                    'propagate': False
                },
                'edu.ohsu.compbio': { 
                    'level': 'DEBUG',
                    'handlers': ['default'],
                    'propagate': False,
                },
            }
        }
if __name__ == '__main__':
    logging.config.dictConfig(stdout_log_config)
    logger = logging.getLogger("edu.ohsu.compbio.cgd.cgd_client2")

    args = _parse_args()

    secret = os.environ.get('CGD_OAUTH_SECRET')
    if not secret:
        raise ValueError("OAuth secret must be set in environment variable CGD_OAUTH_SECRET")

    cgd_client = CgdClient2(args.base_uri, secret)

    if args.input_file:
        cgd_client.post(args.path, args.input_file)
        args.input_file.close()
    else:
        cgd_client.get(args.path)

    if args.simple:        
        cgd_client.log_simple_response()

    if args.output_file:
        cgd_client.write(args.output_file)
        args.out.close()