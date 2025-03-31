'''
Created on Mar 27, 2025

@author: pleyte
'''
import argparse
import json
import logging.config
import os
from urllib.parse import urljoin

import requests
from requests_oauthlib import OAuth1


VERSION = '0.0.0.1'
class CgdApi(object):
    '''
    Interface with CGD's API
    '''
    def __init__(self, service_base, oauth_secret):
        '''
        Constructor
        '''
        self.logger = logging.getLogger(__name__)
        self._service_base = service_base
        self._oauth_secret = oauth_secret
        self._parameters = {}

    def _get_auth(self):
        return OAuth1(client_key = 'cgd-key', client_secret = self._oauth_secret, signature_method="HMAC-SHA1")        

    def _get_uri(self, path_template):
        """
        Build the URI        
        """        
        # Substitute each {name} in the template with the corresponding value from the _parameters dict
        path = path_template.format(**self._parameters)
        
        # the urljoin function needs the service base to end with a slash (eg http://localhost/cgd/) and the path must not begin with a slash.        
        return urljoin(self._service_base, path)
    
    def get(self, path_template):
        '''
        Send a GET request to CGD
        '''
        uri = self._get_uri(path_template)        
        response = requests.get(uri, auth=self._get_auth())
        self.logger.info(f"GET {uri} result: {response}")
        
        if response.status_code != 200:
            raise ValueError(f"Request failed: {response}")
        
        return response

    def post(self, path_template, in_file, in_file_json):
        """
        Send a file using POST
        """
        uri = self._get_uri(path_template)
        
        data = in_file.read() if in_file else None 
        json_data = json.loads(in_file_json.read()) if in_file_json else None
        
        response = requests.post(uri, data=data, json=json_data, auth=self._get_auth())

        self.logger.info(f"GET {uri} result: {response}")

        if response.status_code != 200:
            raise ValueError(f"Request failed: {response.content}")

        return response

    def log_simple_response(self, response):
        '''
        Log a simple response from CGD. A simple response is json with a 'message' attribute and an 'errors' list. 
        '''
        self.logger.info(json.loads(response.content))

    def set_parameter_variables(self, name_value_pairs: list):
        """
        Create a dict object containing the named parameters for substitution in the uri
        The `name_value_pairs` parameter is constructed by argparser and is a list of name=value pairs; one for each --variable parameter        
        """
        if(name_value_pairs):
            for pair in name_value_pairs:
                if '=' not in pair[0]:
                    raise ValueError(f"Expected name=value but found: {pair[0]}")
                name, value = pair[0].split('=')
                self._parameters[name] = value

    def write(self, output, response):
        """
        """
        output.write(response.content)
        self.logger.info(f"Wrote response to {output.name}")
        
def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Communicate with CGD')
    
    parser.add_argument('--service_base', 
                        help='Base URI (eg http://localhost:8080/cgd/)',
                        required=True)

    parser.add_argument('--path', 
                        help='API path. The first character should not be a slash',
                        required=True)
    
    parser.add_argument('--in_file', 
                        help='Send file to CGD',
                        type=argparse.FileType('r'),
                        required=False)
    
    parser.add_argument('--in_file_json', 
                        help='Send json file to CGD',
                        type=argparse.FileType('r'),
                        required=False)
    
    parser.add_argument('--out_file', 
                        help='Write CGD''s response to file',
                        type=argparse.FileType('wb'),
                        required=False)
    
    parser.add_argument('--variable', 
                        help='Name=value pair',
                        required=False,
                        action='append',
                        nargs=1)

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    args = parser.parse_args()

    if args.in_file and args.in_file_json:
        raise ValueError("The --in_file and --in_file_json parameters are mutually exclusive")
    
    if not args.service_base.endswith("/"): 
        logger.warning(f"The service base should have a trailing slash: {args.service_base}")
        
    if args.path.startswith("/"):
        logger.warning(f"The path should not have a leading slash: {args.path}")
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
    logger = logging.getLogger("edu.ohsu.compbio.cgd.cgd_api")

    args = _parse_args()

    secret = os.environ.get('CGD_OAUTH_SECRET')
    if not secret:
        raise ValueError("OAuth secret must be set in environment variable CGD_OAUTH_SECRET")

    cgd_api = CgdApi(args.service_base, secret)
    cgd_api.set_parameter_variables(args.variable)
    
    if args.in_file or args.in_file_json:
        response = cgd_api.post(args.path, args.in_file, args.in_file_json)
    else:
        response = cgd_api.get(args.path)
    
    if args.out_file:
        cgd_api.write(args.out_file, response)
    else:
        cgd_api.log_simple_response(response)
