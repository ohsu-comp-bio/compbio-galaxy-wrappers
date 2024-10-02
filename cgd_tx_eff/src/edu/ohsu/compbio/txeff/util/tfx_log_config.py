'''
Created on Feb. 14, 2023
Logger configurations that can be passed to the logging.config.dictConfig configurator.
@author: pleyte
'''

class TfxLogConfig(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        # Configuration for the Transcript Effects tool 
        self.log_config = { 
            'version': 1,
            'disable_existing_loggers': False,
            'formatters': {
                'standard': { 
                    'format': '%(asctime)s %(levelname)s: %(name)s::%(module)s:%(lineno)s: %(message)s',
                    'datefmt': '%Y-%m-%d %H:%M:%S'
                },
            },
            'handlers': {
                'default': {
                    'formatter': 'standard',
                    'class': 'logging.FileHandler',
                    'filename': 'cgd_tx_eff.log', 
                    'mode': 'a'
                },
            },
            'loggers': { 
                '': {  # root logger
                    'level': 'WARNING',
                    'handlers': ['default'],
                    'propagate': False
                },
                'edu.ohsu.compbio': { 
                    'level': 'INFO',
                    'handlers': ['default'],
                    'propagate': False
                },
            } 
        }
        
        # Configuration for helper utilities in this source tree 
        self.utility_config = { 
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
                    'stream': 'ext://sys.stdout',
                },
            },
            'loggers': { 
                '': {  # root logger
                    'level': 'WARNING',
                    'handlers': ['default'],
                    'propagate': False
                },
                'edu.ohsu.compbio': { 
                    'level': 'DEBUG',
                    'handlers': ['default'],
                    'propagate': False
                },
            } 
        }
        
        # Configuration for when you want to see logging on stdout 
        self.stdout_config = { 
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
                    'level': 'WARNING',
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