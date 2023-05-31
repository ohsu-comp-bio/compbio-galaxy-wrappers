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
                    'format': '%(levelname)s: %(name)s::%(module)s:%(lineno)s: %(message)s'
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
                    'handlers': ['default'],
                    'level': 'DEBUG',
                    'propagate': False
                },
                'edu.ohsu.compbio.txeff': { 
                    'level': 'DEBUG',
                    'propagate': False,
                    'handlers': ['default'],
                },
            } 
        }
        
        # Configuration for helper utilties in this source tree 
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
                    'handlers': ['default'],
                    'level': 'DEBUG',
                    'propagate': False
                },
                'edu.ohsu.compbio.txeff': { 
                    'level': 'DEBUG',
                    'propagate': False,
                    'handlers': ['default'],
                },
            } 
        }
        
        
        