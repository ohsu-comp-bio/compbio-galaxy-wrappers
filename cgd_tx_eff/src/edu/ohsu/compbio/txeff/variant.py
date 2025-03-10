'''
A variant defined by chromosome, postion, ref and alt

Created on May 2, 2022

@author: pleyte
'''

class Variant(object):
    '''
    A simple variant  
    '''
    def __init__(self, chromosome: str, position: int, ref: str, alt: str, _id = None):
        self.chromosome = chromosome
        self.position = int(position)
        self.reference = ref
        self.alt = alt
        self._id = _id

    def __str__(self):
        '''
        String representation of the AnnovarVariantFunction object
        '''
        return f'[Variant: {self.chromosome}-{self.position}-{self.reference}-{self.alt}]'
    
    def __eq__(self, other):
        
        return self.chromosome == other.chromosome  and \
            self.position == other.position and         \
            self.reference == other.reference and                  \
            self.alt == other.alt
    
    def __hash__(self):
        return hash((self.chromosome, self.position, self.reference, self.alt))