'''
Created on Aug. 1, 2024

Provide mapping from refseq accession to chromosome number. 
@author: pleyte
'''

refseq_to_ncbi = {
    'NC_012920.1': 'MT', 
    'NC_000001.10': '1',
    'NC_000002.11': '2',     
    'NC_000003.11': '3',
    'NC_000004.11': '4', 
    'NC_000005.9': '5',
    'NC_000006.11': '6', 
    'NC_000007.13': '7',
    'NC_000008.10': '8', 
    'NC_000009.11': '9',
    'NC_000010.10': '10', 
    'NC_000011.9': '11',
    'NC_000012.11': '12', 
    'NC_000013.10': '13',
    'NC_000014.8': '14', 
    'NC_000015.9': '15',
    'NC_000016.9': '16', 
    'NC_000017.10': '17',
    'NC_000018.9': '18', 
    'NC_000019.9': '19',
    'NC_000020.10': '20', 
    'NC_000021.8': '21',
    'NC_000022.10': '22', 
    'NC_000023.10': 'X',
    'NC_000024.9': 'Y'
}

ncbi_to_refseq = {
    'MT': 'NC_012920.1', 
    '1': 'NC_000001.10',
    '2': 'NC_000002.11',     
    '3': 'NC_000003.11',
    '4': 'NC_000004.11', 
    '5': 'NC_000005.9',
    '6': 'NC_000006.11', 
    '7': 'NC_000007.13',
    '8': 'NC_000008.10', 
    '9': 'NC_000009.11',
    '10': 'NC_000010.10', 
    '11': 'NC_000011.9',
    '12': 'NC_000012.11', 
    '13': 'NC_000013.10',
    '14': 'NC_000014.8', 
    '15': 'NC_000015.9',
    '16': 'NC_000016.9', 
    '17': 'NC_000017.10',
    '18': 'NC_000018.9', 
    '19': 'NC_000019.9',
    '20': 'NC_000020.10', 
    '21': 'NC_000021.8',
    '22': 'NC_000022.10', 
    'X': 'NC_000023.10',
    'Y': 'NC_000024.9' 
}

def get_refseq(chromosome_in_ncbi_format: str) -> str:
    '''
    Return refseq chromosome that is mapped to the requested NCBI formatted chromosome. If the requested chromosome is not 
    mapped to anything then an error is raised.  
    ''' 
    if chromosome_in_ncbi_format.startswith('chr'):
        chromosome_in_ncbi_format = chromosome_in_ncbi_format.replace('chr', '', 1)
    
    refseq_chromosome = ncbi_to_refseq.get(chromosome_in_ncbi_format)
    
    if not refseq_chromosome:
        raise ValueError('Unrecognized RefSeq chromosome' + chromosome_in_ncbi_format)
    
    return refseq_chromosome 

def get_ncbi(chromosome_in_refseq_format: str, include_chr_prefix = False) -> str:
    '''
    Return NCBI chromosome that is mapped to the requested RefSeq formatted chromosome. If the requested chromosome is not 
    mapped to anything then an error is raised. 
    '''
    ncbi_chromosome = refseq_to_ncbi.get(chromosome_in_refseq_format)
    
    if not ncbi_chromosome:
        raise ValueError('Unrecognized NCBI chromosome ' + chromosome_in_refseq_format)
    
    return 'chr'+ncbi_chromosome if include_chr_prefix else ncbi_chromosome

def get_ordinal(chromosome: str): 
    """
    Return the chromosome as an integer 
    """
    if chromosome == 'X':
        return 23
    elif chromosome == 'Y':
        return 24
    elif chromosome == 'MT':
        return 25
    else:
        return int(chromosome)
    