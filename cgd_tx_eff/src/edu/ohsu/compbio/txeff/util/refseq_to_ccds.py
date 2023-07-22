'''
Created on Aug. 24, 2022

@author: pleyte
'''

import csv
import logging.config
import sys

from argparse import ArgumentParser, RawDescriptionHelpFormatter, FileType
from collections import defaultdict
from BCBio import GFF
from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig

__version__ = 0.1

class RefseqToCcds(object):
    '''
    This class handles operations related to mapping from a RefSeq id to a CCDS id.
    The mappings can be found in NCBI's gff file: https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
    
    Running the refseq_to_ccdcs.py module from the command line takes the NCBI gff file and writes out a csv file with only the mappings, which significantly improves performance.      
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.logger = logging.getLogger(__name__)

    def get_mappings_from_gff(self, gff_file_name: str):
        '''
        Read the GFF file and return all refseq-to-ccds mappings 
        '''
        self.logger.info(f"Reading {gff_file_name}")
        refseq_ccds_map = defaultdict(set)
        
        limit_info = dict(gff_type=['CDS'])
        with open(gff_file_name, 'r') as gff_file:
            for rec in GFF.parse(gff_file, limit_info=limit_info):
                for feature in rec.features:
                    # There are two types of records that provide refseq to CCDS mappings:
                    if feature.id.startswith('rna-'):                        
                        # Remove the 'rna-' prefix and the occasional "-n" suffix (eg rna-NM_123.2-1)
                        ref_seq_id_parts = feature.id.split('-')                    
                        ref_seq_id = ref_seq_id_parts[1]
                        
                        if feature.sub_features:
                            for sub_feature in feature.sub_features:        
                                # jDebug: change this; use   _getDbxrefCcds(l) instead                      
                                dbXrefs = sub_feature.qualifiers.get('Dbxref')
                                for dbxref in dbXrefs:
                                    if(dbxref.startswith('CCDS')):
                                        ccds_id = dbxref.split(':')[1]
                                        refseq_ccds_map[ref_seq_id].add(ccds_id)

                    elif feature.id.startswith('cds-'):
                        assert len(feature.qualifiers.get('Parent')) == 1
                        parent = feature.qualifiers.get('Parent')[0]
                        if parent.startswith('rna-'):
                            ref_seq_id = self._get_ref_seq_id(parent)
                            if ref_seq_id:
                                dbxrefs = feature.qualifiers.get('Dbxref')
                                dbxref_ccds_id = self._get_ccds_id(dbxrefs) 
                                if dbxref_ccds_id:
                                    refseq_ccds_map[ref_seq_id].add(dbxref_ccds_id)
                                                
        return refseq_ccds_map
    
    def _get_ref_seq_id(self, s: str):        
        '''
        Extract the refseq id from a string like "rna-NM_001005484.2". Sometimes the refseq id is the GFF row id and there is an 
        extra number like this "rna-NM_178129.5-2". 
        Sometimes the refseq isn't a refseq id, as in the case of cds-YP_003024026.1 and "rna-ND1". If the refseq id doesn't 
        start with "NM_" then None is returned. 
        '''
        ref_seq_id_parts = s.split('-')
        ref_seq_id = ref_seq_id_parts[1]
        
        return ref_seq_id if ref_seq_id.startswith('NM_') else None
        
    def _get_ccds_id(self, dbxrefs: list):
        '''
        Extract the CCDS id from a string like Dbxref="CCDS:CCDS30547.1,GeneID:79501,Genbank:NP_001005484.2,HGNC:HGNC:14825;Name=NP_001005484.2;gbkey=CDS;gene=OR4F5"        
        '''
        l = []
        for dbxref in dbxrefs:
            if(dbxref.startswith('CCDS')):
                # Extract the ccds id from the string that looks like 'CCDS:CCDS72675.1'
                ccds_id = dbxref.split(':')[1]
                l.append(ccds_id)
                
        if len(l) == 0:
            return None
        elif len(l) == 1:
            assert l[0].startswith("CCDS"), f"Doesn't start with CCDS: {l[0]}"
            return l[0]
        else:
            raise Exception(f"More than one CCDS id found when only zero or one expected: {dbxrefs}")
                
    def get_mappings_from_csv(self, csv_file_name: str):
        '''
            Read the csv file created by write_mappings, and return thm in a map object.
        '''
        refseq_ccds_map = defaultdict()
        
        with open(csv_file_name) as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                if refseq_ccds_map.get(row['refseq_id']) is not None:
                    raise Exception(f"RefSeq id {row['refseq_id']} is already mapped to {refseq_ccds_map.get(row['refseq_id'])}. Cannot add additional mapping to {row['ccds_id']}")
                
                refseq_ccds_map[row['refseq_id']] = row['ccds_id']
        
        self.logger.info(f"Read RefSeq-CCDS map of size {len(refseq_ccds_map)} from CSV")
        
        return refseq_ccds_map
        
        
    def write_mappings(self, refSeq_to_ccds: map, csv_file_name: str):
        '''
        Write the RefSeq-to-CCDS map to CSV file. The map seems to be a RefSeq id mapped to a list containing exactly one CCDS id. 
        '''
        self.logger.info(f"Writing to {csv_file_name}")
        
        fields = ['refseq_id', 'ccds_id']
        with open(csv_file_name, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(fields)
            
            for refseq_id, ccds_ids in refSeq_to_ccds.items():
                if(len(ccds_ids) != 1):
                    raise Exception(f"RefSeq ids are expected to map to exactly one CCDS id, but {len(ccds_ids)} found: {ccds_ids}")
                for value in refSeq_to_ccds.get(refseq_id):
                    csv_writer.writerow([refseq_id, value])

        self.logger.info(f"Wrote {len(refSeq_to_ccds)} refseq-to-cccds mappings")

def _main():
    logging.config.dictConfig(TfxLogConfig().stdout_config)
    if logging.root.handlers[0].stream: 
        output = str(logging.root.handlers[0].stream.name)
    else:
        output = logging.root.handlers[0].baseFilename
        
    print(f"Log level={logging.root.getEffectiveLevel()}, output={output}")
    
    parser = ArgumentParser(description='', formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-c", "--csv", help="CSV file to write out RefSeq-to-CCDS mappings", type=FileType('w'), required=True)
    parser.add_argument("-g", "--gff", help="GFF file from which to read RefSeq-to-CCDS mappings", type=FileType('r'), required=True)
    parser.add_argument('-V', '--version', action='version', version=__version__)

    # Process arguments
    args = parser.parse_args()
    
    refSeqToCcds = RefseqToCcds()
    
    # Read mappings from gff
    refSeq_to_ccds_map = refSeqToCcds.get_mappings_from_gff(args.gff.name)
    
    # Write mappings to csv 
    refSeqToCcds.write_mappings(refSeq_to_ccds_map, args.csv.name)
    
    # Read mappings from csv for validation
    csv_mappings = refSeqToCcds.get_mappings_from_csv(args.csv.name)
    assert len(refSeq_to_ccds_map.keys()) == len(csv_mappings.keys()), f"The size of the GFF and CSV maps should be the same but {len(refSeq_to_ccds_map.keys())} != {len(csv_mappings.keys())}" 
    
    return
    
if __name__ == "__main__":    
    sys.exit(_main())
        