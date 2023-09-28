'''
Created on Aug. 24, 2022

@author: pleyte
'''

import csv
import logging.config
import sys

from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig
from argparse import ArgumentParser, RawDescriptionHelpFormatter, FileType
from collections import defaultdict
from BCBio import GFF

import Bio
from sys import getsizeof

__version__ = '0.6.6'
GFF_RELEASE_105 = '105'
GFF_RELEASE_105_20220307 = '105.20220307'

class RefseqToCcds(object):
    '''
    This class handles operations related to mapping from a RefSeq id to a CCDS id.
    The mappings can be found in NCBI's gff3 file: https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
    
    Running the refseq_to_ccdcs.py module from the command line takes the NCBI gff file and writes out a csv file with only the mappings, which significantly improves performance.      
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.logger = logging.getLogger(__name__)

    def get_mappings_from_gff(self, source_release: str, gff_file_name: str):
        '''
        Read the GFF file and return all refseq-to-ccds mappings 
        '''
        self.logger.info(f"Reading {gff_file_name}")
        refseq_ccds_map = defaultdict(set)
        
        # There are two versions of the GFF and and they have to be parsed slightly differently 
        if source_release == GFF_RELEASE_105_20220307:
            limit_info = dict(gff_type=['CDS'])
        elif source_release == GFF_RELEASE_105:
            limit_info = dict(gff_type=['exon', 'CDS'])  
        else: 
            raise Exception(f"Unknown source_release version: '{source_release}")
        
        with open(gff_file_name, 'r') as gff_file:
            for rec in GFF.parse(gff_file, limit_info=limit_info):                
                for feature in rec.features:
                    ref_seq_id, ccds_id = self._get_ids_from_record(source_release, feature)
                    if ref_seq_id and ccds_id:
                        if ref_seq_id in refseq_ccds_map and refseq_ccds_map[ref_seq_id] != ccds_id:
                            raise Exception(f"Duplicate mapping encountered: {ref_seq_id} to {refseq_ccds_map[ref_seq_id]} and {ccds_id}")
                        
                        refseq_ccds_map[ref_seq_id] = ccds_id
                    
        return refseq_ccds_map
    
    def _get_ids_from_record(self, source_release:str, feature: Bio.SeqFeature.SeqFeature):
        '''
        There are mutliple GFF file versions and various record types that can be parse, and this 
        function figures out to handle each line. If the feature contains a RefSeq-to-CCDS mapping then the
        Refseq and CCDS ids are returned. If the line is not of any use then a Null-Null tuple is returned.   
        ''' 
        # The different types of features in the GFF can be identified by the prefix of their ids

        if source_release == GFF_RELEASE_105 and (feature.id.startswith('rna')): 
            return self._get_ids_from_rel105_rna_record(feature)

        elif source_release == GFF_RELEASE_105_20220307 and feature.id.startswith('rna-'):
            return self._get_ids_from_rna_record(feature)

        elif source_release == GFF_RELEASE_105_20220307 and feature.id.startswith('cds-'):
            return self._get_ids_from_cds_record(feature)
    
        else:
            return None, None

    def _get_ids_from_rel105_rna_record(self, feature: Bio.SeqFeature.SeqFeature):
        '''      
        Parse refseq and ccds ids from the older (release 105 GFF) feature.  
        '''
        refseq_ids = set()
        ccds_ids = set()
        
        # Collect all the refseq and ccds ids from each of the sub-features. 
        for sub_feature in feature.sub_features:
            if sub_feature.qualifiers['gbkey'][0] == 'mRNA':
                # There should only be one, but it is stored in a list
                refseq_ids.update(sub_feature.qualifiers['transcript_id'])
            elif sub_feature.qualifiers['gbkey'][0] == 'CDS':
                ccds_ids.add(self._get_ccds_id(sub_feature.qualifiers['Dbxref']))
                    
        # It is possible for a feature to have a refseq id and not a ccds id.
        if len(refseq_ids) == 1 and len(ccds_ids) == 0:
            self.logger.debug(f"Refseq {refseq_ids.pop()} does not have a corresponding CCDS id")
            return None, None
        elif len(ccds_ids) == 1 and len(refseq_ids) == 0:
            raise Exception(f"CCDS {ccds_ids.pop()} does not have a corresponding RefSeq id")
        if len(refseq_ids) > 0 and len(ccds_ids) > 0:
            assert len(refseq_ids) == 1, f"More than one RefSeq id found for feature {feature.id}"
            assert len(ccds_ids) == 1, f"More than one CCDS id found for feature {feature.id}"
            return refseq_ids.pop(), ccds_ids.pop()
        else:
            # No RefSeq or CCDS ids found
            return None, None

    def _get_ids_from_cds_record(self, feature: Bio.SeqFeature.SeqFeature):
        '''
        Return the refseq and corresponding CCDS id when the gff feature has an id starting with 'cds-'. 
        For this type of feature teh CCDS id is in feature.qualifiers['Dbxref'].
        '''
        ref_seq_id = None
        dbxref_ccds_id = None

        parent = feature.qualifiers.get('Parent')[0]
        if parent.startswith('rna-'):
            ref_seq_id = self._get_ref_seq_id(parent)
            if ref_seq_id:
                dbxrefs = feature.qualifiers.get('Dbxref')
                dbxref_ccds_id = self._get_ccds_id(dbxrefs)

        return ref_seq_id, dbxref_ccds_id
    
    def _get_ids_from_rna_record(self, feature: Bio.SeqFeature.SeqFeature):
        '''
        Return the refseq and corresponding CCDS id when the gff feature has an id starting with 'rna' (eg rna-NM_001005484.2).
        For this type of feature the CCDS id is in feature.sub_features[].qualifiers['Dbxref']
        '''
        ref_seq_id = None
        dbxref_ccds_ids = set()
        
        if feature.sub_features:
            ref_seq_id = self._get_ref_seq_id(feature.id)
            for sub_feature in feature.sub_features:
                dbxrefs = sub_feature.qualifiers.get('Dbxref')
                dbxref_ccds_ids.add(self._get_ccds_id(dbxrefs))

        assert len(dbxref_ccds_ids) == 1, f"Refseq {ref_seq_id} has more than one CCDS: {dbxref_ccds_ids}"
        ccds_id = dbxref_ccds_ids.pop()
        return ref_seq_id, ccds_id

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
                
    def get_mappings_from_csv(self, csv_file_name: str) -> dict:
        '''
            Read the csv file created by write_mappings, and return them in a map object.
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

    def write_mappings(self, refSeq_to_ccds: dict, csv_file_name: str):
        '''
        Write the RefSeq-to-CCDS map to CSV file. The map seems to be a RefSeq id mapped to a list containing exactly one CCDS id. 
        '''
        self.logger.info(f"Writing to {csv_file_name}")
        
        fields = ['refseq_id', 'ccds_id']
        
        # Using python csv results in a file with dos line endings and then galaxy doesn't recognise the file as a CSV.        
        with open(csv_file_name, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(fields)
            
            for refseq_id, ccds_id in refSeq_to_ccds.items():
                csv_writer.writerow([refseq_id, ccds_id])
    
    def _merge_maps(self, older_mapping:dict, newer_mapping:dict):
        '''
        Combine the two dicts by adding all the entries in newer_mapping to older_mapping, overwriting any that are already there. 
        '''
        older_size = len(older_mapping)
        newer_size = len(newer_mapping)
        
        older_mapping.update(newer_mapping)
        newest_size = len(older_mapping)
        self.logger.info(f"Combined older map (size {older_size}) with newer map (size {newer_size}) creating new map with (size {newest_size})")

        return older_mapping

def _get_app_operation(args):
    if args.gff and args.source_release and args.csv and not args.merge_old_csv and not args.merge_new_csv:
        logging.debug("Parameters indicate user wants to load GFF and create a CSV file")
        return 'CREATE_MAPPINGS'
    elif args.merge_old_csv and args.merge_new_csv and args.csv and not args.gff and not args.source_release:
        logging.debug("Parameters indicate user wants to merge two CSV files")
        return 'COMBINE_CSV'
    elif args.summarize:
        logging.debug(f"Parameter indicate user wants to see summary of {args.summarize.name}")
        return 'SUMMARIZE'
    else:
        print("This script can export RefSeq to CCDS mappings, or it can combine two CSV files that contain mappings. The parameters you used do not make it clear which operation is desired.")
        print("To read a GFF and create CSV use: --gff GFF_FILE --source_release (105|105.20220307) --csv OUTPUT_CSV_FILE")
        print("To merge two CSV files use: --merge_old_csv OLD_CSV --merge_new_csv NEW_CSV --CSV OUTPUT_CSV_FILE")
        print("To print a summary of a CSV use: --summary --csv CSV_FILE")
        return None
        
def _create_mappings(source_release, gff_file, csv_file):
    '''
    Parse the GFF and write out all the RefSeq-to-CCDS mappings
    '''
    refSeqToCcds = RefseqToCcds()
    
    # Read mappings from gff
    refSeq_to_ccds_map = refSeqToCcds.get_mappings_from_gff(source_release, gff_file)
        
    # Write mappings to csv 
    refSeqToCcds.write_mappings(refSeq_to_ccds_map, csv_file)
        
    # Read mappings from csv for validation
    csv_mappings = refSeqToCcds.get_mappings_from_csv(csv_file)
    assert len(refSeq_to_ccds_map.keys()) == len(csv_mappings.keys()), f"The size of the GFF and CSV maps should be the same but {len(refSeq_to_ccds_map.keys())} != {len(csv_mappings.keys())}" 
     
def _combine_csv_files(older_csv_file, newer_csv_file, output_csv_file):
    '''
    Take all the mappings from two csv files and combine them into one. Sometimes the old source and the new source have different CCDS id for 
    a RefSeq id, when that happens the CCDS from the newer source is selected.    
    '''
    refSeqToCcds = RefseqToCcds()
    older_mapping = refSeqToCcds.get_mappings_from_csv(older_csv_file)
    newer_mapping = refSeqToCcds.get_mappings_from_csv(newer_csv_file)
    merged_mapping = refSeqToCcds._merge_maps(older_mapping, newer_mapping)
    refSeqToCcds.write_mappings(merged_mapping, output_csv_file)

def _print_summary(csv_file: str):
    refSeqToCcds = RefseqToCcds()
    d = refSeqToCcds.get_mappings_from_csv(csv_file)
    print(f"Number of RefSeq-to-CCDS mappings: {len(d)}")
    
    size_bytes = getsizeof(d)
    size_bytes += sum(map(getsizeof, d.values())) + sum(map(getsizeof, d.keys()))
    size_mb = round(size_bytes / (1024 * 1024), 2)
    print(f"Size of dictionary in memory: {size_mb}MB")

def _main():
    logging.config.dictConfig(TfxLogConfig().stdout_config)
    if logging.root.handlers[0].stream: 
        output = str(logging.root.handlers[0].stream.name)
    else:
        output = logging.root.handlers[0].baseFilename
        
    print(f"Log level={logging.root.getEffectiveLevel()}, output={output}")
    
    parser = ArgumentParser(description='', formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-c", "--csv", help="CSV file to write out RefSeq-to-CCDS mappings", type=FileType('w'), required=False)
 
    # These parameters are used when converting GFF to CSV
    parser.add_argument("-g", "--gff", help="GFF file from which to read RefSeq-to-CCDS mappings", type=FileType('r'), required=False)
    parser.add_argument("-r", "--source_release", help="Specify the GFF release version (found in the GFF's header annotation-source value)", action="store", required=False)
 
    # These parameters are used when merging two CSVs
    parser.add_argument("-m0", "--merge_old_csv", help="Previous revision of mappings in CSV", type=FileType('r'), required=False)
    parser.add_argument("-m1", "--merge_new_csv", help="", type=FileType('r'), required=False)

    # Display summary of a mapping CSV
    parser.add_argument('-s', '--summarize', help="Display summary of CSV mapping file", type=FileType('r'), required=False)
    
    # Version    
    parser.add_argument('-V', '--version', action='version', version=__version__)


    # Process arguments
    args = parser.parse_args()
    operation = _get_app_operation(args)
    
    if operation == 'CREATE_MAPPINGS':
        _create_mappings(args.source_release, args.gff.name, args.csv.name)
    elif operation == 'COMBINE_CSV':
        _combine_csv_files(args.merge_old_csv.name, args.merge_new_csv.name, args.csv.name)
    elif operation == 'SUMMARIZE':
        _print_summary(args.summarize.name)
    else:
        return 0
    
if __name__ == "__main__":    
    sys.exit(_main())
        
