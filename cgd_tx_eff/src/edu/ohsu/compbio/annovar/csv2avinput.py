'''
Read a CSV file and write out an Annovar input file. In an Annovar input file "The first five space or tab delimited fields are Chromosome ("chr" prefix is optional), 
Start, End, Reference Allele, Alternative Allele. The rest of the columns are completely optional." (See https://annovar.openbioinformatics.org/en/latest/user-guide/input)
Created on Jul 11, 2022

@author: pleyte
'''
import argparse
import logging.config
import csv
from edu.ohsu.compbio.txeff.util.tfx_log_config import TfxLogConfig

VERSION = '0.0.1'

def _parse_args():
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='Read a CSV and write out a tab-delimited Annovar input file')

    parser.add_argument('-c', '--csv',  
                        help='Input CSV with genotypes',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('-o', '--out_file',  
                        help='Annovar "avinput" file to create',
                        type=argparse.FileType('w'),
                        required=True)    
    
    args = parser.parse_args()
    
    return args

def _main():
    '''
    main function
    '''
    logging.config.dictConfig(TfxLogConfig().utility_config)
    
    args = _parse_args()

    with open(args.csv.name) as csv_file, open(args.out_file.name, "w") as tsv_out_file:
        csv_reader = csv.DictReader(csv_file)
        tsv_writer = csv.writer(tsv_out_file, delimiter = '\t')
        
        duplicates = 0
        total = 0 
        
        keys = set()
        key_maker = lambda x: "-".join([x['chromosome'], str(x['position_start']), str(x['position_end']), x['reference_base'], x['variant_base']])
        
        rowId = 1
        for row in csv_reader:
            key = key_maker(row)
            if key in keys:
                logging.debug(f"Skipping duplicate variant {key}")
                duplicates = duplicates + 1
            else:
                tsv_writer.writerow([row['chromosome'], row['position_start'], row['position_end'], row['reference_base'], row['variant_base'],
                                     row['chromosome'], row['position_start'], rowId, row['reference_base'], row['variant_base'],
                                     "1.0", "PON", "x=y", "GT", "0/0", ''])
                keys.add(key)
                rowId = rowId + 1
   
            total = total + 1
            
    logging.info(f"duplicate count={duplicates}")
    logging.info(f"unique count={rowId}")
    logging.info(f"total={total}")
    logging.info(f"Wrote {args.out_file.name}")
            
if __name__ == '__main__':
    _main()