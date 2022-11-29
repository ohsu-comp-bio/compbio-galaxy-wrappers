'''
Read a CSV file and write out an Annovar input file. In an Annovar input file "The first five space or tab delimited fields are Chromosome ("chr" prefix is optional), 
Start, End, Reference Allele, Alternative Allele. The rest of the columns are completely optional." (See https://annovar.openbioinformatics.org/en/latest/user-guide/input)
Created on Jul 11, 2022

@author: pleyte
'''
import argparse
import logging
import csv

VERSION = '0.0.1'

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

stream_handler = logging.StreamHandler()
logging_format = '%(levelname)s: [%(filename)s:%(lineno)s - %(funcName)s()]: %(message)s'

stream_format = logging.Formatter(logging_format)
stream_handler.setFormatter(stream_format)
stream_handler.setLevel(logging.DEBUG)
logger.addHandler(stream_handler)


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
                logger.debug(f"Skipping duplicate variant {key}")
                duplicates = duplicates + 1
            else:
                tsv_writer.writerow([row['chromosome'], row['position_start'], row['position_end'], row['reference_base'], row['variant_base'],
                                     row['chromosome'], row['position_start'], rowId, row['reference_base'], row['variant_base'],
                                     "1.0", "PON", "x=y", "GT", "0/0", ''])
                keys.add(key)
                rowId = rowId + 1
   
            total = total + 1
            
    logger.info(f"duplicate count={duplicates}")
    logger.info(f"unique count={rowId}")
    logger.info(f"total={total}")
    logger.info(f"Wrote {args.out_file.name}")
            
if __name__ == '__main__':
    _main()