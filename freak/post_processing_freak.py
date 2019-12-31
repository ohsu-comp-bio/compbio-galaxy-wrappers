#!/usr/bin/python3
import sys
import pandas as pd
import argparse

'''
Post processing functions

Should follow the format shown in the 'identity' function, where the Freak 
output file and reference file and taken as inputs and returned as an output 
tuple. Each post processing function may make alterations to either DataFrame 
before returning, but the effects of these operations should be considered 
when determining the order in which they should be invoked.
'''
def identity(freak, ref):
    '''
    An "identity" function which return the Freak output and reference file 
    without altering them. Used as a default function in the argparser.
    '''
    return freak, ref

def remove_duplicates(freak, ref):
    '''
    Removes rows from the Freak output file which have bases that do not match 
    the reference file, returns the modified Freak file and the reference as 
    a tuple
    '''
    result = freak.copy()
    # Create an array of tuples with the chromosome and start index from the 
    # freak output file
    locations = set([(row['rname'],row['start (1-based, inclusive)']) for i, row in freak.iterrows()])
    # Interate through the array...
    for chrm, idx in locations:
        # Find the variant base specified in the reference file
        var = ref[(ref['Chromosome'].str.slice(start=3)==chrm) 
                  & (ref['Position Start']==idx)]['Variant Base'].values[0]
        # Select the rows that either do not occur as this location, or which 
        # have the variant base from the reference file
        result = result[(result['variant']==var)
                        | ~(result['start (1-based, inclusive)']==idx)
                        | ~(result['rname']==chrm)]
    return result, ref

# Instantiate the argument parser
parser = argparse.ArgumentParser(description='Post processing for freak output')

# Add arguments to the arg parser, and then parse
parser.add_argument('freak_out', action='store')
parser.add_argument('variant_ref', action='store')
parser.add_argument('-d','--no-dups',action='store_const',const=remove_duplicates,default=identity,dest='do_dupes')
parser.add_argument('-o','--output',action='store',default=None,dest='output')
args = parser.parse_args()

# Read in Freak's output and the reference 
freak_out = pd.read_csv(args.freak_out, sep='\t', header=0)
variant_ref = pd.read_csv(args.variant_ref, sep='\t', header=0)


# Apply all of the appropriate filters 
freak_out, variant_ref = args.do_dupes(freak_out, variant_ref)


# Output the results of the post processing
if args.output:
    out = open(args.output, mode='w')
else:
    out = sys.stdout
print(freak_out.to_csv(sep='\t',index=False), file=out)
