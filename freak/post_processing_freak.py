import sys
import pandas as pd
import argparse

# Define function for filtering
def identity(freak, ref):
    return freak, ref

def remove_duplicates(freak, ref):
    '''
    result = freak.copy()
    # Inefficient, find a better way
    print('In remove dupes', file=sys.stderr)
    for i, row in ref.iterrows():
        chrm = row['Chromosome'][4:]
        start = row['Position Start']
        base = row['Variant Base']
        result = result[(result['rname']!=chrm)
                        | (result['start (1-based, inclusive)']!=start)
                        | (result['variant']==base)]
    '''
    result = freak.copy()
    locations = set([(row['rname'],row['start (1-based, inclusive)']) for i, row in freak.iterrows()])
    print(locations)
    for chrm, idx in locations:
        print(f'{chrm},{idx}')
        var = ref[(ref['Chromosome'].str.slice(start=3)==chrm) & (ref['Position Start']==idx)]['Variant Base'].values[0]
        print(var)
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
