'''
Created on Sep. 21, 2023

Read the transcript effects from json and render them in a text table like this: 

```
Variant: 1-11026381-TAA-T
Tfx                0                        1                        2
-----------------  -----------------------  -----------------------  -----------------------
id                 2654332                  2654332                  2654332
chromosome         1                        1                        1
position           11026381                 11026381                 11026381
reference          TAA                      TAA                      TAA
alt                T                        T                        T
variantEffect
variantType        intronic                 intronic                 intronic
aminoAcidPosition
basePosition       372                      -46                      372
exon
gene               C1orf127                 C1orf127                 C1orf127
cDot               c.372+20_372+21del       c.-46+20_-46+21del       c.372+20_372+21del
pDot1              p.?                      p.?                      p.?
pDot3              p.?                      p.?                      p.?
splicing
cdnaTranscript     NM_001170754.1           NM_001366227.1           CCDS53267.1
proteinTranscript  NP_001164225.1           NP_001353156.1           NP_001164225.1
sequenceVariant    g.11026394_11026395del   g.11026394_11026395del   g.11026394_11026395del
referenceContext   AATCTCCCTTTAAAAAAAAAAAA  AATCTCCCTTTAAAAAAAAAAAA  AATCTCCCTTTAAAAAAAAAAAA
```

@author: pleyte
'''
from _collections import defaultdict
from argparse import ArgumentParser
import argparse
import json

from tabulate import tabulate


class PrintTfx(object):
    def print(self, in_stream, genotype:str = None):
        variant_transcripts = defaultdict(list)
        
        for rec in json.load(in_stream):
            if not genotype or genotype == self._get_genotype(rec):
                variant_transcripts[self._get_genotype(rec)].append(rec)
        
        for variant, transcripts in variant_transcripts.items():
            self._print(variant, transcripts)
            print("")

    def _get_genotype(self, row: dict):
        return f"{row['chromosome']}-{row['position']}-{row['reference']}-{row['alt']}"
        
    def _print(self, variant: str, transcripts: list):
        '''
        Transcripts is a list of dict objects where each dict is a transcript read from the input json  
        '''        
        print(f"Variant: {variant}")
        self._print_tfx(transcripts)        
        
    def _print_tfx(self, transcripts: list):
        '''
        Display the list of transcripts in a table
        '''
        # Column headings
        headers = ['Tfx'] + list(range(0,len(transcripts)))
        
        values_matrix = []

        for key in transcripts[0].keys():
            tfx_row = [key] + [x.get(key) for x in transcripts]
            values_matrix.append(tfx_row)
            
        print(tabulate(values_matrix, headers=headers))

if __name__ == '__main__':
    parser = ArgumentParser(description='Read json with transcript effects and display the effects in a table.')
    parser.add_argument("--in",
                        dest="input",
                        type=argparse.FileType('r'),
                        required=True,
                        help="json file having transcript effects")
    
    parser.add_argument("-g", "--genotype", help="Limit output to a specific variant (eg 1-123-A-C).", type=str)
    
    args = parser.parse_args()
    
    print_tfx = PrintTfx()
    PrintTfx().print(args.input, args.genotype)
    
pass