'''
Created on Sep. 21, 2023

Read the transcript effects from vcf and render them in a text table like this: 

```
Variant: 16-89985662-GGACTAT-G
Tfx                      0                       1
-----------------------  ----------------------  ----------------------
TFX_BASE_POSITION        -5                      -5
TFX_EXON
TFX_GENE                 MC1R                    MC1R
TFX_G_DOT                g.89985662_89985667del  g.89985662_89985667del
TFX_HGVSC                c.-5_1del               c.-5_1del
TFX_HGVSP                p.Met1?                 p.Met1?
TFX_SPLICE
TFX_AMINO_ACID_POSITION
TFX_TRANSCRIPT           NM_002386.3             CCDS56011.1
TFX_VARIANT_TYPE         exonic                  exonic
TFX_PROTEIN_TRANSCRIPT   NP_002377.4             NP_002377.4
TFX_VARIANT_EFFECT       startloss               startloss
```

This script is a little awkward to run because it takes as it's parameter the line from a vcf. The best way to execute it is: 
```sh
grep -v '#' tfx.vcf | while read -r line ; do
  python $TFX_HOME/src/edu/ohsu/compbio/txeff/util/print_tfx.py --tfx "$line"
done
```


@author: pleyte
'''
import re
from argparse import ArgumentParser
from tabulate import tabulate

class PrintTfx(object):
    def print(self, vcf_line:str):
        genotype = self._get_genotype(vcf_line)
        info = self._get_info(vcf_line)
        records = self._parse(info)
        
        print(f"Variant: {genotype}")
        self._print_tfx(records)
        
    def _print_tfx(self, records: dict):
        '''
        Print a table of tfx values with labels (eg TFX_TRANSCRIPT) starting each row, and each column being 
        an index in the array of values. 
        '''
        array_length = len(next(iter(records.values())))
        headers = ['Tfx'] + list(range(0,array_length))
        values_matrix = []

        for key, values in records.items():
            tfx_row = [key] + values
            values_matrix.append(tfx_row)
            
        print(tabulate(values_matrix, headers=headers))
        
    def _parse(self, info:str) -> dict:
        '''
        Retrieve the TFX labels and value arrays from a vcf line and return them in a dict object
        ''' 
        
        info_fields = info.split(';')
        array_size = 0
        type_values = dict()
        
        for info in info_fields:
            if(info.startswith('TFX')):
                if(info == 'TFX_VARIANT_EFFECT'):
                    print(123)
                
                label, values = self._parse_individual(info)
                
                if array_size == 0:
                    array_size = len(values)
                elif len(values) != array_size:
                    raise ValueError(f"Inconsistent number of values found in {label} len={len(values)}, expected={array_size}")
                    
                type_values[label] = values
        
        return type_values
                
    def _parse_individual(self, section: str) -> (str, list):
        '''
        Return the TFX type and an array of values from a string like "TFX_BASE_POSITION=951:951:951:951"
        '''
        label = re.search('TFX_[A-Z_]+', section).group(0)
        delimited_values = section.split("=",1)[1]
        return label, delimited_values.split(':')
        
    def _get_genotype(self, vcf_line:str):
        '''
        Return the genotype defined in a VCF line
        '''
        sections = vcf_line.split("\t")
        chromosome = sections[0]
        pos = sections[1]
        ref = sections[3]
        alt = sections[4]
        return f"{chromosome}-{pos}-{ref}-{alt}"
    
    def _get_info(self, vcf_line):
        '''
        Return the INFO section of the VCF line. 
        '''
        sections = vcf_line.split("\t")
        return sections[7]
        
        
if __name__ == '__main__':
    parser = ArgumentParser(description='Accepts the delimited TFX string that is found in the INFO section of a VCF and prints it to screen in a table.')
    parser.add_argument("-t", "--tfx", help="Delimited transcript effects string", type=str, required=True)
    
    args = parser.parse_args()
    PrintTfx().print(args.tfx)
    
pass