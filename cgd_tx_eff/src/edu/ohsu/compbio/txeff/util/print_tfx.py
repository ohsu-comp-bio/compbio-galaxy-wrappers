'''
Created on Sep. 21, 2023

Read the transcript effects from vcf and render them in a text table like this: 

```
Variant: 16-89985662-GGACTAT-G
Splicing: False, Reference Context: TTCTCACAATTTGCAGAAACT
Tfx                      0                       1
-----------------------  ----------------------  ----------------------
TFX_BASE_POSITION        -5                      -5
TFX_EXON
TFX_GENE                 MC1R                    MC1R
TFX_G_DOT                g.89985662_89985667del  g.89985662_89985667del
TFX_HGVSC                c.-5_1del               c.-5_1del
TFX_HGVSP                p.Met1?                 p.Met1?
TFX_AMINO_ACID_POSITION
TFX_TRANSCRIPT           NM_002386.3             CCDS56011.1
TFX_VARIANT_TYPE         exonic                  exonic
TFX_PROTEIN_TRANSCRIPT   NP_002377.4             NP_002377.4
TFX_VARIANT_EFFECT       startloss               startloss
```

@author: pleyte
'''
from argparse import ArgumentParser
from tabulate import tabulate
import vcfpy

class PrintTfx(object):
    def print(self, vcf_filename: str, genotype:str = None):
        
        vcf_reader = vcfpy.Reader.from_path(vcf_filename)
        
        for row in vcf_reader:
            if not genotype or genotype == self._get_genotype(row):
                self._print(row)
    
    def _get_genotype(self, row: vcfpy.record.Record):
        return f'{row.CHROM}-{row.POS}-{row.REF}-{row.ALT[0].value}'
        
    def _print(self, row):
        '''
        '''
        splicing, ref_context, info_field_values = self._get_info_field_values(row)
        
        print(f"Variant: {self._get_genotype(row)}")
        print(f"Splicing: {splicing}, Reference Context: {ref_context}")
        
        if not info_field_values:
            # Not all variants have transcript effects.  
            print("No transcript effects found")
        else:        
            self._print_tfx(info_field_values)
        
        print("")

    def _get_info_field_values(self, row: vcfpy.record.Record):    
        '''
        Returns the transcript effects values.
        '''
        tfx_field_values = dict()
        tfx_splicing = None
        tfx_reference_context = None
        tfx_array_length = None
         
        for field, info_value in row.INFO.items():
            if not field.startswith('TFX'):
                continue
            
            # Unlike the other tfx fields these may be empty or only have a single value.             
            if field == 'TFX_SPLICE':
                value = self._get_single_value(info_value)
                if value:
                    tfx_splicing = value
                else:
                    tfx_splicing = False
                continue  
            elif field == 'TFX_REFERENCE_CONTEXT':
                tfx_reference_context = self._get_single_value(info_value)                
                continue

            # vcfpy wraps each info value in a list, but there is only ever one item in the list.
            assert len(info_value) == 1, "The INFO value returned by vcfpy should only have a single value"
            
            # transcript effects values are separated by a ':'. 
            tfx_field_values[field] = info_value[0].split(':')
            
            # Make sure all the delimited tfx fields (except splice and ref context) have the same number of values
            if not tfx_array_length:
                tfx_array_length = len(tfx_field_values[field])
            elif tfx_array_length != len(tfx_field_values[field]):
                raise ValueError(f'Variant {self.get_genotype(row)} field {field} does not have the expected number of values: {tfx_array_length}')
        
        return tfx_splicing, tfx_reference_context, tfx_field_values
          
    def _get_single_value(self, info_value):
        '''
        The value of an INFO field is returned as a list. This function extracts that value but returns
        None if the list is empty.  
        '''
        if len(info_value) == 0:
            return None
        elif len(info_value) == 1:
            return info_value[0]
  
    def _print_tfx(self, records: dict):
        '''
        Print a table of tfx values with labels (eg TFX_TRANSCRIPT) starting each row, and each column being 
        an index in the array of values. 
        '''
        # Number of values in the list for every key is the same.  
        array_length = len(next(iter(records.values())))
        
        # Column headings
        headers = ['Tfx'] + list(range(0,array_length))
        
        values_matrix = []

        for key, values in records.items():
            # The first column is the label, the rest are the list of values 
            tfx_row = [key] + values
            values_matrix.append(tfx_row)
            
        print(tabulate(values_matrix, headers=headers))
                
if __name__ == '__main__':
    parser = ArgumentParser(description='Read VCF with transcript effects and display the effects in a table.')
    parser.add_argument("-f", "--vcf", help="VCF file having transcript effects INFO fields.", type=str, required=True)
    parser.add_argument("-g", "--genotype", help="Limit output to variant matching the specified genotype.", type=str)
    
    args = parser.parse_args()
    
    print_tfx = PrintTfx()
    PrintTfx().print(args.vcf, args.genotype)
    
pass