""" de_dup.py
    Script that takes two text files with lists of genes and returns a text file with a list of the genes from the medical list with duplicates removed
    Eleanor Campbell
"""

import sys

def build_list(file_name):
    
    file = open(file_name, 'r')
    
    my_list = []
    
    for line in file:
        entry = line.strip()
        entry = entry.strip(' ').upper()
        my_list.append(entry)
        
    my_list.sort()

    file.close()
    
    return my_list

def de_dup(med_list, hpo_list):
    
    new_list = []
    
    for gene in med_list:
        if gene not in hpo_list:
            new_list.append(gene)
    
    return new_list
            
    
def write_list(my_list, new_name):
    
    file = open(new_name, 'w')
    
    for entry in my_list:
        file.write(entry + '\n')
        
    file.close()
    
    
def main(med_genes, hpo_genes):
    med_gene_list = build_list(med_genes)
    hpo_gene_list = build_list(hpo_genes)
    de_duped_list = de_dup(med_gene_list, hpo_gene_list)
    
    new_name = 'output.txt'
    #hpo_genes.split('.')[0] + '_de_duped.txt'
    
    write_list(de_duped_list, new_name)
        
if __name__ == "__main__":
    args = sys.argv
    main(args[1], args[2])
    