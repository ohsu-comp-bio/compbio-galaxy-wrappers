""" gene_filter.py
    Script that takes a text file with a list of genes and a tsv file with CNVs and returns a tsv file with a list of CNVs that contain at least one of the genes. 
    Eleanor Campbell
"""

import sys

def build_list(file_name):
    """Build a list of genes from a txt file name with one gene per line"""
    #open the file
    file = open(file_name, 'r')
    #initiate list
    my_list = []
    #iterate over each line the file
    for line in file:
        #strip off the end line character
        entry = line.strip().upper()
        #strip off any extra white space and turn all letters uppercase to make it not case sensitive
        entry = entry.strip(' ').upper()
        #check if the gene is now empty
        if entry != '':
            #add the clean uppcase gene to the list
            my_list.append(entry)
    #sort the list
    my_list.sort()
    #close the file
    file.close()
    #return the list
    return my_list

def new_split(string_list, delimeter):
    """Split an already split string again by a new delimeter, takes a list (the output of the first split) and a new delimeter"""
    #Since the genes are also split by semi colon, check each 'gene' to see if it has a semi colon in it
    for gene in string_list:
        if delimeter in gene:
            #if it has a semicolon then remove that string from the list
            string_list.remove(gene)
            #split the string by the semicolon
            new_genes = gene.split(delimeter)
            #add each of the new splits back into the list
            for new_gene in new_genes:
                string_list.append(new_gene)

def gene_in_line(gene_list, line, gene_index):
    """returns true or false depending on if any of the genes in this line are in the gene list. parameters: gene_list of type list, line of type string, and gene_index, where the gene index is the position of the gene list in the TSV"""
    #turn the line uppercase to make the search non-case sensitive
    line_upper = line.upper()
    #split the line by tab and take the entry at the gene index
    line_genes = line_upper.split('\t')[gene_index]
    #split the line by commas to creat a list of genes
    line_genes_list = line_genes.split(',')
    #now split on the semicolon
    new_split(line_genes_list, ';')
    #check each gene in the line's gene list to see if it is in the given gene list
    for gene in line_genes_list:
        #strip down the gene to the essentials
        gene = gene.strip('\n')
        gene = gene.strip(' ')
        gene = gene.strip('"')
        #check if the gene is in the gene list 
        if gene in gene_list:
            #if the gene is in the gene list then simply return true
            return True
    #if the each of the genes are check and none force a true return then return False
    return False
    
def filter_tsv_by_list(tsv_file_name, gene_list):
    """Takes a tsv and a list and writes a new tsv with only those tsv lines who include genes from the gene list"""
    #open the tsv
    tsv = open(tsv_file_name, 'r')
    #create the new tsv, with a name based on the old tsv name
    filtered_tsv_name = 'filtered.tsv'
    filtered_tsv = open(filtered_tsv_name, 'w')
    #Get the first line from the old tsv as a header
    first_line = tsv.readline()
    #find the index of the GENES column
    header_array = first_line.split('\t')
    #initiate the count
    count = 0
    #initiate a gene index
    gene_index = 0
    #loop over each of the elements in the header to find GENES
    for header in header_array:
        #strip off the new line character in case it is there
        header = header_array[count].strip('\n')
        header = header.strip('"')
        #check if the header for this column is GENES
        if header == 'GENE':
            #if the header name is GENES then make the gene_index equal to this count
            gene_index = count
        #move the count forward    
        count = count + 1    
    #Add the first line from the old tsv to the new tsv as a header    
    filtered_tsv.write(first_line)
    #Loop over each line in the old tsv
    for line in tsv:
        #check if the line has a valid gene in its gene list
        if gene_in_line(gene_list, line, gene_index):
            #if so add it to the new tsv
            filtered_tsv.write(line)
    #close the files
    tsv.close()
    filtered_tsv.close()
    
    
def main(args):
    if len(args) == 3:
        gene_list = build_list(args[2])
    elif len(args) == 4:
        list_1 = build_list(args[2])
        list_2 = build_list(args[3])
        gene_list = list_1 + list_2
    else:
        gene_list = []
        print("Please provide an appropriate number of arguments: \npython gene_filter.py cnv_file.tsv gene_list.txt \nor \npython gene_filter.py cnv_file.tsv gene_list_1.txt gene_list_2.txt")
    filter_tsv_by_list(args[1], gene_list)
    print(gene_list)
    
        
if __name__ == "__main__":
    args = sys.argv
    main(args)
    
    
