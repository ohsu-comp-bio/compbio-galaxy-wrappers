""" gene_filter_v2.py
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
        #entry = entry.strip(' ').upper()
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
        gene = gene.strip()
        #gene = gene.strip(' ')
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
    #tsv_file_name.split('.')[0] + '_filtered.tsv'
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
        header = header_array[count].strip()
        #check if the header for this column is GENES
        if header == 'GENES':
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
    #return the new file name
    return filtered_tsv_name

def find_frequency(line, type_index, del_freq_index, dup_freq_index):
    """Find the frequency of the variant given the variant line in the CNV file, the index of the type of CNV, the index of the deletion frequency and the index of the duplication frequency """
    #Split the line by tab to get each column alone
    entry = line.split('\t')
    #get the entry type with the type index
    entry_type = entry[type_index]
    if entry_type == 'DUP':
        #if the entry is a duplication then get the frequency with the duplication index
        frequency = entry[dup_freq_index]
    elif entry_type == 'DEL':
        #if the entry is a deletion then get the frequency with the deletion index
        frequency = entry[del_freq_index]
    else:
        #if the type is not one the given types then return a frequency of -1 so that it will be called and looked a more closely
        frequency = -1        
    return frequency

def get_index(header_line, header_name):
    """gets the index of the given header from a the header line"""
    #Split the hedaer line into its individual headers
    entry = header_line.split('\t')
    #initiate an index count
    count = 0
    #loop over each header in the line
    for header in entry:
        #if the header matches the desired headername then return the count
        if header.strip() == header_name:
            return count
        else:
            #increase the count if the header is not a match
            count = count + 1
            
def filter_by_freq(tsv_name):
    """Filter a tsv by frequency and write two tsvs one with a fitler of .05 and one with a fitler of .01"""
    #open the tsv to be filtered
    tsv = open(tsv_name, 'r')
    #create names for the new tsvs by adding _u1 and _u5 to the end of the original tsv name
    filtered_5_tsv_name = 'u5.tsv'
    #tsv_name.split('.')[0] + '_u5.tsv'
    filtered_1_tsv_name = 'u1.tsv'
    #tsv_name.split('.')[0] + '_u1.tsv'
    #create the new tsvs
    tsv_5 = open(filtered_5_tsv_name, 'w')
    tsv_1 = open(filtered_1_tsv_name, 'w')
    
    #get the header line from the original and write it to the new tsvs
    first_line = tsv.readline()
    tsv_1.write(first_line)
    tsv_5.write(first_line)
    #get the relevant indecies from the hedaer line
    type_index = get_index(first_line, 'DEL/DUP')
    del_freq_index = get_index(first_line, 'POP DEL AF')
    dup_freq_index = get_index(first_line, 'POP DUP AF')
    #loop over the lines in the tsv
    for line in tsv:
        #get the frequency
        frequency = float(find_frequency(line, type_index, del_freq_index, dup_freq_index))
        #if the frequency is small enough then write it to the right files
        if frequency <= .01:
            tsv_1.write(line)
            tsv_5.write(line)
        elif frequency <= .05:
            tsv_5.write(line)
    #close the files        
    tsv.close()
    tsv_1.close()
    tsv_5.close()
 
    
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
    new_tsv_name = filter_tsv_by_list(args[1], gene_list)
    filter_by_freq(new_tsv_name)
        
if __name__ == "__main__":
    args = sys.argv
    main(args)
    
    