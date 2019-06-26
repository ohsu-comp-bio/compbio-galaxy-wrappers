"""
    Filters the all variants file from MOON
"""

import sys
import argparse

def supply_args():
    """Populates the arguments"""
    #Initiate the parser
    parser = argparse.ArgumentParser(description='Filters and populates All Rare Variants file from Moon.')
    #Add the desired arguments
    #The all rare variants file
    parser.add_argument('variants', help='All Rare Variants file from Moon')
    #The optional gene list
    parser.add_argument('-geneList', help='Optional list of genes to further filter on.')
    #The data files for the populator
    parser.add_argument('-omim', help='OMIM data file.')
    parser.add_argument('-orphConv', help='Orphanet conversion file.')
    parser.add_argument('-orphPrev', help='Orphanet prevelence file.')
    #Filter off option
    parser.add_argument('-filterOff', action='store_const', const=True, help='Turns off the filter function.')
    #Gene header name
    parser.add_argument('-header', default='Gene', help='Header on the gene column of the file.')
    #create the object with the arguments as atributes
    args=parser.parse_args()
    return args

def build_array(tsv_file_name):
    """Takes the tsv file and turns it into an array"""
    #open the tsv file
    tsv_file = open(tsv_file_name, 'r')
    
    #start the array 
    tsv_array =  []
    
    #populate the array
    #loop over each line
    for line in tsv_file:
        #turn the line into an array by splitting the lines along the tabs
        array_line = line.split('\t')
        #add the line array to the tsv array
        tsv_array.append(array_line)
        
    #close the file
    tsv_file.close()
    
    return tsv_array

def find_index(header_line, header_name):
    """Find the index of a given header in a given header line"""
    header_array = header_line.split('\t')
    count = 0
    for entry in header_array:
        if entry.strip() == header_name:
            return count
        else:
            count += 1
    print('Could not find \"' + header_name + '\" column. Check to make sure your data is labeled correctly, or change the Gene Column Header under Non-Moon Specifications.')
    return -1

def filter_by_nums(array, index, max_num):
    """Filters an array by a numerical filter. takes an array, the index of the filter value, and the maximum passing value"""
    filtered_array = []
    for entry in array:
        try:
            entry_num = float(entry[index])
            if entry_num <= max_num:
                filtered_array.append(entry)
        except ValueError:
            pass
    return filtered_array
        
def filter_by_string(array, index, invalid_string):
    """Filters an array by a string type filter. Takes an array, the index of the filter value and the string that means the entry should be removed"""
    filtered_array = []
    for entry in array:
        entry_value = entry[index]
        if entry_value != invalid_string:
            filtered_array.append(entry)
    return filtered_array
            
def add_by_string(array_to, array_from, index, valid_string):
    """Adds entrys to an array if they have the valid string and removes the added ones from the original array"""
    for entry in array_from:
        entry_value = entry[index]
        if entry_value == valid_string:
            array_to.append(entry)
            array_from.remove(entry)

            
#Gene list filter
def build_list(file_name):
    """Build a list of genes from a txt file name with one gene per line"""
    #open the file
    file = open(file_name, 'r')
    #initiate list
    my_list = []
    #iterate over each line the file
    for line in file:
        #strip off the end line character
        entry = line.strip()
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
            #add the new splits back into the list
            string_list = string_list + new_genes
            

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
    
def filter_tsv_by_list(tsv_file_name, gene_list, gene_header = 'Gene'):
    """Takes a tsv and a list and writes a new tsv with only those tsv lines who include genes from the gene list"""
    #open the tsv
    tsv = open(tsv_file_name, 'r')
    #create the new tsv, with a name based on the old tsv name
    filtered_tsv_name = 'RelGenes.tsv'
    #tsv_file_name.split('.')[0] + '_RelGenes.tsv'
    filtered_tsv = open(filtered_tsv_name, 'w')
    #Get the first line from the old tsv as a header
    first_line = tsv.readline()
    #find the index of the GENES column
    #header_array = first_line.split('\t')
    #Find Gene index
    gene_index = find_index(first_line, gene_header)
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
            
def filter_moon_tsv(moon_tsv):
    """Filters the Moon tsv"""
    moon_file = open(moon_tsv, 'r')

    filtered_moon_file_name = 'filtered.tsv'
    #moon_tsv.split('.')[0] + '_filtered.tsv'
    filtered_moon_file = open(filtered_moon_file_name, 'w')
    
    header_line = moon_file.readline()
    filtered_moon_file.write(header_line)
    moon_file.close()
    
    
    moon_array = build_array(moon_tsv)
    #create array of variants we definately want to take a closer look at
    white_list_array = []
    
    #Add all the moon variants that were not filtered out
    moon_filter_index = find_index(header_line, 'Filter')
    add_by_string(white_list_array, moon_array, moon_filter_index, '')
    #Add all clinVar path and likely path
    clinVar_index = find_index(header_line, 'ClinVar')
    add_by_string(white_list_array, moon_array, clinVar_index, 'Pathogenic')
    add_by_string(white_list_array, moon_array, clinVar_index, 'Likely pathogenic')
    #Add all de novo variants
    de_novo_index = find_index(header_line, 'De Novo')
    add_by_string(white_list_array, moon_array, de_novo_index, 'TRUE')
    #Add all compound het variants
    comp_het_index = find_index(header_line, 'Compound Het.')
    add_by_string(white_list_array, moon_array, comp_het_index, 'TRUE')
    
    #For the variants not yet added to the white list filter out those we don't need to look at
    #gnomad frequency over 1%
    gnomad_freq_index = find_index(header_line, 'gnomAD frequency')
    moon_array = filter_by_nums(moon_array, gnomad_freq_index, .01)
    #diploid frequency over 5%
    diploid_freq_index = find_index(header_line, 'Diploid frequency')
    moon_array = filter_by_nums(moon_array, diploid_freq_index, .05)
    #More than 5 homozygotes in gnomad
    homozygotes_index = find_index(header_line, 'gnomAD homs')
    moon_array = filter_by_nums(moon_array, homozygotes_index, 5)
    #No associated disorder
    disorder_index = find_index(header_line, 'Disorder')
    moon_array = filter_by_string(moon_array, disorder_index, 'No associated disorder')
    #No coding or splicing effect and not a VUS on clinVar
    add_by_string(white_list_array, moon_array, clinVar_index, 'VUS')
    effect_index = find_index(header_line, 'Effect')
    moon_array = filter_by_string(moon_array, effect_index, '')
    
    #put the arrays back together
    total_array = white_list_array + moon_array
    #flag low read variants
    reads_index = find_index(header_line, 'Depth')
    
    for entry in total_array:
        count = 1
        
        for value in entry:
            if count < len(entry):
                filtered_moon_file.write(value + '\t')
                count += 1
            else:
                try:
                    if float(entry[reads_index]) < 10:
                        value = value.strip('\n') + '\t LOW READ DEPTH \n'
                except ValueError:
                    pass
                filtered_moon_file.write(value)
            
    filtered_moon_file.close()
    
    return filtered_moon_file_name
    
#Populator    
def build_dictionary(data_file_name, search_index, info_index, filter_on = False, filter_index = 32, filter_by = 'OMIM'):
    """Builds a dictionary of search term and info entries.
    Returns a dictionary.
    Parameters: File - tab delimited, index of the search term, index of the information, 
    Optional: turn on the filter to filter by the value in another term. filter on (true/false), index of the term to filter by, term needed to pass the filter.
    """
    #open the file with the needed data
    data_file = open(data_file_name, 'r')
    
    #Start the data dictionary
    data_dictionary = {}
    
    #populate the dictionary
    #loop over each entry
    for line in data_file:
        #split the entry into an array where each spot is one peice of information from that entry
        data_entry = line.split('\t')
        #check to be sure the entry has all the needed information
        if len(data_entry) > max(search_index, info_index):
            #extract the needed information from the array and assign to variables
            search = data_entry[search_index].strip()
            info = data_entry[info_index].strip()
            #.replace(',',';')
            #check if using the filter and the entry is long enough to use the filter
            if filter_on and len(data_entry) > filter_index:
                #get the info from the entry at the filter index
                filter_info = data_entry[filter_index].strip(' ')
                #if the filter does not pass then remove the info from the dictionary
                if filter_info != filter_by:
                    info = ''
            #check if the search term is already in the dictionary
            if search in data_dictionary:
                #if the search is in the dictionary then append the new value to the list and the info is not already in the dictionary then add it to the list
                if info not in data_dictionary[search] and info != '':
                    data_dictionary[search].append(info)
                    #check if the empty info is in the dictionary and remove it
                    if '' in data_dictionary[search]:
                        data_dictionary[search].remove('')
            #if the search term is not in the dictionary add it with its info
            else:
                    data_dictionary[search] = [info]
        #if the entry is not long enough to have both the search term and the info term then print a message to the command line so the problem entry can be examined 
        #ignore the comment lines and that one OMIM comment line is keeps poping up anyway. 
        elif data_entry[0].startswith('#') or data_entry[0].startswith('\xef\xbb\xbf# Copyright'):
            pass
        else:
            print("Problem with entry:")
            print(data_entry)
    #remove the empty entry from the dictionary        
    if '' in data_dictionary:
        del data_dictionary['']
    #close the file
    data_file.close()
    
    return data_dictionary

def add_info(search_array, info_dictionary, search_index):
    """adds information from the data dictionary to the array made from the tsv
        Parameters: two tiered array made from a tsv, dictionary with string keys and list values, index of the term to search on in the array/tsv
    """
    #loop over each entry in the array from the tsv file/ array
    #initiate count
    count = 0
    while count < len(search_array):
        #find the next entry
        entry = search_array[count]
        #check to make sure the entry is long enough to include the seach term
        if len(entry) > search_index:
            #identify the searching item in question
            search = entry[search_index].strip(' ')
        else:
            #if the entry is not long enough then give a dumy search term
            search = "None"
        #check that the location is in the dictionary
        if search in info_dictionary:
            #retrieve the desired information from the dictionary
            info = info_dictionary[search]
            #loop over each of the entries in info, popping them out after each one
            while len(info) > 1:
                #loop through the info array and add new line to the array
                #find the index of the current entry and add one to get just below
                index = count + 1
                #form new a line just below the current entry with the current peice of info
                search_array.insert(index, [info.pop(0), '\n'])
                #advance the count of the search array so that the new line is skipped
                count = count + 1
            #with only one entry left add it to the entry in-line    
            entry.insert(0, info[0])
        #if the entry is deemed not a real entry then add a new empty column to the array to keep things orderly
        else:
            entry.insert(0, '')
        #advance the count to the next entry line
        count = count + 1
    
def write_tsv_from_array(two_tier_array, file_name):
    """creates a tsv from a two tiered array
        Parameters: two tiered array, name for new tsv file
    """
    #create a new file
    new_tsv = open(file_name, 'w')
    
    #loop over each line in the array
    for line in two_tier_array:
        #loop over each entry in the line
        for entry in line:
            #for each entry write it to the new file and add a comma at the end
            if line.index(entry) < len(line)-1:
                new_tsv.write(entry + '\t')
            #for the final entry don't add a comma to the end.
            else:
                new_tsv.write(entry)
        
    #close the file
    new_tsv.close()
    
    return file_name

def add_data_from(array_name, data_file_name, data_search_index, tsv_search_index, data_to_add_index, filter_on = False, filter_index = 32, filter_by = 'OMIM'):
    """Returns a tsv like the given tsv but with new data from a given data file
        parameters: tsv file name that needs data added
            tab delimited file name with the desired information
            index of the search term in the data file
            index of the search term in the tsv
            index of the data in data file
            Optional: if the data dictionary needs to be filtered then:
            filter on (true/false)
            index of the filter term in the data file
            term needed to pass the filter
    """
    #build the dictionary
    data_dict = build_dictionary(data_file_name, data_search_index, data_to_add_index, filter_on, filter_index, filter_by)
    #add the infomation from the dictionary to the array
    add_info(array_name, data_dict, tsv_search_index)
    
    #return the new tsv file name
    return array_name
    
def clean(tsv_array, out_put_file_name):
    """takes a tsv file created by the populate function and returns a tsv file that has column headers"""
    
    #add the headers 
    tsv_array[0][0] = 'Inheritance/Disorder: OMIM'
    tsv_array[0][1] = 'Prevelance 2016'
    tsv_array[0][2] = 'Disorder: Orphanet'
    tsv_array[0][3] = 'Orphanum'
   
    #turn the array back into a tsv named the input name        
    new_tsv = write_tsv_from_array(tsv_array, out_put_file_name)
    #return the new tsv
    return new_tsv
    
def populate(tsv_file_name, OMIM_data, Orphadata_1, Orphadata_2, gene_header = 'Gene'):
    """Takes a tsv file and populates it from OMIM's genemap2.txt and orphanet 
        adds: orphanet number based on the gene symbol
        adds: disease, OMIM number and prevelace based on orphanet number
        adds: OMIM inheritance based on OMIM number
        removes entries with only OMIM numbers to clean up the tsv
    """
    #Gene to Orphanum
    #Orphanum to Disease
    #Orphanum to prevelance
    #Gene to Inheritance/OMIM Disorder
    
    
    #CHECK THESE NUMBERS AND DATA FILES FIRST IN CASE OF NOT WORKING
    
    #Indecies
    
    #orphanet
    #gene to disease database
    gene_orphanet = 0
    orphanum = 2
    #omim_number_orphanet = 33
    disease_orphanet = 1
    #omim_validation_orphanet = 32
    #prevelance database
    prevelance_orphanet = 27
    validation_orphanet = 36
    orphanum_prev = 5
    
    #OMIM
    omim_number_omim = 5
    inheritance_omim = 12
    OMIM_gene = 8
    
    #tsv
    #gene_tsv = gene_index
    
    #data files
    #Orphanet cross referencer. Gene to disorder file 
    #Orphadata_1 = 'orphanet_data_symbol_to_disorder.txt'
    #Orphanet prevelance. 
    #Orphadata_2 = 'orphanet_data_1.txt'
    #OMIM data. genemap2.txt
    #OMIM_data = 'genemap2.txt'
    
    
    array_1 = build_array(tsv_file_name)
    #get the gene index 
    first_line = ''
    for entry in array_1[0]:
        first_line = first_line + entry + ('\t')
    gene_tsv = find_index(first_line, gene_header)
    #get the Orphanet number from the gene symbol using the orphanet gene to disease reference
    array_2 = add_data_from(array_1, Orphadata_1, gene_orphanet, gene_tsv, orphanum)
    #update the tsv indecies
    orphanum_tsv = 0
    gene_tsv += 1
    #get the disease from the orphanet number using the orphanet gene to disease reference
    array_3 = add_data_from(array_2, Orphadata_1, orphanum, orphanum_tsv, disease_orphanet)
    #update the needed tsv indecies
    orphanum_tsv += 1
    gene_tsv += 1
    #get the prevelance from the orphanet number using the Orphanet database, filter to make sure only validated prevelance data is retrieved
    array_4 = add_data_from(array_3, Orphadata_2, orphanum_prev, orphanum_tsv, prevelance_orphanet, filter_on = True, filter_index = validation_orphanet, filter_by = 'Validated')
    #update needed indices
    gene_tsv += 1
    #get the inheritance/disorder from the gene using the OMIM data
    array_5 = add_data_from(array_4, OMIM_data, OMIM_gene, gene_tsv, inheritance_omim)
        
    #Create a name for the output tsv by adding 'populated' to the end of the name

    new_tsv_name = 'populated.tsv'
    #tsv_file_name.split('.')[0] + '_populated.tsv'
    clean(array_5, new_tsv_name)
           
"""    
def filter_manual(moon_tsv):
    #Seeks user input for filtering a moon file
    moon_file = open(moon_tsv, 'r')
    filtered_moon_file_name = moon_tsv.split('.')[0] + '_filtered.tsv'
    filtered_moon_file = open(filtered_moon_file_name, 'w')
    
    header_line = moon_file.readline()
    filtered_moon_file.write(header_line)
    moon_file.close()
    
    moon_array = build_array(moon_tsv)
    #create array of variants we definately want to take a closer look at
    white_list_array = []
    
    #initiate a while counter
    again = True
    while again:
        mode = input('Filter type (Add/Remove/Under):')
        if mode.lower() == 'add':
            header = input('What is the header on the colunm you want to search by?')
            valid_string = input('What is the entry type that you want to add?')
            index = find_index(header_line, header)
            add_by_string(white_list_array, moon_array, index, valid_string)
        elif mode.lower() == 'remove':
            header = input('What is the header on the colunm you want to search by?')
            invalid_string = input('What is the entry type that you want to remove?')
            index = find_index(header_line, header)
            moon_array = filter_by_string(moon_array, index, invalid_string)
        elif mode.lower() == 'under':
            header = input('What is the header on the colunm you want to filter on?')
            max_num_str = input('What is the maximum value you want to keep?')
            max_num = float(max_num_str)
            index = find_index(header_line, header)
            moon_array = filter_by_nums(moon_array, index, max_num)
        else:
            print('To pick entries you want to keep regardless of future filters that may remove them use ADD. \nTo remove entries based on a word colunm use REMOVE, for instance remove all entries with \'No associated disorder\' under Disorder. \nTo remove entries based on a number type column use UNDER, for instance you want to keep all entries with a diploid frequency under .5')
        another = input('Do you want to add another filter? (Yes/No)')
        if another.lower() == 'no':
            again = False
        elif another.lower() == 'yes':
            again = True
        else:
            again = False
            confirm = input('That was not a Yes or No, so it will be taken as a No. \nIf you would like to add another filter please type Yes now.')
            if confirm.lower() == 'yes':
                again = True
    #put the arrays back together
    total_array = white_list_array + moon_array
    #flag low read variants
    reads_index = find_index(header_line, 'Depth')
    
    for entry in total_array:
        count = 1
        
        for value in entry:
            if count < len(entry):
                filtered_moon_file.write(value + '\t')
                count += 1
            else:
                try:
                    if float(entry[reads_index]) < 10:
                        value = value.strip('\n') + '\t LOW READ DEPTH \n'
                except ValueError:
                    pass
                filtered_moon_file.write(value)
            
    filtered_moon_file.close()
    
    return filtered_moon_file_name        
            
            
"""    
    
def main():
    args = supply_args()
    tsv_file = args.variants
    moon_filtered = tsv_file
    if not args.filterOff:
        moon_filtered = filter_moon_tsv(tsv_file)
    if args.geneList:
        gene_list = build_list(args.geneList)
        moon_filtered = filter_tsv_by_list(moon_filtered, gene_list, args.header)
    if args.omim and args.orphConv and args.orphPrev:
        omim_data = args.omim
        orphanet_conv_data = args.orphConv
        orphanet_prev_data = args.orphPrev
        populate(moon_filtered, args.omim, args.orphConv, args.orphPrev, args.header)
    """    
    if len(args) == 1:
        file_name = input('File name: ')
        mode = input('Mode (Filter/Gene/Populate): ')
        if mode.lower() == 'gene':
            gene_list_file = input('Gene list:')
            gene_header = input('Header of Gene Column: ')
            gene_list = build_list(gene_list_file)
            filter_tsv_by_list(file_name, gene_list, gene_header)
        elif mode.lower() == 'populate':
            gene_index_str = input('Index of Gene Column: \n(Column A has index 0 \nMoon index is 7, Galaxy Filtered Variant files index is 0)\n')
            gene_index = int(gene_index_str)
            populate(file_name, gene_index)
        elif mode.lower() == 'filter':
            filter_manual(file_name)
            
    else:        
        tsv_file = args[1]
        moon_filtered = filter_moon_tsv(tsv_file)
        if len(args) > 2:
            gene_list = build_list(args[2])
            moon_filtered = filter_tsv_by_list(moon_filtered, gene_list)
        populate(moon_filtered)
    """
    
if __name__ == "__main__":
    main()