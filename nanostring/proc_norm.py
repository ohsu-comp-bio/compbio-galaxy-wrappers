#!/usr/bin/env python

## USAGE:
#   Nanostring output RCC readouts (html format) to raw data matrix
## NOTES:
#   TODO: break out tests
#       checking sample metrics
#       check for unique sample names in sample sheet
#
## ARGUMENTS:
#   samplesheet: tsv of RCC filename in one col and sample name, 'ANTIBODY_REFERENCE.csv'
#   rcc files
#   abfile: antibody file to rename antibody names.
## RETURN:
#   rawdata.txt: tab-sep file with table of Antibody x Sample
#
import os
import re
import argparse
import pandas
import xml.etree.ElementTree as ET
#norm_geomean import statements
import math
import csv
import copy
from functools import reduce

VERSION="1.0.0"

#Global variables for norm_geomean
controls = ["MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468control", "MDA468+EGF"]
mouseAb = ["Ki-67", "Pan-Keratin", "S6 Ribosomal", "p53", "MmAb-IgG1"]

#make the control help statement for norm geomean ctrl input
control_help = "Marks all samples to be used as control. Otherwise controls are "
for sample in controls:
    control_help += sample

def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='Nanostring output RCC readouts (html format) and '
                                                 'converts to raw data tsv')
    parser.add_argument('rcc_files', type=str, nargs='+', help='raw RCC files')
    parser.add_argument('--samplesheet', type=argparse.FileType('r'), help='samplesheet.txt')
    parser.add_argument('--abfile', nargs='?', type=argparse.FileType('r'),
                        help='ANTIBODY_REFERENCE.csv')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    #For norm_geomean
    parser.add_argument('-ctrl', action="store_true", help=control_help)
    args = parser.parse_args()
    return args

def parseRCC(rcc_file):
    """
    Take RCC file output df
    Args:
        path to RCC
    Returns:
        dataframe
    Examples:
    """
    dict_rcc = {}
    with open(rcc_file, 'r') as xml_handle:
        complete_xml = "<root>" + str(xml_handle.read()) + "</root>"
        root = ET.fromstring(complete_xml)
        for child in root:
            dict_rcc[child.tag] = child.text

        counts = [line.strip().split(',') for line in dict_rcc['Code_Summary'].strip().split('\n')]
        header = [line.strip().split(',') for line in dict_rcc['Header'].strip().split('\n')]
        samp_attrib = [line.strip().split(',') for line in
                       dict_rcc['Sample_Attributes'].strip().split('\n')]
        lane_attrib = [line.strip().split(',') for line in
                       dict_rcc['Lane_Attributes'].strip().split('\n')]

        dfcounts = pandas.DataFrame(counts[1:len(counts)], columns=counts[0])
        df_samp_attrib = pandas.DataFrame(header, samp_attrib).T
        df_lane_attrib = pandas.DataFrame(lane_attrib[1:len(lane_attrib)], columns=lane_attrib[0])

    return(dfcounts, df_samp_attrib, df_lane_attrib)

def parse_Ab_ref(abfile):
    """
    get antibody shortnames.
    Args:
        ab_handle
    Returns:
        dict
    Examples:
    """
    ab_name ={}
    for line in abfile:
        if not line.lstrip().startswith('#'):
            items = line.strip().split(",")
            ab_name[items[1]] = items[0]
    return(ab_name)

def parse_samplesheet(samplesheet):
    """
    get sample names, expecting to get more complex with this.
    ie. make directory structure to store data, etc
    Args:
        samplesheet
    Returns:
        dict
    Examples:
    """
    ss_dict = {}
    with samplesheet:
        count = 0
        for line in samplesheet:
            items = line.strip().split("\t")
            ss_dict[count] = items[1]
            count += 1
    return(ss_dict)


# norm_geomean functions

#geometric mean
def gm_mean(vector):
    """Get the geometric mean of a vector of numbers
    Args:
        vector
    Returns:
        float
    """
    log_vector = []
    for num in vector:
        log_vector.append(math.log(float(num)))
    num_sum = math.fsum(log_vector)
    g = math.exp(num_sum/len(vector))
    return g

#get output
def get_output(input_file):
    """Format the raw data input file into usable arrays of Samples, positive and negative
    Args:
        input_file: the raw data file
    Returns:
        dict with keys: pos, neg and samples. Values are two dimensional arrays.
    """
    #read in the file. separator is a \t 
    with open(input_file, 'r') as raw_file:
        reader = csv.reader(raw_file, delimiter="\t")
        file_array = []
        for row in reader:
            file_array.append(row)
    #get the name column, by checking for Name in header
    count = 0
    name_col = 1
    for header in file_array[0]:
        if header.strip() == "Name":
            name_col = count
        else:
            count += 1
    #get the first sample column by checking for an integer entry
    sample_col = 3
    count = 0
    for entry in file_array[1]:
        try:
            num = int(entry)
            sample_col = count
            break
        except ValueError:
            count += 1
    
    #iterate over the array, adding each line to the appropriate new array.
    samples = []
    pos = []
    neg = []
    header = True
    for row in file_array:
        if header:
            pos.append([row[name_col]] + row[sample_col:])
            samples.append([row[name_col]] + row[sample_col:])
            neg.append([row[name_col]] + row[sample_col:])
            header = False
        elif "POS" in row[name_col]:
            pos.append([row[name_col]] + row[sample_col:])
        elif "NEG" in row[name_col]:
            neg.append([row[name_col]] + row[sample_col:])
        else:
            samples.append([row[name_col]] + row[sample_col:])
            
    #samples are a matrix with rows are genes and columns are samples
    output_dict = {"pos":pos, "neg":neg, "samples":samples}
    return output_dict
    
            
#geomean norm
def geomean_norm(samples, controls, pos, mm, rb):
    """Normalizes the data
    Args:
        samples: two dimensional samples matrix from get_output
        controls: list of control samples
        pos: two dimensinal pos matrix from get_output
        mm: list of mouse antibodies
        rb: list of rabbit antibodies
    returns:
        list: list of four two dimensional arrays from after each of the three normalization steps and the log step.
    """
    #samples are matrix of sampls from get_output
    #controls are list of control names (columns) based on the ctrl_flag
    #pos is matrix of positive controls, formatted like the sample matrix
    #mm is a list of mouse antibodies
    #rb is a list of rabbit antibodies
    
    num_controls = 0.0
    sum_of_sums = 0.0
    header = True
    for row in pos:
        if not header:
            count = 0
            num_controls = 0.0
            for col in row:
                if pos[0][count].strip() in controls:
                    sum_of_sums += float(col)
                    num_controls += 1
                count += 1
        header = False
    mean_of_sums = sum_of_sums/num_controls
    
    for sample_index in range(1, len(pos[0])):
        sample_sum = 0
        for row_index in range(1, len(pos)):
            sample_sum += float(pos[row_index][sample_index])
        cf = mean_of_sums/sample_sum
        for row_index in range(1, len(samples)):
            samples[row_index][sample_index] = float(samples[row_index][sample_index]) * cf
            
    samples_1 = copy.deepcopy(samples)
    
    geomeans = []    
    for sample_index in range(1, len(samples[0])):
        if samples[0][sample_index].strip() in controls:
            sample_col = []
            for row in samples:
                sample_col.append(row[sample_index])
            sample_col.pop(0)
            sample_geomean = gm_mean(sample_col)
            geomeans.append(sample_geomean)
            
    grand_mean = sum(geomeans)/len(geomeans)
            
    for sample_index in range(1, len(samples[0])):
        sample_col = []
        for row in samples:
            sample_col.append(row[sample_index])
        sample_col.pop(0)    
        sample_geomean = gm_mean(sample_col)
        cf = grand_mean/sample_geomean
        header = True
        for row in samples:
            if header:
                header = False
            else:
                row[sample_index] = row[sample_index] * cf
    
    samples_2 = copy.deepcopy(samples)
    
    mouse_igg_index = 0
    rabbit_igg_index = 0
    count = 0
    for row in samples:
        if row[0].strip() == "MmAb-IgG1":
            mouse_igg_index = count
        elif row[0].strip() == "RbAb-IgG":
            rabbit_igg_index = count
        count += 1
        
        
    for sample_index in range(1, len(samples[0])):
        mouse_igg = samples[mouse_igg_index][sample_index]
        rabbit_igg = samples[rabbit_igg_index][sample_index]
        
        for row in samples:
            if row[0].strip() in mm:
                row[sample_index] = row[sample_index]/mouse_igg
            elif row[0].strip() in rb:
                row[sample_index] = row[sample_index]/rabbit_igg
                
    samples_3 = copy.deepcopy(samples)
    
    for row in samples:
        for index in range(1, len(samples[0])):
            try:
                row[index] = math.log(row[index], 2)
            except TypeError:
                pass
                
                
    output_list = [samples_1, samples_2, samples_3, samples]            
    return output_list

def run_norm_geomean(input_file, ctrl):
    """Runs the entire norm geomean process from the raw data file to writing the normalized output files.
    Args:
        input_file: the raw data file output from the processing steps
        ctrl: Boolean indicating whether all samples should be used as controls or just the default ones.
    Returns:
        none
    """
    #process the input file
    raw_data = get_output(input_file)
    #set controls
    global controls
    if ctrl:
        #set controls to all the samples
        controls = raw_data["samples"][0][1:]
        
    #set the mm to all the rows that are from the mouseAb list
    #set the rb to all the other rows
    mm = []
    rb = []
    for row in raw_data["samples"]:
        if row[0] in mouseAb:
            mm.append(row[0])
        elif row[0] != "Name":
            rb.append(row[0])
    
    #get the normalized data
    norm_dat = geomean_norm(raw_data["samples"], controls, raw_data["pos"], mm, rb)
    
    #get the prefix for the name for the output file
    outfile_name = input_file.split('_rawdata')[0]
    
    #write the normalized data to a csv file with the name ending in _NORMALIZED.tsv
        
    normfile_name = "1_ERCC_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer_1 = csv.writer(csvfile, delimiter="\t")
        norm_writer_1.writerows(norm_dat[0])
    
    normfile_name = "2_GEOMEAN_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer_2 = csv.writer(csvfile, delimiter="\t")
        norm_writer_2.writerows(norm_dat[1])
        
    normfile_name = "3_IGG_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer_2 = csv.writer(csvfile, delimiter="\t")
        norm_writer_2.writerows(norm_dat[2])
        
    normfile_name = "4_LOG_2_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer = csv.writer(csvfile, delimiter="\t")
        norm_writer.writerows(norm_dat[3])



def main():

    args = supply_args()
    if len(args.rcc_files) < 12:
        lrcc = len(args.rcc_files)
        print("Only ", lrcc, " RCC files specified!\nContinuing with ", lrcc, " files.")

    sampleid = parse_samplesheet(args.samplesheet)
    print(sampleid)

    rcc_counts_dict = {}
    samp_attrib_dict = {}
    lane_attrib_dict = {}

    outputDir = os.path.dirname(os.path.abspath(args.rcc_files[1]))

    for file in args.rcc_files:
        if file.endswith(".RCC"):
            #get sample number
            #samp_number = re.sub(".RCC", "", file.split("_")[-1])
            dfrcc, df_samp_attrib, df_lane_attrib = parseRCC(file)
            samp_number = int(df_lane_attrib.columns[1])
            #rename column name with sample id
            dfrcc.rename(columns={'Count': sampleid[samp_number]}, inplace=True)
            rcc_counts_dict[sampleid[samp_number]] = dfrcc

            #df_lane_attrib.rename(columns={df_lane_attrib.columns[1]: sampleid[file]}, inplace=True)
            df_lane_attrib = df_lane_attrib.append(pandas.Series(['SampleName', sampleid[samp_number]],
                                            index=df_lane_attrib.columns),ignore_index=True)
            lane_attrib_dict[df_lane_attrib.columns[1]] = df_lane_attrib

            samp_attrib_dict[sampleid[samp_number]] = df_samp_attrib

            #batch_name=re.sub(r'_%s+.RCC' % samp_number, "", file)

    #merge dictionary of dataframes together
    raw_data = reduce(lambda x, y: pandas.merge(x, y, on=['CodeClass', 'Name', 'Accession']),
                      rcc_counts_dict.values())

    lane_attrib = reduce(lambda x, y: pandas.merge(x, y, on=['ID']),
                         lane_attrib_dict.values())



    # Change long name to something else if necessary
    if args.abfile:
        with args.abfile:
            ab_name = parse_Ab_ref(args.abfile)
        for name in raw_data["Name"].iteritems():
            if not re.search("POS|NEG", name[1]):
                #print(name)
                raw_data['Name'][name[0]] = ab_name[name[1].strip().split("|")[0]]

    raw_data.to_csv(outputDir + "/rawdata.txt", sep='\t', index=False)

    #test to see if samp_attribs are all the same
    last_value = None
    first = True
    for value in samp_attrib_dict.values():
        if first:
            last_value = value
            first = False
            with open(outputDir + "/run_metrics.txt", 'w') as met:
                met.write(value.to_csv(sep="\t"))
                met.write(lane_attrib.to_csv(sep="\t"))
        #if samp_attrib_dict.values()[v].equals(samp_attrib_dict.values()[v+1]) != True:
        if value.equals(last_value) != True:
            print("Samples are not from one batch. Sample Attributes differ")
            print(samp_attrib_dict.keys()[v])
        else:
            print("Samples are from one batch. OK")
        last_value = value
    
        

    #run the norm_geomean
    raw_data_file = outputDir + "/rawdata.txt"
    run_norm_geomean(raw_data_file, args.ctrl)

if __name__ == "__main__":
    main()