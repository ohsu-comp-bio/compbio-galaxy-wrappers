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
import re
import argparse
import pandas
import xml.etree.ElementTree as ET
import math
import csv
import copy
from functools import reduce

VERSION="1.0.0"

def define_controls():
    """
        Made this function to avoid incorporating global variables
        And this may be more complex in the future.
    Args: no args
    Return: list of controls, and control help string
    """
    controls = ["MCF7","HCC1954","BT474","HeyA8","MDA468 control","MDA468control", "MDA468+EGF"]

    #make the control help statement for norm geomean ctrl input
    control_help = "Marks all samples to be used as control. Otherwise controls are "+ ",".join([str(i) for i in controls])
    return(controls, control_help)

def supply_args():
    """
    Input arguments
    """
    controls, control_help=define_controls()
    parser = argparse.ArgumentParser(description='Nanostring output RCC readouts (html format) and '
                                                 'converts to raw data tsv')
    parser.add_argument('rcc_files', type=str, nargs='+', help='raw RCC files')
    parser.add_argument('--samplesheet', type=str, help='samplesheet.txt')
    parser.add_argument('--abfile', type=argparse.FileType('r'),help='ANTIBODY_REFERENCE.csv')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    #For norm_geomean
    parser.add_argument('-ctrl', action="store_true", default=False, help=control_help)
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
        df_lane_attrib = pandas.DataFrame(lane_attrib[1:len(lane_attrib)], columns=lane_attrib[0])
        #df_samp_attrib = pandas.DataFrame(header, samp_attrib).T
        #df_samp_attrib = pandas.DataFrame(header[0:1], samp_attrib).T
    return(dfcounts, header, samp_attrib, df_lane_attrib)

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
    mouseAb = []
    rabbitAb =[]
    for line in abfile:
        if not line.lstrip().startswith('#'):
            items = line.strip().split(",")
            ab_name[items[1]] = items[0]
            if items[4] == "mouse":
                mouseAb.append(items[0])
            else:
                rabbitAb.append(items[0])
    return(ab_name, rabbitAb, mouseAb)

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
    with open(samplesheet, 'r') as ss:
        count = 0
        for line in ss:
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

def df2dictarrays(raw_data):
    """Format the raw data input file into usable arrays of Samples, positive and negative
    Args:
        rawdata pandas dataframe
    Returns:
        dict with keys: pos, neg and samples. Values are two dimensional arrays.
    """
    # iterate over the dataframe, adding each line to the appropriate new array.
    samples = []
    pos = []
    neg = []
    pos.append(list(raw_data.columns))
    samples.append(list(raw_data.columns))
    neg.append(list(raw_data.columns))
    for index,row in raw_data.iterrows():
        if "POS" in row.Name:
            pos.append(row.to_list())
        elif "NEG" in row.Name:
            neg.append(row.to_list())
        else:
            samples.append(row.to_list())

    # samples are a matrix with rows are antibodies and columns are samples
    output_dict = {"pos": pos, "neg": neg, "samples": samples}
    return output_dict
            
#geomean norm
def geomean_norm(samples, controls, pos, mouseAb, rabbitAb):
    """Normalizes the data
    Args:
        samples: two dimensional samples matrix from get_output
        controls: list of control samples
        pos: two dimensinal pos matrix from get_output
        mouseAb: list of mouse antibodies
        rabbitAb: list of rabbit antibodies
    returns:
        list: list of four two dimensional arrays from after each of the three normalization steps and the log step.
    """
    #samples are matrix of sampls from get_output
    #controls are list of control names (columns) based on the ctrl_flag
    #pos is matrix of positive controls, formatted like the sample matrix
    #mouseAb is a list of mouse antibodies
    #rabbitAb is a list of rabbit antibodies
    
    #num_controls = 0.0
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
    #    for row in samples:
    #        if row[0].strip() in mouseAb:
    #            row[sample_index] = row[sample_index] / mouse_igg
    #        elif row[0].strip() in rabbitAb:
    #            row[sample_index] = row[sample_index] / rabbit_igg
        for row in samples:
            if row[0].strip() in mouseAb:
                if row[sample_index]-mouse_igg <= 0:
                    row[sample_index] = 0
                else:
                    row[sample_index] = row[sample_index] - mouse_igg
            elif row[0].strip() in rabbitAb:
                if row[sample_index]-rabbit_igg <= 0:
                    row[sample_index] = 0
                else:
                    row[sample_index] = row[sample_index] - rabbit_igg

    samples_3 = copy.deepcopy(samples)

    for row in samples:
        for index in range(1, len(samples[0])):
            try:
                row[index] = math.log(row[index]+1, 2)
            except TypeError:
                pass
                
                
    output_list = [samples_1, samples_2, samples_3, samples]            
    return output_list

def write_norms(norm_dat):
    """Runs the entire norm geomean process from the raw data file to writing the normalized output files.
    Args:
        #norm_dat: list of 4 arrays
    Returns:
        none
    """

    #write the normalized data to a csv file with the name ending in _NORMALIZED.tsv
        
    normfile_name = "1_ERCC_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer_1 = csv.writer(csvfile, delimiter="\t")
        norm_writer_1.writerows(norm_dat[0])
    
    normfile_name = "2_GEOMEAN_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer_2 = csv.writer(csvfile, delimiter="\t")
        norm_writer_2.writerows(norm_dat[1])

    normfile_name = "3_IGG_SUBTRACTED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer_2 = csv.writer(csvfile, delimiter="\t")
        norm_writer_2.writerows(norm_dat[2])
        
    normfile_name = "4_LOG_2_NORMALIZED.tsv"
    with open(normfile_name, 'w') as csvfile:
        norm_writer = csv.writer(csvfile, delimiter="\t")
        norm_writer.writerows(norm_dat[3])



def main():
    #avoiding global variables for controls and this might get more complicated in the future
    controls, control_help = define_controls()

    args = supply_args()
    if len(args.rcc_files) < 12:
        lrcc = len(args.rcc_files)
        print("\nOnly "+ str(lrcc) + " RCC files specified!\nContinuing with "+ str(lrcc) + " files.\n")

    sampleid = parse_samplesheet(args.samplesheet)
    print(sampleid)

    rcc_counts_dict = {}
    samp_attrib_dict = {}
    header_dict={}
    lane_attrib_dict = {}

    #outputDir = os.path.dirname(os.path.abspath(args.rcc_files[1]))
    for file in args.rcc_files:
        if file.endswith(".RCC"):
            #get sample number
            #samp_number = re.sub(".RCC", "", file.split("_")[-1])
            dfrcc, header, samp_attrib, df_lane_attrib = parseRCC(file)

            samp_number = int(df_lane_attrib.columns[1])
            #rename column name with sample id for rcc files
            dfrcc.rename(columns={'Count': sampleid[samp_number]}, inplace=True)
            rcc_counts_dict[sampleid[samp_number]] = dfrcc

            #df_lane_attrib.rename(columns={df_lane_attrib.columns[1]: sampleid[file]}, inplace=True)
            df_lane_attrib = df_lane_attrib.append(pandas.Series(['SampleName', sampleid[samp_number]],
                                            index=df_lane_attrib.columns),ignore_index=True)
            lane_attrib_dict[df_lane_attrib.columns[1]] = df_lane_attrib

            samp_attrib_dict[sampleid[samp_number]] = samp_attrib
            header_dict[sampleid[samp_number]] = header

    #Check header and samp_attribs

    for idx in range(1, len(header_dict)):
        if header_dict.values()[idx-1]!=header_dict.values()[idx]:
            print("RCC header are not equal")
            #print(header_dict.values()[idx])

    # test to see if samp_attribs are all the same, all RCC files are from same run
    for idx in range(1, len(samp_attrib_dict)):
        if samp_attrib_dict.values()[idx-1]!=samp_attrib_dict.values()[idx]:
            print("RCC Sample Attributes are not equal")
            #print(samp_attrib_dict.values()[idx])
            if samp_attrib_dict.values()[idx][0][1:]==samp_attrib_dict.values()[idx][0][1:]:
                print("Actually only ID don't match, is this an old RCC file or have chaned RLFs")
            else:
                print("RCC Sample Attributes IDs and values are not equal. Stopping")
                print("Samples are not from one batch. Sample Attributes differ")
                exit("Error, Check your RCC files")
        else:
            print("Samples attribs are okay and are from one batch. OK")





    #merge dictionary of dataframes together
    raw_data = reduce(lambda x, y: pandas.merge(x, y, on=['CodeClass', 'Name', 'Accession']),
                      rcc_counts_dict.values())

    lane_attrib_combined = reduce(lambda x, y: pandas.merge(x, y, on=['ID']),
                         lane_attrib_dict.values())


    with open("run_metrics.txt", 'w') as met:
        for item in header_dict.values()[0]:
            met.write('\t'.join(item))
            met.write('\n')
        for item in samp_attrib_dict.values()[0]:
            met.write('\t'.join(item))
            met.write('\n')
        met.write(lane_attrib_combined.to_csv(sep="\t"))


    # Change long name to something else if necessary
    with args.abfile:
        ab_name, rabbitAb, mouseAb = parse_Ab_ref(args.abfile)

    # Write raw data file
    for name in raw_data["Name"].iteritems():
        if not re.search("POS|NEG", name[1]):
            #print(name)
            raw_data['Name'][name[0]] = ab_name[name[1].strip().split("|")[0]]

    raw_data.to_csv("rawdata.txt", sep='\t', index=False)


    #Prep for dictionary of arrays
    raw_data_short = raw_data.drop(['CodeClass','Accession'], 1)
    split_data = df2dictarrays(raw_data_short)

    # set controls
    if args.ctrl:
        # set controls to all the samples
        controls = split_data["samples"][0][1:]
    else:
        controls

    norm_dat = geomean_norm(split_data["samples"], controls, split_data["pos"], mouseAb, rabbitAb)
    write_norms(norm_dat)

if __name__ == "__main__":
    main()
