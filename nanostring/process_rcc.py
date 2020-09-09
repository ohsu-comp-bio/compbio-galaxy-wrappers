#!/usr/bin/env python

## USAGE:
#   Nanostring output RCC readouts for processing whole slide assay nCounter data
#   (html format) to raw data matrix
## NOTES:
#   TODO: break out tests
#
## ARGUMENTS:
#   samplesheet: tsv of RCC filename in one col and sample name, 'ANTIBODY_REFERENCE.csv'
#   rcc files
#   abfile: antibody file to rename antibody names.
## RETURN:
#   rawdata.txt: tab-sep file with table of Antibody x Sample
#   run_metrics.txt
#
import re
import argparse
import pandas
import xml.etree.ElementTree as ET
from functools import reduce

VERSION = "1.0.0.0"
def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='Nanostring output RCC readouts (html format) and '
                                                 'converts to raw data tsv')
    parser.add_argument('rcc_files', type=str, nargs='+', help='raw RCC files')
    parser.add_argument('--samplesheet', type=str, help='samplesheet.txt')
    parser.add_argument('--abfile', type=argparse.FileType('r'), help='ANTIBODY_REFERENCE.csv')
    parser.add_argument('--omit_blank', action='store_true', help='Omitt Blank Samples')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
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
        # df_samp_attrib = pandas.DataFrame(header, samp_attrib).T
        # df_samp_attrib = pandas.DataFrame(header[0:1], samp_attrib).T
    return (dfcounts, header, samp_attrib, df_lane_attrib)

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

def make_unique(seq): # Order preserving
  '''
  Makes values unique if they are duplicated
  '''
  seen = set(seq)
  if len(seq)==len(seen):
      print("all values unique")
      return seq
  uniq_seq = []
  for x in seq:
      if x in uniq_seq:
          i = 1
          newx = x + "-" + str(i)
          while newx in uniq_seq:
              i = i + 1
              newx = x + "-" + str(i)
          uniq_seq.append(x + "-" + str(i))
      else:
          uniq_seq.append(x)
  return uniq_seq

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
        for line in ss:
            items = line.strip().split("\t")
            rcc_int = int(items[0].split(".")[0][-2:])
            ss_dict[rcc_int] = items[1]
    #Fix duplicated sample names by making them unique
    if len(ss_dict.values())==len(set(ss_dict.values())):
        return (ss_dict)
    else:
        newval = make_unique(ss_dict.values())
        ss_keys = ss_dict.keys()
        ss_uniq = dict(zip(ss_keys, newval))
        return (ss_uniq)

def main():
    args = supply_args()
    if len(args.rcc_files) < 12:
        lrcc = len(args.rcc_files)
        print("\nOnly " + str(lrcc) + " RCC files specified!\nContinuing with " + str(
            lrcc) + " files.\n")

    sampleid = parse_samplesheet(args.samplesheet)
    print(sampleid)

    rcc_counts_dict = {}
    samp_metrics_dict = {}
    for file in args.rcc_files:
        #if file.endswith(".RCC"): #only works for non-galaxy instance
            # get sample number
            # samp_number = re.sub(".RCC", "", file.split("_")[-1])
        dfrcc, header, samp_attrib, df_lane_attrib = parseRCC(file)

        samp_number = int(df_lane_attrib.columns[1])
        # rename column name with sample id for rcc files
        dfrcc.rename(columns={'Count': sampleid[samp_number]}, inplace=True)
        rcc_counts_dict[sampleid[samp_number]] = dfrcc

        #sample_metrics add sample name
        df_lane_attrib = df_lane_attrib.append(
            pandas.Series(['SampleName', sampleid[samp_number]],
                          index=df_lane_attrib.columns), ignore_index=True)
        samp_attrib_df = pandas.DataFrame(samp_attrib + header, columns=df_lane_attrib.columns)

        df_samp_metrics = pandas.concat([samp_attrib_df, df_lane_attrib])

        df_samp_metrics = df_samp_metrics.set_index('ID')
        samp_metrics_dict[df_samp_metrics.columns[0]] = df_samp_metrics

    # merge dictionary of dataframes together
    raw_data = reduce(lambda x, y: pandas.merge(x, y, on=['CodeClass', 'Name', 'Accession']),
                      rcc_counts_dict.values())

    samp_metrics_df = pandas.concat(samp_metrics_dict.values(), axis=1)


    # Change long name to something else if necessary
    with args.abfile:
        ab_name, rabbitAb, mouseAb = parse_Ab_ref(args.abfile)

    # Write raw data file
    for name in raw_data["Name"].iteritems():
        if not re.search("POS|NEG", name[1]):
            # print(name)
            raw_data['Name'][name[0]] = ab_name[name[1].strip().split("|")[0]]
    print(raw_data.columns)

    #omit blanks if flag is true
    idx_of_samples = []
    if args.omit_blank == True:
        for label, content in raw_data.items():
            lastvalue = content.values[-1]
            try:
                int(lastvalue)
                colsum = content.astype('int').sum()
                if colsum > 0:
                    idx_of_samples.append(label)
            except ValueError:
                pass
        raw_data = raw_data[['CodeClass', 'Name', 'Accession'] + idx_of_samples]
        smet_boolidx = samp_metrics_df.loc['SampleName'].isin(idx_of_samples)
        samp_metrics_df = samp_metrics_df[smet_boolidx.index[smet_boolidx]]

    samp_metrics_df.to_csv("run_metrics.txt", sep='\t', index=True)
    raw_data.to_csv("rawdata.txt", sep='\t', index=False)

if __name__ == "__main__":
    main()
