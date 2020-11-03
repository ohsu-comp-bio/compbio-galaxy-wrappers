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

import argparse
import pandas
import xml.etree.ElementTree as ET
import json
import os

VERSION = "0.0.0.1"

def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='Nanostring output RCC readouts (html format) and '
                                                 'converts to raw data tsv')
    parser.add_argument('rcc_files', type=str, nargs='+', help='raw RCC files')
    parser.add_argument('--samplesheet', type=str, help='samplesheet.txt')
    parser.add_argument('--pkc_dir', type=str, help='Assay pkc file')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

def parse_pkc(pkc_file):
    with open(pkc_file, 'r') as xml_handle:
        pkc = json.load(xml_handle)
    targets = pkc['Targets']
    abs = [i['DisplayName'] for i in targets]
    panel_name = pkc['Name']

    pg_dict = {key: [] for key in abs}
    for item in pkc['ProbeGroups']:
        for ab in item['Targets']:
            if ab in item['Targets']:
                pg_dict[ab].append(item['Name'])

    return panel_name, targets, pg_dict


def parseRCC_DSP(rcc_file):
    """
    Take RCC file output and dictionary of locations
    Args:
        path to RCC
        dictpkc: dictionary of pkc locations
    Returns:
        dataframe
    Examples:
    """
    dict_rcc = {}
    with open(rcc_file, 'r') as xml_handle:
        xml_root = "<root>" + str(xml_handle.read()) + "</root>"
        root = ET.fromstring(xml_root)
        for child in root:
            dict_rcc[child.tag] = child.text


        header = [line.strip().split(',') for line in dict_rcc['Header'].strip().split('\n')]
        samp_attrib = [line.strip().split(',') for line in
                       dict_rcc['Sample_Attributes'].strip().split('\n')]
        lane_attrib = [line.strip().split(',') for line in
                       dict_rcc['Lane_Attributes'].strip().split('\n')]
        fullcounts = {}
        for line in dict_rcc['Code_Summary'].strip().split('\n'):
            fullcounts[line.strip().split(',')[1]] = line.strip().split(',')[3]

        df_lane_attrib = pandas.DataFrame(lane_attrib[1:len(lane_attrib)], columns=lane_attrib[0])
        # df_samp_attrib = pandas.DataFrame(header, samp_attrib).T
        # df_samp_attrib = pandas.DataFrame(header[0:1], samp_attrib).T
    return (fullcounts, header, samp_attrib, df_lane_attrib)


def main():
    args = supply_args()
    rcc_counts_dict = {}
    samp_metrics_dict = {}
    for file in args.rcc_files:
        if file.endswith(".RCC"):
            #combines all coordinates in RCC files into one dictionary file
            counts_dict, header, samp_attrib, df_lane_attrib = parseRCC_DSP(file)
            samp_number = int(df_lane_attrib.columns[1])

            samp_attrib_df = pandas.DataFrame(samp_attrib + header, columns=df_lane_attrib.columns)
            df_samp_metrics = pandas.concat([samp_attrib_df, df_lane_attrib])

            df_samp_metrics = df_samp_metrics.set_index('ID')
            samp_metrics_dict[df_samp_metrics.columns[0]] = df_samp_metrics

            rcc_counts_dict[samp_number] = counts_dict

    counts = {}
    md_counts = {}
    pkc_files = os.listdir(args.pkc_dir)
    # 1) iterate over rccs
    for rccidx in rcc_counts_dict.keys():
        for letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
            roiid = str(rccidx) + letter
            counts.update({roiid:{}}) # empty dictionary in dictionary
    # 2) iterate over pkc files
            for pkc_file in pkc_files:
                panel_name, targets, pg_dict = parse_pkc(os.path.join(pkc_dir, pkc_file))
    # 3) iterate over each antibody/probe in pkc
                for ab in pg_dict.keys():
    # 4) find target, complex nested dict easier to find and reassign, first occurence
                    target = next((target for target in targets if target['DisplayName']==ab),None)
                    md_counts[ab] = [panel_name, pg_dict[ab], target['CodeClass'], target['AnalyteType']]
                    counts[roiid].update({ab: rcc_counts_dict[rccidx][target['DSP_ID'][letter]]})
    md = pandas.DataFrame.from_dict(md_counts, columns=["PanelName", "ProbeGroup","CodeClass", "AnalyteType"],
                                orient='index')
    pcounts = pandas.DataFrame.from_dict(counts)
    pandas.DataFrame.to_csv(pcounts, "rawdata.txt", sep="\t")

if __name__ == "__main__":
    main()