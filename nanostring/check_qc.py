#!/usr/bin/env python

## USAGE:
#
#
## ARGUMENTS:
#
## RETURN:
#
#

import math
import argparse
import pandas
import statistics
import re
import json

VERSION="0.0.1.0"
def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='check basic QC from RCC files')
    parser.add_argument('--rawdata', type=str, help='raw_data.txt')
    parser.add_argument('--runmetrics', type=str, help='run_metrics.txt')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


def checkFOV(fov_count, fov_counted):
    """
    Registered FOVs. FOV Counted (Fields of View successfully counted)/ FOV
    Count (Fields of View attempted) >=0.75
    :param fov_count: double
    :param fov_counted:double
    :return: dict
    """
    fovs = float(fov_count)/float(fov_counted)
    if fovs >= 0.75:
        qc = "PASS"
    else:
        qc = "FAIL"
    fov_dict = {"value": round(fovs,4),
                "threshold" : ">0.75",
                "qc": qc}
    return fov_dict

def check_binding_density(bd):
    """
    Binding density threshold 0.1-0.8 SPRINT, 0.1-2.25 for max
    :param bd: double
    :return: dict qc="PASS/FAIL"
    """
    bd=float(bd)
    if bd >= 0.1 and bd <= 2.25:
        qc = "PASS"
    else:
        qc = "FAIL"

    bd_dict = {"value": bd,
               "threshold": "0.1-2.25",
               "qc": qc}
    return bd_dict

def check_posE_limit_of_detection(list_neg, pos_e):
    """
    To determine whether POS_E is detectable, it can be compared to the counts for the negative
    control probes. Every nCounter Gene Expression assay is manufactured with eight negative control
    probes that should not hybridize to any targets within the sample. Counts from these will
    approximate general non-specific binding of probes within the samples being run. The counts of
    POS_E should be higher than two times the standard deviation above the mean of the negative control.
    :param list_neg: list_neg
    :return: dict of mean, sd, string qc="PASS/FAIL"
    """
    list_neg = [float(i) for i in list_neg]
    pos_e = float(pos_e)
    mean = sum(list_neg) / len(list_neg)
    res = statistics.pstdev(list_neg)
    m2sd = sum([mean, 2*res])
    if float(pos_e) > m2sd:
        qc = "PASS"
    else:
        qc = "FAIL"

    pose_ld_dict = {"value": pos_e,
                    "threshold": ' >' + str(round(mean/len(list_neg),1)) +'+' + str(round(2*res,1)),
                    "qc": qc}

    return pose_ld_dict

def check_pos_linearity(dat, samp):
    """
    Six synthetic DNA control targets are included with every nCounter Gene Expression assay.
    Their concentrations range linearly from 128 fM to 0.125 fM, and they are referred to as
    POS_A to POS_F, respectively. These positive controls are typically used to measure the
    efficiency of the hybridization reaction, and their step-wise concentrations also make
    them useful in checking the linearity performance of the assay. An R2 value is calculated
    from the regression between the known concentration of each of the positive controls and
    the resulting counts from them (this calculation is performed using log2-transformed values).
    :param dat: dataframe
    :return: dict
    """
    pos_name = [item for item, item in enumerate(list(dat.index)) if "POS" in item]
    posf_idx = [i for i, item in enumerate(pos_name) if "POS_F" not in item] #omit pos_f
    pos_name.pop(posf_idx[0])
    rexp = re.compile(r"((?<=\()(.*?)(?=\)))")
    conc = []
    for pos in pos_name:
        val=float(re.search(rexp, pos).group())
        val=math.log(sum([val,1]))
        conc.append(val)

    list_pos = dat.loc[pos_name, samp].tolist()
    list_pos = [math.log(sum([float(i),1])) for i in list_pos]
    #list_pos = [float(i) for i in list_pos]
    cor_df = pandas.DataFrame({'conc':conc, 'pos': list_pos})
    corval= (cor_df.corr().iloc[0]) * (cor_df.corr().iloc[0])

    if round(corval.loc['pos'],1) >= 0.90:
        qc = "PASS"
    else:
        qc = "FAIL"

    pos_linear_dict = {"value": round(corval.loc['pos'],1), "threshold": ">0.9","qc": qc}

    return pos_linear_dict

def main():
    args = supply_args()

    raw_data = pandas.read_csv(args.rawdata, sep="\t", index_col=False)
    run_metrics = pandas.read_csv(args.runmetrics, sep="\t", index_col=0)
    run_metrics.columns = run_metrics.loc["SampleName"].values.tolist()
    samps = raw_data.columns.drop(['CodeClass', 'Name', 'Accession'])

    qc = {}
    for i in range(0,len(samps)):
        qc_samp={}
        samp = samps[i]
        met = run_metrics[[samp]]
        dat = raw_data[['Name', samp]]
        dat.set_index('Name', inplace=True)

        bd_dict = check_binding_density(met._get_value("BindingDensity", samp))
        fov_dict = checkFOV(met._get_value("FovCount", samp), int(met._get_value("FovCounted", samp)))

        #limit of detection
        list_neg = dat.loc[[item for item, item in enumerate(list(dat.index)) if "NEG" in item],samp].tolist()
        pose = dat._get_value("POS_E(0.5)",samp)
        pose_ld_dict = check_posE_limit_of_detection(list_neg, pose)

        #ercc linearity
        pos_line_dict = check_pos_linearity(dat, samp)
        #make samp dict
        qc_samp["Binding Density"]=bd_dict
        qc_samp["Fov Ratio"] = fov_dict
        qc_samp["Limit Of Detection"] = pose_ld_dict
        qc_samp["POS ERCC R2"] = pos_line_dict
        qc[samp]=qc_samp

    with open('qc_metrics.json', 'w') as json_file:
        json.dump(qc, json_file, indent=4, sort_keys=True)


if __name__ == "__main__":
    main()
