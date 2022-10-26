#!/usr/bin/env python

VERSION = '0.0.1'


def calc_bkgd_est(base_depth, total_depth):
    try:
        bkgd = int(base_depth) / int(total_depth)
    except ZeroDivisionError:
        bkgd = 0.0
    return float(bkgd)


class DepthOfCoverageReader:
    """
    Ingest a GATK Depth of Coverage tsv file (with base counts)

    Locus	Total_Depth	Average_Depth_sample	Depth_for_RDM-PON10-1	RDM-PON10-1_base_counts
    3:38182641	8455	8455	8455	A:24 C:14 G:2 T:8415 N:0
    9:5073770	3327	3327	3327	A:8 C:1 G:3317 T:1 N:0
    """
    def __init__(self, doc_file):
        self.doc_file = doc_file
        self.doc = self.read()

    def read(self):
        doc_dict = {}
        with open(self.doc_file, 'r') as f:
            header_line = next(f)
            for line in f:
                lsplit = line.strip().split('\t')
                if lsplit[0] not in doc_dict:
                    base_counts_col = lsplit[4]
                    base_depth_list = []
                    for val in base_counts_col.split(' '):
                        base_depth_list.append(int(val.split(':')[1]))
                    base_dict = {'total_depth': int(lsplit[1]),
                                 'off_target': int(lsplit[1]) - max(base_depth_list),
                                 'A': base_depth_list[0],
                                 'C': base_depth_list[1],
                                 'G': base_depth_list[2],
                                 'T': base_depth_list[3],
                                 'N': base_depth_list[4]}
                    doc_dict[lsplit[0]] = base_dict
        return doc_dict
