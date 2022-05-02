#!/usr/bin/env python

VERSION = '0.0.1'


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
        with open(self.doc_file, 'r') as infile:
            i = 0
            for line in infile:
                i += 1
                if i > 1:
                    lsplit = line.split('\t')
                    if lsplit[0] not in doc_dict:
                        base_counts = lsplit[4]
                        base_dict = {'depth': lsplit[1],
                                     'A': base_counts.split(' ')[0].split(':')[1],
                                     'C': base_counts.split(' ')[1].split(':')[1],
                                     'G': base_counts.split(' ')[2].split(':')[1],
                                     'T': base_counts.split(' ')[3].split(':')[1],
                                     'N': base_counts.split(' ')[4].split(':')[1]}
                        doc_dict[lsplit[0]] = base_dict
        return doc_dict
