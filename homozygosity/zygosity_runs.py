#!/usr/bin/env python

import argparse

VERSION = '1.1.0'


def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', help='input file')
    args = parser.parse_args()
    return args


def is_hom(line):
    line = line.split("\t")
    genotype = line[2].strip()
    for nucl in ["A", "T", "C", "G"]:
        if genotype == nucl + "/" + nucl:
            return True
    return False


def get_runs(input_file):
    with open(input_file, 'r') as geno_file:
        chrom = "CHROM"
        start = 0
        count = 0
        runs = []
        sum_homs = 0.0
        sum_hets = 0.0
        homs = 0.0
        hets = 0.0
        chrom_ratios = []
        for line in geno_file:
            if not chrom == line.split("\t")[0]:
                if not chrom == "CHROM":
                    if homs == 0:
                        homs = 1
                    chrom_ratio = float(hets)/float(homs)
                    chrom_ratios.append((chrom, chrom_ratio))
                homs = 0
                hets = 0
                count = 0
                start = 0
                chrom = line.split("\t")[0]
                
            if is_hom(line):
                homs = homs + 1
                sum_homs = sum_homs + 1
                count = count + 1
            else:
                hets = hets + 1
                sum_hets = sum_hets + 1
                if count >= 50:
                    # and not chrom == "X"
                    end = line.split("\t")[1]
                    runs.append((chrom, count, start, end))
                count = 0
                start = line.split("\t")[1]
            
    return [runs, chrom_ratios, float(sum_hets), float(sum_homs)]

    
def main():
    args = supply_args()
    stuff = get_runs(args.input_file)
    if stuff[3] == 0:
        stuff[3] = 1
    ratio = stuff[2]/stuff[3]
    print("Overall Ratio: " + str(ratio))
    for chrom in stuff[1]:
        print("Chr" + str(chrom[0]) + " Ratio: " + str(chrom[1]))
    stuff[0].sort(key=lambda tup: tup[1])
    print("Homozygosity Runs Over 50:")
    for run in stuff[0]:
        print(str(run[1]) + " at chr" + str(run[0]) + ":" + str(run[2]) + "-" + str(run[3]))
    
    
if __name__ == "__main__":
    main()
