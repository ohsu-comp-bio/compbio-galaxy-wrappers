#!/usr/bin/env python

import argparse

VERSION = '1.3.0'


def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('input_file', help='input file')
    #parser.add_argument('--cut_off', default=50)
    parser.add_argument('--flag_bp', default=9000000)
    parser.add_argument('--flag_len', default=40)
    parser.add_argument('--list_bp', default=3000000)
    parser.add_argument('--list_len', default=25)
    parser.add_argument('--list_num', default=15)
    args = parser.parse_args()
    return args


def is_hom(line):
    line = line.split("\t")
    genotype = line[2].strip()
    for nucl in ["A", "T", "C", "G"]:
        if genotype == nucl + "/" + nucl:
            return True
    return False


def get_runs(input_file, flag_bp, flag_len, list_bp, list_len):
    with open(input_file, 'r') as geno_file:
        chrom = "CHROM"
        start = 0
        count = 0
        runs = []
        flags = []
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
                if count >= list_len and not chrom == "X":
                    # and not chrom == "X"
                    end = line.split("\t")[1]
                    bp = int(end) - int(start)
                    #print(bp)
                    if bp >= list_bp:
                        runs.append((chrom, count, start, end, bp))
                    if count >= flag_len and bp >= flag_bp:
                        flags.append((chrom, count, start, end, bp))
                count = 0
                start = line.split("\t")[1]
            
    return [runs, chrom_ratios, float(sum_hets), float(sum_homs), flags]

    
def main():
    args = supply_args()
    with open(args.input_file) as f:
        first_line = f.readline()
        sample_gt = first_line.split("\t")[2]
        sample_name = sample_gt.split(".")[0]
        print(sample_name)
    
    stuff = get_runs(args.input_file, int(args.flag_bp), int(args.flag_len), int(args.list_bp), int(args.list_len))
    if stuff[3] == 0:
        stuff[3] = 1
    ratio = stuff[2]/stuff[3]
    print("Overall Ratio: " + str(ratio))
    for chrom in stuff[1]:
        print("Chr" + str(chrom[0]) + " Ratio: " + str(chrom[1]))
    stuff[0].sort(key=lambda tup: tup[1])
    print("Homozygosity Runs Over " + str(args.list_len) + ": ")
    for run in stuff[0]:
        print(str(run[1]) + " at chr" + str(run[0]) + ":" + str(run[2]) + "-" + str(run[3]) + " of length: " + str(run[4]))
    print("Flagged Runs:")
    for run in stuff[4]:
        print(str(run[1]) + " at chr" + str(run[0]) + ":" + str(run[2]) + "-" + str(run[3]) + " of length: " + str(run[4]))
        
    flag_file = open("flag_file.json", "w")
    
    if len(stuff[0]) >= int(args.list_num):
        flag_file.write('{"homozygosity_flag": 1}')
    elif len(stuff[4]) > 0:
        flag_file.write('{"homozygosity_flag": 1}')
    else:
        flag_file.write('{"homozygosity_flag": 0}')
    
    flag_file.close()
    
    
if __name__ == "__main__":
    main()