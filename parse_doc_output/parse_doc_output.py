#!/usr/bin/env python

### This will take read depth data from GATK DepthOfCoverage and will calculate Q30%
### as well as depth of coverages at varying depths.  DepthOfCoverage should be run
### twice, once with the option mbq=30 and once with mbq=0.  This will allow
### an easy calculation of these metrics per interval.
### Inputs: 2 DepthOfCoverage depth outputs and an interval list in BED format
### Run match_intervals_genes.py to create interval list with associated gene names.
### Output: tsv with chrom, interval start, interval stop, gene, q30%,
###         and columns of binned coverage values.
### John Letaw 032915 

import argparse
import collections

### Set bins globally so we don't have to pass the list over and over again.

BINS = [200, 100, 50, 25, 10, 1, 0]

def binDepths(depth, offset):
    """
    Input: a depth value
    Output: index of which bin to add depth
    """
    for bin in BINS:
        if int(depth) >= bin:
            return BINS.index(bin)+offset

def importIntervals(handle):
    """
    Input: file handle of genomic intervals
    Output: [chrom, start, stop, gene, q0, q30, d200, d100, d50, d25, d10, d1, d0]
    """

    print("Reading intervals file.")
    interval_list = []
    gene_dict = collections.OrderedDict()
    with handle as intervals:
        for intval in intervals:

            intval = intval.rstrip('\n')
            chrom = intval.split('\t')[1]
            start = int(intval.split('\t')[2])
            stop = int(intval.split('\t')[3])
            gene = intval.split('\t')[0]

            interval_list.append([chrom, start, stop, gene, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            region = (stop - start + 1)
            if gene not in gene_dict:
                gene_dict[gene] = [chrom, start, stop, region, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            else:
                gene_dict[gene][3] += region
                gene_dict[gene][2] = stop

    return interval_list, gene_dict

def populateIntervals(int_list, gene_dict, coord_dict_0, coord_dict_30):
    for intval in int_list:
        for i in range(intval[1], intval[2]):
            curr_key = intval[0] + ':' + str(i)
            if curr_key in coord_dict_0:
                curr_val_0 = coord_dict_0[curr_key]
                intval[binDepths(curr_val_0, 6)] += 1
                intval[4] += curr_val_0
                gene_dict[intval[3]][binDepths(curr_val_0, 6)] += 1
                gene_dict[intval[3]][4] += curr_val_0
            else:
                intval[12] += 1
                gene_dict[intval[3]][12] += 1
            
            if curr_key in coord_dict_30:
                curr_val_30 = coord_dict_30[curr_key]
                intval[5] += curr_val_30
                gene_dict[intval[3]][5] += curr_val_30

    return int_list, gene_dict

def importDepth(handle):
    """
    Input: file handle of depths from DepthOfCoverage, index of where to store counts, interval list
    Output: updated interval list
    """

    coord_dict = {}

    print("Reading DepthOfCoverage per locus output.")
    with handle as depth:
        next(depth)  ### Skip header line.
        i = 0
        for coord in depth:
            
            coord = coord.rstrip('\n').split('\t')
            dcoord = int(coord[0].split(':')[1])
            ddepth = int(coord[1])

            coord_dict[coord[0]] = int(coord[1])

    return coord_dict


def sumDepths(depths):
    """
    Input: Arbitrary list of depths
    Output: sum of the depths
    """
    sum = 0
    for d in depths:
        sum += d
    return sum

def formatPercent(x, y):    
    """
    Input: Two float values.
    Output: Formatted percent string.
    """
    if y != 0.0:
        return str("{0:.1f}".format(x/y*100))
    return "0.0"
    
def writeIntervalQC(handle, int_list):
    """
    Input: file handle for writing results, interval list
    Output: written file
    """

    print("Writing Interval QC.")
    handle.write("Chromosome\tStart\tStop\tGene\tAvgD\tQ30\tD200\tD100\tD50\tD25\tD10\n")
    for interval in int_list:
        handle.write(interval[0] + '\t' + str(interval[1]) + '\t' + str(interval[2]) + '\t' + interval[3] + '\t')
        isize = float(interval[2] - interval[1] + 1)
        handle.write(str("{0:.1f}".format(interval[4]/isize)) + '\t')
        handle.write(formatPercent(interval[5], interval[4]))
        for i in range(7, len(interval)-1):
            handle.write('\t' + formatPercent(sumDepths(interval[6:i]), sumDepths(interval[6:])))
        handle.write('\n')
    handle.close()

def writeGeneQC(handle, gene_dict):
    
    print("Writing Gene QC.")
    handle.write("Chromosome\tStart\tStop\tGene\tAvgD\tQ30\tD200\tD100\tD50\tD25\tD10\n")
    for gene in gene_dict:
        handle.write(gene_dict[gene][0] + '\t' + str(gene_dict[gene][1]) + '\t' + str(gene_dict[gene][2]) + '\t' + gene + '\t')
        isize = float(gene_dict[gene][3])
        handle.write(str("{0:.1f}".format(gene_dict[gene][4]/isize)) + '\t')
        handle.write(formatPercent(gene_dict[gene][5], gene_dict[gene][4]))
        for i in range(7, len(gene_dict[gene])-1):
            handle.write('\t' + formatPercent(sumDepths(gene_dict[gene][6:i]), sumDepths(gene_dict[gene][6:])))
        handle.write('\n')
    handle.close()

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("interval_list", help="List of intervals in form chr:start-stop")
    parser.add_argument("depth_0", help="DepthOfCoverage output at mbq=0")
    parser.add_argument("depth_30", help="DepthOfCoverage output at mbq=30")
    parser.add_argument("output_int", help="Output Interval TSV")
    parser.add_argument("output_gene", help="Output Gene TSV")
    args = parser.parse_args()

    coord_dict_0 = importDepth(open(args.depth_0, 'rU'))
    coord_dict_30 = importDepth(open(args.depth_30, 'rU'))
    interval_list, gene_dict = importIntervals(open(args.interval_list, 'rU'))
    interval_list, gene_dict = populateIntervals(interval_list, gene_dict, coord_dict_0, coord_dict_30)
    
    writeIntervalQC(open(args.output_int, 'w'), interval_list)
    writeGeneQC(open(args.output_gene, 'w'), gene_dict)

if __name__ == "__main__":
    main()
