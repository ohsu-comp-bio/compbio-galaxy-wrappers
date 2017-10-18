#!/usr/bin/env python

###
# Create the necessary coverage metrics for Richards and Corless lab variant calling workflows.
# This script will produce a gene, exon, and probe level coverage metrics sheet.  There are slight
# differences in the creation of gene level metrics between the TruSightOne gene panels workflow
# and the Agilent CRE exome workflow.  Check --help for usage.
### John Letaw 10/27/15

import argparse

DEPTHS = [200, 100, 50, 25, 10, 0]
#DEPTHS = [2000, 500, 100, 20, 10, 0]

def createDocDict(handle):

    """ Create dictionary to hold coverage values per coordinate.  Input file originates
    from GATK's DepthOfCoverage walker.
    DepthOfCoverage will be run on an interval list from the capture kit manufacturer.
    
    Input: handle: file handle
    Output: doc_dict: {"chrom:coord": depth} """

    print("Creating depth of coverage dictionary.")

    doc_dict = {}
    with handle as doc:
        i = 0
        next(doc) # Skip header line.
        for line in doc:
            if i % 1000000 == 0:
                print(i)
            line = line.rstrip('\n').split('\t')
            doc_dict[line[0]] = int(line[1])
            i += 1
 
    return doc_dict


def createExonDict(handle, doc_dict_0, doc_dict_30, exon_buff):
    
    """ Create dictionary from exon source files:
    
    Exons file looks like:
    1       11872   12227   OR4F5 ENSE00002234632
    OFS = '\t'
    
    Inputs: handle: file handle 
            doc_dict_0: doc_dict 
            doc_dict_30: doc_dict
    Output: exon_dict: {ENSE_HGNC = [HGNC, read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]} """

    print("Creating exon dictionary.")

    exon_dict = {}
    with handle as exons:
        i = 0
        next(exons)
        for line in exons:
            if i % 100000 == 0:
                print(i)
            line = line.rstrip('\n').split('\t')

            chrom = line[0]
            start = line[1]
            stop = line[2]
            hgnc = line[3]
            ense = line[4]

            read_count_0 = 0.0
            read_count_30 = 0.0
            this_depth = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]            

            for coord in range(int(start)-exon_buff+1, int(stop)+1+exon_buff):
                comp_coord = line[0] + ':' + str(coord)
                if comp_coord in doc_dict_0 and comp_coord in doc_dict_30:
                    read_count_0 += doc_dict_0[comp_coord]
                    read_count_30 += doc_dict_30[comp_coord]
                    this_depth = calcDepth(doc_dict_0[comp_coord], this_depth) ### Avg. depth metrics based on Q30.
                else:
                    pass ### Add logging here for if exception.
        

            total_range = int(stop) - int(start) + 1 + (2*exon_buff)

            uniq_key = ense + '_' + hgnc
            exon_dict[uniq_key] = [hgnc, read_count_0, read_count_30, total_range]
            exon_dict[uniq_key].extend(this_depth)
            exon_dict[uniq_key].extend([chrom, int(start), int(stop)])
            exon_dict[uniq_key].extend([chrom, start, stop])
            i += 1

    return exon_dict


def writeExonQC(handle, exon_dict):

    """
    Write the QC for exons.
    Input: file handle
           exon_dict (below)
    Output: tsv with QC metrics
    exon_dict now looks like:
    {ENSE_HGNC = [HGNC, read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]}
    """

    print("Writing Exon QC.")

    handle.write("Chromosome\tStart\tStop\tGene\tExon\tAvgD\tQ30\tD200\tD100\tD50\tD25\tD10\n")

    for key in exon_dict:

        read_count_0 = exon_dict[key][1]
        read_count_30 = exon_dict[key][2]
        total_range = exon_dict[key][3]
        hgnc = exon_dict[key][0]

        if exon_dict[key][1] != 0:
            q30 = percent(read_count_30, read_count_0)
        else:
            q30 = 0.0
        avgd = "%.1f" % (read_count_30/total_range)
        d200 = percent(exon_dict[key][4], total_range)
        d100 = percent(exon_dict[key][5], total_range)
        d50 = percent(exon_dict[key][6], total_range)
        d25 = percent(exon_dict[key][7], total_range)
        d10 = percent(exon_dict[key][8], total_range)

        chrom = exon_dict[key][10]
        start = exon_dict[key][11]
        stop = exon_dict[key][12]
        
        handle.write(chrom + '\t' + str(start) + '\t' + str(stop) + '\t' + hgnc + '\t' + key.split('_')[0] + '\t' + str(avgd) + '\t' + str(q30) + '\t' + str(d200) + '\t' + str(d100) + '\t' + str(d50) + '\t' + str(d25) + '\t' + str(d10) + '\n')


def writeIntervalQC(handle, int_dict):

    """
    Write the interval (probe) level QC metrics to TSV.
    Input: file handle
           int_dict (below)
    int_dict now looks like:
    {chrom:start-stop = [read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0]}
    """

    print("Writing Interval QC.")

    handle.write("Chromosome\tStart\tStop\tAvgD\tQ30\tD200\tD100\tD50\tD25\tD10\n")

    sample_total = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    for key in int_dict:
        
        read_count_0 = int_dict[key][0]
        read_count_30 = int_dict[key][1]
        total_range = int_dict[key][2]

        chrom = key.split(':')[0]
        start = key.split(':')[1].split('-')[0]
        stop = key.split(':')[1].split('-')[1].rstrip('\n')
        
        if int_dict[key][1] != 0:
            q30 = percent(read_count_30, read_count_0)
        else:
            q30 = 0.0

        avgd = "%.1f" % (read_count_30/total_range)
        d200 = percent(int_dict[key][3], total_range)
        d100 = percent(int_dict[key][4], total_range)
        d50 = percent(int_dict[key][5], total_range)
        d25 = percent(int_dict[key][6], total_range)
        d10 = percent(int_dict[key][7], total_range)

        handle.write(chrom + '\t' + start + '\t' + stop + '\t' + str(avgd) + '\t' + str(q30) + '\t' + str(d200) + '\t' + str(d100) + '\t' + str(d50) + '\t' + str(d25) + '\t' + str(d10) + '\n')

        sample_total = [x + y for x, y in zip(sample_total, int_dict[key])]
        
    if sample_total[0] != 0:
        q30 = percent(sample_total[1], sample_total[0])
    else:
        q30 = 0.0

    avgd = "%.1f" % (sample_total[1]/sample_total[2])
    d200 = percent(sample_total[3], sample_total[2])
    d100 = percent(sample_total[4], sample_total[2])
    d50 = percent(sample_total[5], sample_total[2])
    d25 = percent(sample_total[6], sample_total[2])
    d10 = percent(sample_total[7], sample_total[2])

    ### Write the final totals row, for sample level metrics.

    handle.write('TOTAL\t\t\t' + str(avgd) + '\t' + str(q30) + '\t' + str(d200) + '\t' + str(d100) + '\t' + str(d50) + '\t' + str(d25) + '\t' + str(d10) + '\n')


def percent(num1, num2):
    
    """
    Find the num1 over num2 and return as percentage.
    Input: num1, num2
    Output: formatted percentage
    """

    return "%.1f" % (num1*100/num2)


def calcDepth(depth, total_depths):

    """
    Calculate the depths at values as defined by DEPTHS globally.
    Input: depth, total_depths
    Output: updated total_depths
    """

    if len(total_depths) != len(DEPTHS):
        raise Exception("total_depths is not the right size, please check code.")

    for value in DEPTHS:
        if depth >= value:
            total_depths[DEPTHS.index(value)] += 1
    return total_depths


def createIntervalDict(handle, doc_dict_0, doc_dict_30):

    """
    Create a dictionary from a probe interval file.
    For instance <chrom:start-stop>
    Inputs: handle, doc_dict_0, doc_dict_30
    Output: int_dict (below)
    {chrom:start-stop = [read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0]}
    """

    print("Creating Probe dictionary.")
    
    int_dict = {}
    with handle as intervals:
        i = 0
        for line in intervals:
            if i % 1000000 == 0:
                print(i)
            nline = line.rstrip('\n').split('-')
            read_count_0 = 0.0
            read_count_30 = 0.0
            this_depth = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            start = int(nline[0].split(':')[1])
            stop = int(nline[1])
            chrom = (nline[0].split(':')[0])

            if len(nline) == 2:
                for coord in range(start, stop+1):
                    comp_coord = chrom + ':' + str(coord)
                    if comp_coord in doc_dict_0 and comp_coord in doc_dict_30:
                        read_count_0 += doc_dict_0[comp_coord]
                        read_count_30 += doc_dict_30[comp_coord]
                        this_depth = calcDepth(doc_dict_0[comp_coord], this_depth)  ### Depth metrics based at Q30.
                total_range = stop - start + 1
                int_dict[line.rstrip('\n')] = [read_count_0, read_count_30, total_range]
                int_dict[line.rstrip('\n')].extend(this_depth)
                i += 1

    return int_dict


def createManifestDict(handle, doc_dict_0, doc_dict_30):

    """
    Create a data structure holding the data from a manifest file, such as from TruSightOne.
    This file lists intervals associated with gene names.
    Inputs: file handle, doc_dict_0, doc_dict_30
    Output: mani_dict (below)
    {HGNC = [read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]}
    """

    print("Creating Manifest dictionary.")
    
    mani_dict = {}
    gene_coords = {}

    with handle as intervals:

        i = 0
        for line in intervals:
            if i % 1000000 == 0:
                print(i)

            nline = line.rstrip('\n').split('\t')
            read_count_0 = 0.0
            read_count_30 = 0.0
            total_range = 0.0
            this_depth = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            start = int(nline[0].split(':')[1].split('-')[0])
            stop = int(nline[0].split('-')[1])
            chrom = nline[0].split(':')[0]
            gene = nline[1]

            all_coord = nline[0]
            
            if gene not in gene_coords:
                gene_coords[gene] = []
            if gene not in mani_dict:
                mani_dict[gene] = [read_count_0, read_count_30, total_range]
                mani_dict[gene].extend(this_depth)
                mani_dict[gene].extend([chrom])
            for coord in range(start, stop+1):
                comp_coord = chrom + ':' + str(coord)
                if coord not in gene_coords[gene]:
                    gene_coords[gene].append(coord)
                    if comp_coord in doc_dict_0 and comp_coord in doc_dict_30:
                        read_count_0 += doc_dict_0[comp_coord]
                        read_count_30 += doc_dict_30[comp_coord]
                        this_depth = calcDepth(doc_dict_0[comp_coord], this_depth) ### Depth metrics at Q30.
            total_range = stop - start + 1
            mani_dict[gene][:2] = [x + y for x, y in zip(mani_dict[gene][:2], [read_count_0, read_count_30])]
            mani_dict[gene][3:9] = [x + y for x, y in zip(mani_dict[gene][3:9], this_depth)]
            i += 1
 
    for gene in gene_coords:
        mani_dict[gene][2] = len(gene_coords[gene])
        mani_dict[gene].extend([min(gene_coords[gene]), max(gene_coords[gene])])

    return mani_dict


def createGeneDict(ref_dict, coords, doc_dict_0, doc_dict_30):


    """
    Exons file looks like:
    1       11872   12227   OR4F5 ENSE00002234632
    OFS = '\t'
    """
    print("Creating gene dictionary.")

    gene_dict = {}

    i = 0
    for hgnc in coords:
        if i % 10000 == 0:
            print(i)

        if hgnc in ref_dict:
            chrom = ref_dict[hgnc][0]
            start = ref_dict[hgnc][1]
            stop = ref_dict[hgnc][2]


            if hgnc not in gene_dict:
                gene_dict[hgnc] = []
                this_depth = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]            
                read_count_0 = 0.0
                read_count_30 = 0.0

                total_range = 0
                for ival in coords[hgnc]:
                    for coord in range(ival[0], ival[1]):
                        comp_coord = chrom + ':' + str(coord)
                        if comp_coord in doc_dict_0 and comp_coord in doc_dict_30:
                            read_count_0 += doc_dict_0[comp_coord]
                            read_count_30 += doc_dict_30[comp_coord]
                            this_depth = calcDepth(doc_dict_0[comp_coord], this_depth)
                        total_range += 1
        
            gene_dict[hgnc] = [read_count_0, read_count_30, total_range]
            gene_dict[hgnc].extend(this_depth)
            gene_dict[hgnc].extend([chrom, start, stop])
        
        i += 1

    return gene_dict


def createRefGenes(handle):
    
    """
    Create reference gene dictionary, matching genes to coordinates.
    Input: file handle
    Ref Genes looks like:
    OR4F5	ENSG00000186092	ENST00000335137	NM_001005484	CCDS30547	1	69091	70008	+
    Output: ref_dict
    {GENE: [chrom, start, stop]}
    """

    ref_dict = {}

    with handle as genes:
        for gene in genes:
            gene = gene.rstrip('\n').split('\t')
            ref_dict[gene[0]] = [gene[5], gene[6], gene[7]]

    return ref_dict


def createCoordDict(exon_dict, exon_buff):

    """
    Coordinate dictionary, connect gene to coordinate.
    Input: exon_dict, exon_buff (basepairs extra padding for genomic coords)
    exon_dict now looks like:
    {ENSE_HGNC = [HGNC, read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]}
    Output: coord_dict
    {GENE: [(start, stop), (start, stop), ...]}

    """

    print("Creating coordinate dictinoary.")

    coord_dict = {}

    for exon in exon_dict:
        hgnc = exon_dict[exon][0]
        start = exon_dict[exon][11] - exon_buff
        stop = exon_dict[exon][12] + exon_buff

        if hgnc not in coord_dict:
            coord_dict[hgnc] = [(start,stop)]
        else:
            coord_dict[hgnc].append((start,stop))
     

    return coord_dict


def mergeCoordDict(coords):

    """
    Find unique coordinate regions, mainly for removing overlapping coordinates in QC generation.
    Input: coords
    {GENE: [(start, stop), (start, stop), ...]}
    Output: coords (updated)
    """

    print("Merging coordinate dictionary intervals.")
    
    for key in coords:
        new_ivals = []
        for ival in sorted(coords[key]):
            if new_ivals == []:
                new_ivals.append(ival)
                continue
            else:
                if ival[0] > new_ivals[-1][1]+1:
                    new_ivals.append(ival)
                elif ival[1] > new_ivals[-1][1]:
                    new_ivals[-1] = (new_ivals[-1][0], ival[1])

        coords[key] = new_ivals

    return coords


def writeGeneQC(handle, gene_dict):

    """
    Write the gene level QC coverage metrics.
    Input: file handle, gene_dict (below)
    gene_dict now looks like:
    {HGNC = [read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]}
    Output: TSV file of gene metrics
    """

    print("Writing gene QC.")

    handle.write("Chromosome\tStart\tStop\tGene\tAvgD\tQ30\tD200\tD100\tD50\tD25\tD10\n")

    for key in gene_dict:

        read_count_0 = gene_dict[key][0]
        read_count_30 = gene_dict[key][1]
        total_range = gene_dict[key][2]
        
        chrom = gene_dict[key][9]
        start = gene_dict[key][10]
        stop = gene_dict[key][11]

        if gene_dict[key][0] != 0:
            q30 = percent(read_count_30, read_count_0)
        else:
            q30 = 0.0
        avgd = "%.2f" % (read_count_30/total_range)
        d200 = percent(gene_dict[key][3], total_range)
        d100 = percent(gene_dict[key][4], total_range)
        d50 = percent(gene_dict[key][5], total_range)
        d25 = percent(gene_dict[key][6], total_range)
        d10 = percent(gene_dict[key][7], total_range)

        handle.write(str(chrom) + '\t' + str(start) + '\t' + str(stop) + '\t' + key + '\t' + str(avgd) + '\t' + str(q30) + '\t' + str(d200) + '\t' + str(d100) + '\t' + str(d50) + '\t' + str(d25) + '\t' + str(d10) + '\n')


# def createManifestGeneDict(int_dict):
    
#     """
#     Create gene level metrics based on a manifest of gene to probe relationships.
#     Should look like:
#     {HGNC = [read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]}

#     int_dict now looks like:
#     {chrom:start-stop = [read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, hgnc]}
#     """

#     print("Creating manifest gene dictionary.")

#     gene_dict = {}

#     for coord in int_dict:

#         chrom = coord.split(':')[0]
#         start = coord.split(':')[1].split('-')[0]
#         stop = coord.split('-')[1]

#         hgnc = int_dict[coord][9]
#         read_count_0 = int_dict[coord][0]
#         read_count_30 = int_dict[coord][1]
#         total_range = int_dict[coord][2]
#         d200 = int_dict[coord][3]
#         d100 = int_dict[coord][4]
#         d50 = int_dict[coord][5]
#         d25 = int_dict[coord][6]
#         d10 = int_dict[coord][7]
#         d0 = int_dict[coord][8]

#         if hgnc not in gene_dict:
#             gene_dict[hgnc] = int_dict[coord][:9]
#             gene_dict[hgnc].extend([chrom, start, stop])
#         else:

# #{HGNC = [read_count_0, read_count_30, total_range, D200, D100, D50, D25, D10, D0, chrom, start, stop]}

#             gene_dict[hgnc][0] += read_count_0
#             gene_dict[hgnc][1] += read_count_30
#             gene_dict[hgnc][2] += total_range
#             gene_dict[hgnc][3] += d200
#             gene_dict[hgnc][4] += d100
#             gene_dict[hgnc][5] += d50
#             gene_dict[hgnc][6] += d25
#             gene_dict[hgnc][7] += d10
#             gene_dict[hgnc][8] += d0
            
#             if start < gene_dict[hgnc][10]:
#                 gene_dict[hgnc][10] = start
#             if stop > gene_dict[hgnc][11]:
#                 gene_dict[hgnc][11] = stop


#     return gene_dict


def main():

    parser = argparse.ArgumentParser(description='Create Coverage metrics for Sequencing Data.')
    parser.add_argument("exon_tsv", help='Exon file with format [CHROM, START, STOP, HGNC, ENSE].')
    parser.add_argument("doc0", help='Output from GATK DepthOfCoverage when minimum quality is 0.')
    parser.add_argument("doc30", help='Output from GATK DepthOfCoverage when minimum quality is 30.')
    parser.add_argument("exonqc", help='Write exon level QC metrics to this file.')
    parser.add_argument("probe_ival", help='File containing probe interval regions in format <chrom:start-stop>')
    parser.add_argument("probeqc", help='Write probe level QC metrics to this file.')
    parser.add_argument("--gene_tsv", help='Gene file with format [GENE, ENSG, ENST, REFSEQ, CCDS, CHROM, START, STOP, STRAND].')
    parser.add_argument("geneqc", help='Write gene level QC metrics to this file.')
    parser.add_argument("qc_type", choices=['agilent', 'tso'], default='agilent', help="Metrics will be tailored to either 'agilent' or 'tso'")
    parser.add_argument("--mani_ival", help='File containing probe interval regions from manifest file first [column=<GENE.CHROM.START.STOP>]')
    parser.add_argument("--exon_buff", default=10, type=int, help='Set the number of bases to use as buffer around exon regions.')

    args = parser.parse_args()

    doc_input_0 = open(args.doc0, 'rU')
    doc_input_30 = open(args.doc30, 'rU')    
    doc_dict_0 = createDocDict(doc_input_0)
    doc_dict_30 = createDocDict(doc_input_30)

    all_exons = open(args.exon_tsv, 'rU')
    exon_dict = createExonDict(all_exons, doc_dict_0, doc_dict_30, args.exon_buff)
    exonqc = open(args.exonqc, 'w')
    writeExonQC(exonqc, exon_dict)
    exonqc.close()

    if args.qc_type == "agilent":

        if args.mani_ival or not args.gene_tsv:
            parser.error("Agilent CRE analysis requires --gene_tsv but does not utilize --mani_ival.")

        all_intervals = open(args.probe_ival, 'rU')
        interval_dict = createIntervalDict(all_intervals, doc_dict_0, doc_dict_30)
        intervalqc = open(args.probeqc, 'w')
        writeIntervalQC(intervalqc, interval_dict)
        intervalqc.close()
    
        coord_dict = createCoordDict(exon_dict, args.exon_buff)
        coord_dict = mergeCoordDict(coord_dict)
        ref_dict = createRefGenes(open(args.gene_tsv, 'rU'))
        gene_dict = createGeneDict(ref_dict, coord_dict, doc_dict_0, doc_dict_30)
        geneqc = open(args.geneqc, 'w')
        writeGeneQC(geneqc, gene_dict)
        geneqc.close()

    else:

        if not args.mani_ival or args.gene_tsv:
            parser.error("TruSightOne analysis requires --mani_ival but does not utilize or --gene_tsv.")

        all_intervals = open(args.probe_ival, 'rU')
        interval_dict = createIntervalDict(all_intervals, doc_dict_0, doc_dict_30)
        intervalqc = open(args.probeqc, 'w')
        writeIntervalQC(intervalqc, interval_dict)
        intervalqc.close()

        mani_intervals = open(args.mani_ival, 'rU')
        mani_dict = createManifestDict(mani_intervals, doc_dict_0, doc_dict_30)
        geneqc = open(args.geneqc, 'w')
        writeGeneQC(geneqc, mani_dict)
        geneqc.close()
        

if __name__ == "__main__":
    main()
