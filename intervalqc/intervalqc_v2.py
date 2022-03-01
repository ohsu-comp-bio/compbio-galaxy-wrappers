#!/usr/bin/env python

import sys
from natsort import natsorted
import argparse

VERSION = '1.0.5'

CHROMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]


def parseBed(handle, dlen):

    """
    :param handle:
    :return all_intervals:

    Take in a BED file, and return a dictionary like:
    
    """
    with handle as probes:
        all_intervals = {}
        for probe in probes:
            probe = probe.rstrip('\n')
            chrom = probe.split('\t')[0]
            start = int(probe.split('\t')[1]) + 1
            stop = int(probe.split('\t')[2])

            if chrom not in all_intervals:
                all_intervals[chrom] = []

            q0 = 0
            q30 = 0
            temp_list = [start, stop, q0, q30] + (dlen * [0])
            all_intervals[chrom].append(temp_list)

    return all_intervals


def parseGff3(handle):
    with handle as gff3:
        all_cds = {}
        gff_id_to_hgnc = {}
        refseq_gene = {}
        refseq_coords = {}

        for x in CHROMS:
            all_cds[x] = []

        for entry in gff3:
            if entry[0] != "#":

                entry_type = entry.split('\t')[2]
                start = int(entry.split('\t')[3])
                stop = int(entry.split('\t')[4])
                info = entry.split('\t')[8].split(';')

                chrom = int(entry.split('\t')[0].split('.')[0].split('_')[1])
                if chrom == 23:
                    chrom = "X"
                elif chrom == 24:
                    chrom = "Y"
                elif chrom == 12920:
                    chrom = "MT"
                chrom = str(chrom)

                if chrom in all_cds:

                    if (entry_type == "mRNA" or entry_type == "transcript"
                            or entry_type == "primary_transcript" or entry_type == "ncRNA"):
                        for annot in info:
                            gff_id = annot.split('=')[0]
                            gff_info = annot.split('=')[1]
                            if gff_id == "ID":
                                curr_gff_id = gff_info
                            elif gff_id == "Name":
                                curr_gff_name = gff_info
                            elif gff_id == "gene":
                                curr_gff_gene = gff_info

                        gff_id_to_hgnc[curr_gff_id] = [curr_gff_name, curr_gff_gene]
                        
                        # Special case for NR_ identified transcripts, like MIR's.
                        if curr_gff_name[:3] == 'NR_':
                            refseq_coords[curr_gff_name] = [[chrom, start, stop]]
                            refseq_gene[curr_gff_name] = curr_gff_gene


                    # Special case for MT genes, there are no associated RefSeq id's.
                    elif entry_type == "gene" and chrom == "MT":
                        for annot in info:
                            gff_id = annot.split('=')[0]
                            gff_info = annot.split('=')[1]
                            if gff_id == "Name":
                                refseq_gene[gff_info] = gff_info
                                if gff_info in refseq_coords:
                                    refseq_coords[gff_info].append([chrom, start, stop])
                                else:
                                    refseq_coords[gff_info] = [[chrom, start, stop]]

                                all_cds[chrom].append([start, stop])
                                break
                                

                    elif entry_type == "CDS":
                        for annot in info:
                            gff_id = annot.split('=')[0]
                            gff_info = annot.split('=')[1]
                            if gff_id == "Parent":
                                curr_gff_parent = gff_info
                                break # No need to continue looking for anything.
                        if curr_gff_parent in gff_id_to_hgnc:
                            refseq = gff_id_to_hgnc[curr_gff_parent][0]
                            hgnc = gff_id_to_hgnc[curr_gff_parent][1]
                            if [start, stop] not in all_cds[chrom]:
                                all_cds[chrom].append([start, stop])
                            refseq_gene[refseq] = hgnc
                            if refseq in refseq_coords:
                                refseq_coords[refseq].append([chrom, start, stop])
                            else:
                                refseq_coords[refseq] = [[chrom, start, stop]]
#                            all_cds[chrom].append([start, stop, refseq, hgnc, 0, 0, 0, 0, 0, 0, 0]) # Q0, Q30, D200, D100, D50, d20, D10
                        else:
                            print(curr_gff_parent + " does not have an associated mRNA entry in the GFF.")
#                            pass


    return all_cds, refseq_gene, refseq_coords


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
            chrom = str(line[0].split(':')[0])
            coord = str(line[0].split(':')[1])
            depth = int(line[1])
            if chrom not in doc_dict:
                doc_dict[chrom] = {}
            doc_dict[chrom][coord] = depth
            i += 1

    return doc_dict


def calcDepth(depth, total_depths, my_depths):

    """
    Calculate the depths at values as defined by DEPTHS globally.
    Input: depth, total_depths
    Output: updated total_depths
    """

    if len(total_depths) != len(my_depths):
        raise Exception("total_depths is not the right size, please check code.")

    for value in my_depths:
        if depth >= int(value):
            total_depths[my_depths.index(value)] += 1
    return total_depths


def writeHeader(report_out, pre_header, my_depths):
    """
    Write the header of a QC tsv.
    """
    for label in pre_header:
        if pre_header.index(label) == 0:
            report_out.write(label)
        else:
            report_out.write('\t' + label)
    report_out.write("\tAVGD\tQ30")
    for depth in my_depths:
        report_out.write("\tD" + depth)
    report_out.write('\n')


def writeReportQC(report_qc, total_bp, report_out, refseq_gene, transcripts):
    """
    Write to file the gene level QC metrics.
    """

    for refseq in natsorted(report_qc):

        if transcripts == False:

            if report_qc[refseq][0] != 0.0:
                q30 = str("{:.4}".format((report_qc[refseq][1]*100)/report_qc[refseq][0]))
            else:
                q30 = "0.00"

            avgd = str("{:.1f}".format(report_qc[refseq][0]/total_bp[refseq]))
            report_out.write('\t'.join([refseq, refseq_gene[refseq], avgd, q30]))

            for value in report_qc[refseq][2:]:
                to_write = str("{:.4}".format((value*100)/total_bp[refseq]))
                report_out.write('\t' + to_write)

            report_out.write('\n')

        # This is exceptionally bad form, but I will need to refactor this later.
        elif refseq in transcripts:

            if report_qc[refseq][0] != 0.0:
                q30 = str("{:.4}".format((report_qc[refseq][1]*100)/report_qc[refseq][0]))
            else:
                q30 = "0.00"

            avgd = str("{:.1f}".format(report_qc[refseq][0]/total_bp[refseq]))
            report_out.write('\t'.join([refseq, refseq_gene[refseq], avgd, q30]))

            for value in report_qc[refseq][2:]:
                to_write = str("{:.4}".format((value*100)/total_bp[refseq]))
                report_out.write('\t' + to_write)

            report_out.write('\n')
            
    report_out.close()


def writeProbeQC(all_intervals, position_refseq, doc_dict_0, doc_dict_30, probe_out, refseq_gene, my_depths, dlen, pf_bases_aligned):

    ## From GFF:
    ## [18461046, 18461153, 'XM_005260650.1', 'POLR3F', 0, 0, 0, 0, 0, 0, 0]
    ## {CHROM: [START, STOP, REFSEQ, HGNC, Q0, Q30, D200, D100, D50, D20, D10]}

    ## From BED:
    ## {CHROM: [START, STOP, Q0, Q30, D200, D100, D50, D20, D10]}

    probe_qc = {}
    total_counts = [0.0, 0.0, 0.0] # seqlen, q0, q30, d200, d100, d50, d20, d10
    total_counts.extend(dlen * [0.0])
    for chrom in all_intervals:
        if chrom in doc_dict_0:
            for probe in all_intervals[chrom]:
                probe_temp = [0.0, 0.0]
                probe_temp.extend(dlen * [0.0])
                temp_refseq = []
                start = probe[0]
                stop = probe[1]
                seqlen = stop - start + 1
                for i in range(int(start), int(stop)+1):
                    if str(i) in doc_dict_0[chrom]:
                        probe_temp[0] += doc_dict_0[chrom][str(i)]
                        probe_temp[1] += doc_dict_30[chrom][str(i)]
                        probe_temp[2:] = calcDepth(doc_dict_0[chrom][str(i)], probe_temp[2:], my_depths)
                    else:
                        pass

                    if str(i) in position_refseq[chrom]:
                        temp_refseq.extend(position_refseq[chrom][str(i)])
                if probe_temp[0] != 0.0:
                    q30 = str("{:.4}".format((probe_temp[1]*100)/probe_temp[0]))
                else:
                    q30 = "0.00"

                to_write = ';'.join(natsorted(set(temp_refseq)))
                write_genes = []
                for refseq in set(temp_refseq):
                    write_genes.append(refseq_gene[refseq])

                avgd = str("{:.1f}".format(probe_temp[0]/seqlen))
                
                # Start position off by one error creates mismatch with templates.  This needs
                # to be fixed, but we will provide old coordinates until it can be fixed.
                # RIGHT
                probe_out.write('\t'.join([chrom, str(start), str(stop), to_write, ';'.join(set(write_genes)), avgd, q30]))
                # WRONG
#                probe_out.write('\t'.join([chrom, str(start-1), str(stop), to_write, ';'.join(set(write_genes)), avgd, q30]))

                for value in probe_temp[2:]:
                    write_depths = str("{:.4}".format((value*100)/seqlen))
                    probe_out.write('\t' + write_depths)
                probe_out.write('\n')

            ### Increment totals.
                total_counts[0] += seqlen
                total_counts[1] += probe_temp[0]
                total_counts[2] += probe_temp[1]
                i = 3
                for value in probe_temp[2:]:
                    total_counts[i] += value
                    i += 1

    ### Write totals.
    # No, actually we don't want to write the totals any more...
#     avgd = str("{:.1f}".format(total_counts[1]/total_counts[0]))
#     q30 = str("{:.4}".format((total_counts[2]*100)/total_counts[1]))
#     # Calculate amplicon efficiency, defined as amplicon base coverage
#     # divided by percent bases aligned to genome.
#     print("Total probe aligned bases: " + str(total_counts[1]))
#     print("Total aligned bases to genome: " + str(pf_bases_aligned))
#     amp_eff = str("{:.4}".format((total_counts[1]*100)/pf_bases_aligned))
#     if float(amp_eff) > 100 or float(amp_eff) < 0:
#         amp_eff = '0'
# # Waiting to include this...
# #    if amp_eff > 100 or amp_eff < 0:
# #        raise Exception("This value should not ever be above 100 or below 0.")
#     probe_out.write('\t'.join(["\t\t\tTOTAL", amp_eff, avgd, q30]))
#     for value in total_counts[3:]:
#         write_depths = str("{:.4}".format((value*100)/total_counts[0]))
#         probe_out.write('\t' + write_depths)
#     probe_out.write('\n')
    
    probe_out.close()


def positionRefSeq(refseq_coords):

    ### Associate positions to the RefSeq ID's overlapping them.

    position_refseq = {}
    for refseq, coords in refseq_coords.items():
        for coord in coords:
            chrom = coord[0]
            if chrom not in position_refseq:
                position_refseq[chrom] = {}
            for i in range(int(coord[1]), int(coord[2])+1):
                if str(i) not in position_refseq[chrom]:
                    position_refseq[chrom][str(i)] = [refseq]
                else:
                    if refseq not in position_refseq[chrom][str(i)]:
                        position_refseq[chrom][str(i)].append(refseq)

    return position_refseq


def create_report_qc(all_intervals, position_refseq, doc_dict_0, doc_dict_30, dlen, my_depths):

    ### Develop the metrics for gene-level QC reporting.
    ### To do this, we iterate over probe coordinates, and check to see if they are in the position_refseq dict.

    report_qc = {}
    total_bp = {}
    j = 0
    for chrom in all_intervals:
        if chrom in position_refseq:
            for probe in all_intervals[chrom]:
                if j % 10000 == 0:
                    print(j)
                for i in range(int(probe[0]), int(probe[1])+1):
                    if str(i) in position_refseq[chrom]:
                        for refseq in position_refseq[chrom][str(i)]:
                            if refseq not in total_bp:
                                total_bp[refseq] = 1
                            else:
                                total_bp[refseq] += 1
                            if refseq not in report_qc:
                                report_qc[refseq] = [0.0, 0.0] + (dlen * [0.0])
                            if str(i) in doc_dict_0[chrom]:
                                report_qc[refseq][0] += doc_dict_0[chrom][str(i)]
                                report_qc[refseq][1] += doc_dict_30[chrom][str(i)]
                                report_qc[refseq][2:] = calcDepth(doc_dict_0[chrom][str(i)], report_qc[refseq][2:], my_depths)
                            else:
                                pass

        else:
            print("Chromosome " + chrom + " not in position_refseq.")
        j += 1

    return report_qc, total_bp


def findHoles(refseq_coords, refseq_gene, doc_dict_0, hole_handle):

    ### Find CDS regions not covered by our probes.
    hole_dict = {}
    hole_handle.write("REFSEQ\tCOVERED\n")
    for refseq in refseq_coords:
        hole_stats = [0.0, 0.0]
        for coord in refseq_coords[refseq]:
            chrom = coord[0]
            start = int(coord[1])
            stop = int(coord[2])
            seqlen = stop - start + 1
            hole_stats[0] += seqlen
            for i in range(start, stop+1):
                if str(i) in doc_dict_0[chrom]:
                    hole_stats[1] += 1

        covered = str("{:.1}".format(hole_stats[1]/hole_stats[0]))
        hole_dict[refseq] = [covered, refseq, refseq_gene[refseq]]

    for value in sorted(hole_dict.values(), key=lambda x: (x[2], x[1])):
        hole_handle.write('\t'.join([value[1], value[2], value[0], '\n']))

    hole_handle.close()


def parsePicardMetrics(filename):
    """
    Read the Picard Alignment Summary Metrics file and
    pull the PF_ALIGNED_BASES value from the PAIR section.
    """
    with open(filename, 'r') as metrics:
        for line in metrics:
            line = line.rstrip('\n').split('\t')
            if line[0] == "PAIR" or line[0] == "UNPAIRED":
                return int(line[7])
    return None 

def transcript_list_parse(filename):
    """
    If a list of RefSeq transcripts is passed, only print QC metrics for the given list.
    """
    transcripts = []
    with open(filename, 'r') as refseq:
        for line in refseq:
            transcripts.append(line.rstrip('\n'))
            
    return transcripts


def main():

### Full-size data set test inputs
#    gff_handle = open("/home/exacloud/lustre1/users/letaw/projs/find_probes_in_gff/ref_GRCh37.p13_top_level.gff3", 'r')
#    probe_handle = open("agilent.bed", 'r')
#    doc_0_handle = open("/home/users/letaw/lustre1/projs/find_probes_in_gff/dataset_19105.dat", 'r')
#    doc_30_handle = open("/home/users/letaw/lustre1/projs/find_probes_in_gff/dataset_19114.dat", 'r')

### DMD-only data set test inputs
    # gff_handle = open("dmd.gff3", 'r')
    # probe_handle = open("dmd_probes.bed", 'r')
    # doc_30_handle = open("dmd_doc_0.tsv", 'r')
    # doc_0_handle = open("dmd_doc_30.tsv", 'r')

### Output files

    parser = argparse.ArgumentParser(description='Parameters associated with coverage metric creation.')
    parser.add_argument('gff', help='GFF3 formatted file.  Prefer RefSeq based.')
    parser.add_argument('bed', help='Probe or interval list that metrics will be created from in BED format.')
    parser.add_argument('doc_q0', help='DepthOfCoverage per locus output at Q0.')
    parser.add_argument('doc_q30', help='DepthOfCoverage per locus output at Q30.')
    parser.add_argument('gene_out', help='Gene metrics output file.')
    parser.add_argument('probe_out', help='Probe metrics output file.')
    parser.add_argument('depth', nargs='+', help='Depth cutoffs.')
    parser.add_argument('picard_metrics', help='Picard Alignment Summary Metrics.')
    parser.add_argument('--transcripts', help='List of RefSeq transcript id\'s to include in output.')

    ### Optional, will be used to find percent coverage of CDS regions by our interval set.
    parser.add_argument('--hole_out', help='Percent covered CDS sequences based on input probe set.')
    args = parser.parse_args()

    gff_handle = open(args.gff, 'r')
    probe_handle = open(args.bed, 'r')
    doc_0_handle = open(args.doc_q0, 'r')
    doc_30_handle = open(args.doc_q30, 'r')
    report_out = open(args.gene_out, 'w')
    probe_out = open(args.probe_out, 'w')

    ### Save length of args.depth.
    dlen = len(args.depth)

    print("Creating dictionary to hold BED probes.")
    all_intervals = parseBed(probe_handle, dlen)
    print("Creating dictionary to hold GFF3 CDS definitions.")
    all_cds, refseq_gene, refseq_coords = parseGff3(gff_handle)
    print("Creating dictionary to hold DepthOfCoverage coordinate/depth pairs at Q0.")
    doc_dict_0 = createDocDict(doc_0_handle)
    print("Creating dictionary to hold DepthOfCoverage coordinate/depth pairs at Q30.")
    doc_dict_30 = createDocDict(doc_30_handle)
    print("Creating positionRefSeq, to associate each coordinate with a RefSeq ID.")
    position_refseq = positionRefSeq(refseq_coords)
    print("Creating report_qc and total_bp.  These will hold data for the report-level QC metrics.")
    report_qc, total_bp = create_report_qc(all_intervals, position_refseq, doc_dict_0, doc_dict_30, dlen, args.depth)
 
    pre_header = ["REFSEQ", "HGNC"]
    print("Writing the header for gene metrics.")
    writeHeader(report_out, pre_header, args.depth) 

    # Create a transcript list if a transcript file has been passed.
    if args.transcripts:
        transcripts = transcript_list_parse(args.transcripts)
    else:
        transcripts = False

    print("Writing report QC metrics.")
    writeReportQC(report_qc, total_bp, report_out, refseq_gene, transcripts)

    pre_header = ["CHROM", "START", "STOP", "REFSEQ", "HGNC"]
    print("Writing the header for probe metrics.")
    writeHeader(probe_out, pre_header, args.depth) 

    print("Import Picard Summary Metrics, and retrieve PF_ALIGNED_BASES.")
    pf_bases_aligned = parsePicardMetrics(args.picard_metrics)
    print("Create and write probe metrics to file.")
    
    writeProbeQC(all_intervals, position_refseq, doc_dict_0, doc_dict_30, probe_out, refseq_gene, args.depth, dlen, pf_bases_aligned)

    if args.hole_out:
        print("Finding holes in probe coverage.")
        hole_handle = open(args.hole_out, 'w')
        findHoles(refseq_coords, refseq_gene, doc_dict_0, hole_handle)


if __name__ == "__main__":
    main()
