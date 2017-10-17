#!/usr/bin/env python

### Add a line to the samples file for the tumor sample.
### D707-D505-S11-16979	D707.D505.S11.16979.normal	Normal
### Where it says normal, we will replace with tumor.  Sample ID should be passed as it is
### with other tools.
### Add an entry to each line of the counts file for the tumor sample.
#AmpliconId,D705.D502.S09.04850.normal,D705.D503.S09.05423.normal,D705.D504.S09.18738.normal,D705.D505.S10.01560.normal,D705.D506.S10.01741.normal,D705.D507.S10.04808.normal,D705.D508.S10.10744.normal,D706.D501.S10.11021.normal,D706.D502.S10.12512.normal,D706.D503.S10.14323.normal,D706.D504.S10.19095.normal,D706.D505.S11.04719.normal,D706.D506.S11.04807.normal,D706.D507.S11.05123.normal,D706.D508.S11.06702.normal,D707.D501.S11.09562.normal,D707.D502.S11.09798.normal,D707.D503.S11.12220.normal,D707.D504.S11.12495.normal,D707.D505.S11.16979.normal
#1950341,930,1244,1134,1045,1119,1097,737,1197,1819,875,1102,920,812,722,724,922,1417,934,586,664
### This script will be run from Galaxy.

VERSION='0.1.0'

import sys


def parseProbeQC(probeqc_file, probes):
    """
    Remove depth field from ProbeQC coverage metrics.  Match depth by genomic coordinates.
    """

    with probeqc_file as probeqc:
        next(probeqc)
        for line in probeqc:
            if "TOTAL" not in line:
                line = line.rstrip('\n').split('\t')
                chrom = sex_chrom_switch(line[0])
                # The "+1" will be removed after fixing templates to correct coordinates.
                # RIGHT
                start = int(line[1])
                # WRONG
#                start = int(line[1])+1
                stop = line[2]
                avgd = line[5].split('.')[0]
        
                coord = chrom + ':' + str(start) + '-' + stop
                if coord not in probes:
                    probes[coord] = [avgd]
                else:
                    probes[coord].append(avgd)

    return probes


def sex_chrom_switch(chrom):
    """
    Convert X and Y chroms to 23 and 24.
    """
    if chrom == "X":
        return "23"
    elif chrom == "Y":
        return "24"
    else:
        return chrom


def parse_cnv_amplicons(handle_amps):
    """
    """
    cnv_amps = {}

    with handle_amps as amplicons:
        next(amplicons)
        for line in amplicons:
            line = line.rstrip('\n').split('\t')
            chrom = sex_chrom_switch(line[2])
            coord = chrom + ':' + line[3] + '-' + line[4]
            cnv_amps[coord] = line[0]

    return cnv_amps


def convert_sample_id(sample_id):
    """
    Convert hyphens to dots in sample_id's.
    """
    new_sample_id = '.'.join(sample_id.split('-')) + '.tumor'
    return new_sample_id


def read_and_write_counts(infile, outfile, sample_id, coords_to_id, probe_depth):
    """
    Read the counts file, and write another counts file.
    """

    with infile as counts:
        for entry in counts:
            new_entry = entry.rstrip('\n')
            new_entry += ','
            if entry.split(',')[0] == "AmpliconId":
                new_entry += convert_sample_id(sample_id)
            else:
                ampl_id = entry.split(',')[0]
                for key in coords_to_id:
                    if coords_to_id[key] == ampl_id:
                        new_entry += probe_depth[key][0]
                        continue

            new_entry += '\n'
            outfile.write(new_entry)

    outfile.close()


def write_sample(infile, outfile, sample_id):
    """
    """
    with infile as samples:
        for line in samples:
            outfile.write(line)
    outfile.write('\t'.join([sample_id, convert_sample_id(sample_id), "Tumor"]))
    outfile.write('\n')

    outfile.close()


def main():

    probes = {}

    handle_amplicons = open(sys.argv[1], 'rU')
    handle_counts = open(sys.argv[2], 'rU')
    handle_samples = open(sys.argv[3], 'rU')

    handle_counts_out = open(sys.argv[4], 'w')
    handle_samples_out = open(sys.argv[5], 'w')

    handle_probes = open(sys.argv[6], 'rU')
    sample_id = sys.argv[7]

    coords_to_id = parse_cnv_amplicons(handle_amplicons)
    probe_depth = parseProbeQC(handle_probes, probes)

    read_and_write_counts(handle_counts, handle_counts_out, sample_id, coords_to_id, probe_depth)
    write_sample(handle_samples, handle_samples_out, sample_id)

if __name__ == "__main__":
    main()
