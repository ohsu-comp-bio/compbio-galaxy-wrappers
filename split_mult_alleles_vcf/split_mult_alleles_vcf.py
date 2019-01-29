#!/usr/bin/env python

# DESCRIPTION: VCF files list multiple alternate alleles in the ALT columns,
#  separated by commas. This script will split entries with multiple alleles
# listed, and place each on a separate line.  This will allow us to import
# this data in to annotation software and properly receive prediction
# scores. Post-normalization (vt, bcftools) of this data is highly recommended.
# USAGE: split_mult_alleles.py -h
# CODED BY: John Letaw

import argparse
import itertools
import pysam

VERSION = '0.4.0'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument(dest='input', help='')
    parser.add_argument(dest='output', help='')
    parser.add_argument(dest='ref', help='')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
                                                               VERSION)
    args = parser.parse_args()
    return args


def write_vcf_line(handle, line, allele, alt, style):

    # 1	2488153	.	A	G	16490.6	.	AB=0.504146;ABP=3.19036;AC=1;AF=0.5;AN=2;AO=608;CIGAR=1X;DP=1206;DPB=1206;DPRA=0;EPP=71.026;EPPR=73.2867;GTI=0;LEN=1;MEANALT=2;MQM=59.9062;MQMR=59.943;NS=1;NUMALT=1;ODDS=3721.51;PAIRED=0;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=23566;QR=23045;RO=597;RPL=535;RPP=765.326;RPPR=684.965;RPR=73;RUN=1;SAF=252;SAP=41.6397;SAR=356;SRF=267;SRP=17.4468;SRR=330;TYPE=snp;technology.ILLUMINA=1	GT:DP:AD:RO:QR:AO:QA:GL:AF:PUMI	0/1:1206:597,608:597:23045:608:23566:-1753.83,0,-1708.68:0.5045643:0.542
    if style == 'total_alt':
        counts = line[9].split(':')[2].split(',')[0] + ',' + line[9].split(':')[2].split(',')[alt]
        geno = ':'.join(line[9].split(':')[3:])
    else:
        my_i = find_index(line[8], ':', 'AD')
        counts = line[9].split(':')[my_i].split(',')[0] + ',' + line[9].split(
            ':')[my_i].split(',')[alt]
        geno = ':'.join(line[9].split(':')[2:])

    for i in range(len(line)):
        if i == len(line)-1:
            handle.write("0/1:" + counts + ":" + geno + '\n')
        elif i == 4:
            handle.write(allele + '\t')
        elif i != 4:
             handle.write(line[i] + '\t')


def find_index(splitting, delim, to_find):
    """
    Find the index of the field in the FORMAT column you care about. Usually
    will be AD.
    """

    for entry in splitting.split(delim):
        if entry == to_find:
            curr_index = splitting.split(delim).index(to_find)
            return curr_index

    return None


def calc_gl_pos(ref, alt):
    """
    F(j/k) = (k*(k+1)/2)+j
    j = allele 1
    k = allele 2
    :return:
    """
    return ((ref*(ref+1)/2)+alt)


def geno_prob_parse(alt_cnt):
    """
    Special case of the GL FORMAT field and multiple alt alleles. This will
    return a dict of genotype tuple to position placement:
    ex: {(0,1):1}
    :return:
    """
    values = ''
    genos = {}
    for x in range(alt_cnt+1):
        values += str(x)
    for entry in itertools.combinations_with_replacement(values, 2):
        if entry not in genos:
            genos[entry] = calc_gl_pos(int(entry[1]), int(entry[0]))
    return sorted(genos.items(), key=lambda x: x[1])


def include_gl(genos, allele):
    """
    Given a genotype number, determine whether this entry should be split
    out and included in the individual entry for this ref/alt combo.
    :return:
    """
    gl_ind = []
    for geno in genos:
        if geno[0] == ('0', '0') or geno[0] == ('0', str(allele)) or geno[0]\
                == (str(allele), str(allele)):
            gl_ind.append(geno[1])

    return gl_ind


def info_break(info, allele):
    """
    Break the INFO field up by allele.
    :return:
    """
    new_info = []
    for entry in info.split(';'):
        label = entry.split('=')[0]
        try:
            value = entry.split('=')[1].split(',')
            if len(value) > 1:
                new_info.append('='.join([label, value[allele]]))
            else:
                new_info.append(entry)
        except IndexError:
            new_info.append(entry)

    return ';'.join(new_info)


def sample_break(format, sample):
    """
    Break the SAMPLE field up by alt allele.
    ex:
    POS=27057923, REF=A, ALT=C,G
    GT:DP:AD:RO:QR:AO:QA:GL
    0/0:312:230,46,35:230:7081:46,35:671,513:0,-22.7179,-573.147,-33.6276,-551.406,-587.363
    AD = [230,46,35], AO = [46,35], QA = [671,513]
    GL = [0,-22.7179,-573.147,-33.6276,-551.406,-587.363]
    :return:
    """
    new_info = {}
    sample = sample.split(':')
    format = format.split(':')
    for entry in format:
        curr_ind = format.index(entry)
        new_info[entry] = sample[curr_ind]

    return new_info


def ind_to_geno(selec, value):
    """
    Convert an index value to a genotype.
    :return:
    """
    if selec.index(value) == 0:
        return '0/0'
    if selec.index(value) == 1:
        return '0/1'
    if selec.index(value) == 2:
        return '1/1'
    return './.'


def new_gt(broke_samp, gl_id, label):
    """
    Redefine the GT field in multiallelic entries.
    :return:
    """
    gl = broke_samp[label].split(',')
    this_gl = [float(gl[x]) for x in gl_id]
    for entry in this_gl:
        if entry == 0.0:
            return ind_to_geno(this_gl, entry)
        else:
            # low_val = this_gl.index(min(this_gl, key=abs))
            low_val = min(this_gl, key=abs)
            return ind_to_geno(this_gl, low_val)

    raise Exception("Why are we not returning anything here? ")


def collect_gls(gl_ind, broke_samp, label):
    """
    Take the gl_ind index list and return the actual associated values from
    the GL field.
    :param gl_ind:
    :return:
    """
    temp = []
    old_gls = broke_samp[label].split(',')
    for value in gl_ind:
        temp.append(old_gls[value])

    return ','.join(temp)

def proc_samp():
    """
    Get the sample section written out correctly for each sample in the VCF.
    :return:
    """
    pass

def main():
    
    args = supply_args()

    handle_in_vcf = open(args.input, 'rU')
    handle_out_vcf = open(args.output, 'w')
    # broke_samp = []

    with handle_in_vcf as vcf:
        for line in vcf:
            if not line.startswith('#'):
                new_line = line.rstrip('\n').split('\t')
                chrom = new_line[0]
                pos = new_line[1]
                rsid = new_line[2]
                ref = new_line[3]
                alts = new_line[4]
                qual = new_line[5]
                filter = new_line[6]
                info = new_line[7]
                samples = new_line[9:]

                if ',' in alts:
                    alt_allele = alts.split(',')
                    genos = geno_prob_parse(len(alt_allele))
                    for i in range(1, len(alt_allele)+1):
                        # Assess alt allele here, if asterisk.
                        if alt_allele[i-1] == '*':
                            extra_base = pysam.FastaFile(args.ref).fetch(chrom, int(pos)-2, int(pos)-1)
                            pos = str(int(pos) - 1)
                            ref = extra_base + ref
                            alt_allele[i-1] = extra_base
                        to_write = [chrom, pos, rsid, ref]
                        gl_ind = include_gl(genos, i)
                        to_write.append(alt_allele[i-1])
                        to_write.extend([qual, filter])
                        to_write.append(info_break(info, i-1))

                        try:
                            format = new_line[8]
                            to_write.append(format)

                            for sample in samples:
                                broke_samp = sample_break(format, sample)
                                # Work up the SAMPLE section.
                                # GT:DP:AD:RO:QR:AO:QA:GL
                                # 0/1:1206:597,608:597:23045:608:23566:-1753.83,0,-1708.68:0.5045643:0.542
                                new_samp = []
                                for field in format.split(':'):
                                    if field == 'GT':
                                        if 'GL' in broke_samp:
                                            new_field = new_gt(broke_samp, gl_ind, 'GL')
                                        elif 'PL' in broke_samp:
                                            new_field = new_gt(broke_samp, gl_ind, 'PL')
                                        else:
                                            if broke_samp['GT'] != '1/1':
                                                new_field = '0/1'
                                            else:
                                                new_field = '1/1'
                                    elif field == 'AD' or field == 'F1R2' or field == 'F2R1' or field == 'MBQ' or field == 'MFRL':
                                        this_samp_ad_ref = broke_samp[field].split(',')[0]
                                        this_samp_ad_alt = broke_samp[field].split(',')[i]
                                        new_field = ','.join([this_samp_ad_ref, this_samp_ad_alt])
                                    elif field == 'AO' or field == 'QA' or field == 'AF':
                                        new_field = broke_samp[field].split(',')[i-1]
                                    elif field == 'GL' or field == 'PL':
                                        new_field = collect_gls(gl_ind, broke_samp, field)
                                    else:
                                        new_field = broke_samp[field]

                                    new_samp.append(new_field)
                                to_write.append(':'.join(new_samp))
                        except:
                            pass
                        handle_out_vcf.write('\t'.join(to_write))
                        handle_out_vcf.write('\n')
                else:
                    handle_out_vcf.write(line)
            else:
                handle_out_vcf.write(line)

    handle_out_vcf.close()

if __name__ == "__main__":
    main()
