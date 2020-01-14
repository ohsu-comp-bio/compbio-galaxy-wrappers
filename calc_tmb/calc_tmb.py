#!/usr/bin/env python

import argparse
import json
import vcf

VERSION = '0.4.4'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--infile', nargs='+', help='Input VCF')
    parser.add_argument('--outfile', help='Output JSON')
    parser.add_argument('--variants', help='Output list of variants')
    parser.add_argument('--gnomad_af', default=0.00001, type=float, help='GNOMAD VAF to filter list upon.  Values greater than this will not be utilized.')
    parser.add_argument('--m2_tlod', default=40.0, type=float, help='M2 TLOD score to filter list upon.  Values less than this will not be utilized.')
    parser.add_argument('--min_ab', default=0.05, type=float, help='Variant allele frequency (balance) to filter list upon.  Values less than this will not be utilized.')
    parser.add_argument('--min_ab_pon_ov', default=0.15, type=float, help='Variant allele frequency (balance) to filter list upon when FILTER contains PON_OV.  Values less than this will not be utilized.')
    parser.add_argument('--min_dp_fb', default=500, type=int, help='Minimum read depth for FreeBayes variants used in calculation.  Values less than this will not be utilized.')
    parser.add_argument('--min_dp_m2', default=250, type=int, help='Minimum read depth for Mutect2 variants used in calculation.  Values less than this will not be utilized.')
    parser.add_argument('--seq_space', default=0.61, type=float, help='Total genomic size of sequence targets, in Mb.')
    parser.add_argument('--send_placeholder', type=float, help="Lab isn't interested in our TMB calculation, but we need to send something so they can edit the value in CGD.")
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

class VcfRec(object):
    def __init__(self, rec):

        self.rec = rec
        self.chrom = str(rec.CHROM)
        self.coord = str(rec.POS)
        self.ab = float(rec.samples[0]['AF'])
        self.dp = int(rec.samples[0]['DP'])
        self.filt = str(rec.FILTER)

        try:
            self.gnomad = float(rec.INFO['gnomad.AF'][0])
        except:
            self.gnomad = 0.000000001

        try:
            self.snpeff = [x.split('|')[1] for x in rec.INFO['ANN']]
        except:
            self.snpeff = None

        try:
            self.tlod = float(rec.INFO['TLOD'][0])
        except:
            self.tlod = None

        try:
            self.hgnc = [x.split('|')[3] for x in rec.INFO['ANN']]
        except:
            self.hgnc = None

        self.bad_gene = self._check_gene()


    def _check_gene(self):
        """
        Check to see if the gene is in the list of ignored genes.
        :return:
        """
        gene_ignore = ('HLA-A', 'HLA-B', 'HLA-C')
        for entry in self.hgnc:
            if entry in gene_ignore:
                return True
        return False




    def check_snpeff(self):
        """
        Check that a term we care about is in the SnpEff string.
        # These terms should only be included if they are also missense.  Since we already check for missense, we don't need to include these.
        # spec_terms = ('structural_interaction_variant', 'protein_protein_contact')
        :return:
        """
        snpeff_terms = ('missense_variant','frameshift_variant','stop_gained','missense_variant&splice_region_variant',
                        'splice_acceptor_variant&intron_variant','splice_donor_variant&intron_variant',
                        'disruptive_inframe_deletion','conservative_inframe_deletion','conservative_inframe_insertion',
                        'frameshift_variant&splice_region_variant','frameshift_variant&stop_gained',
                        'disruptive_inframe_insertion','start_lost','stop_lost',
                        '5_prime_UTR_premature_start_codon_gain_variant',
                        'conservative_inframe_insertion&splice_region_variant','frameshift_variant&start_lost',
                        'splice_acceptor_variant&splice_region_variant&intron_variant',
                        'splice_donor_variant&splice_region_variant&intron_variant',
                        'stop_gained&conservative_inframe_insertion',
                        'stop_gained&disruptive_inframe_insertion&splice_region_variant')

        for entry in snpeff_terms:
            if entry in self.snpeff:
                return True

class WholeVcf(object):
    """

    """
    def __init__(self, filename, gnomad_af, m2_tlod, min_ab, min_ab_pon_ov, min_dp_fb, min_dp_m2):
        self.vcf_reader = vcf.Reader(open(filename, 'r'))
        self.gnomad_af = gnomad_af
        self.m2_tlod = m2_tlod
        self.min_ab = min_ab
        self.min_ab_pon_ov = min_ab_pon_ov
        self.min_dp_fb = min_dp_fb
        self.min_dp_m2 = min_dp_m2
        self.to_write = []
        self.tmb_cnt = 0.0
        self._find_matches()

    def _find_matches(self):
        """
        Find the variants that match the criteria.
        :return:
        """
        for record in self.vcf_reader:
            entry = VcfRec(record)
            if entry.gnomad:
                if entry.gnomad < self.gnomad_af:
                    if entry.check_snpeff():
                        if entry.tlod:
                            if entry.tlod > self.m2_tlod:
                                if 'PON' not in entry.filt:
                                    if not entry.bad_gene:
                                        if 'm2' in entry.filt:
                                            if entry.dp > self.min_dp_m2:
                                                if 'PON_OV' in entry.filt:
                                                    if entry.ab > self.min_ab_pon_ov:
                                                        self.tmb_cnt += 1
                                                        self.to_write.append(record)
                                                else:
                                                    if entry.ab > self.min_ab:
                                                        self.tmb_cnt += 1
                                                        self.to_write.append(record)
                                        elif 'fb' in entry.filt or 'm2_fb' in entry.filt:
                                            if entry.dp > self.min_dp_fb:
                                                if 'PON_OV' in entry.filt:
                                                    if entry.ab > self.min_ab_pon_ov:
                                                        self.tmb_cnt += 1
                                                        self.to_write.append(record)
                                                else:
                                                    if entry.ab > self.min_ab:
                                                        self.tmb_cnt += 1
                                                        self.to_write.append(record)

class VcfCollector(object):
    """

    :param object:
    :return:
    """
    def __init__(self, seq_space, vcfs):
        self.seq_space = seq_space
        self.vcfs = vcfs
        self.tmb_cnt = self._collect_tmbs()
        self.tmb = self._calc_tmb()
        self.vars = self._collect_vars()

    def _calc_tmb(self):
        """
        Perform the TMB calculation
        :return:
        """
        return '{:.1f}'.format(self.tmb_cnt / self.seq_space)

    def _collect_tmbs(self):
        """
        Get all of the TMB values and add them together.
        :return:
        """
        tmb = 0.0
        for entry in self.vcfs:
            tmb += entry.tmb_cnt
        return tmb

    def _collect_vars(self):
        """
        Get all the variants from some number of vcfs.
        :return:
        """
        vars = []
        for entry in self.vcfs:
            vars.extend(entry.to_write)
        return vars

class Writer(object):
    """

    """
    def __init__(self, vcf_reader, outfile, vcf_out, tmb, vars):
        self.outfile = outfile
        self.vcf_out = vcf_out
        self.tmb = tmb
        self.vcf_reader = vcf_reader
        self.vars = vars

    def write_json_out(self):
        """
        Prepare output json file.
        :return:
        """
        outfile = open(self.outfile, 'w')
        out_metric = {'tmb': self.tmb}
        json.dump(out_metric, outfile)
        outfile.close()

    def write_vcf_out(self):
        """
        Write the VCF output.
        :return:
        """
        vcf_writer = vcf.Writer(open(self.vcf_out, 'w'), self.vcf_reader)
        for record in self.vars:
            vcf_writer.write_record(record)


def main():
    args = supply_args()
    vars = [WholeVcf(vcf, args.gnomad_af, args.m2_tlod, args.min_ab, args.min_ab_pon_ov, args.min_dp_fb, args.min_dp_m2) for vcf in args.infile]
    all_vcfs = VcfCollector(args.seq_space, vars)
    if not args.send_placeholder:
        writer = Writer(vars[0].vcf_reader, args.outfile, args.variants, all_vcfs.tmb, all_vcfs.vars)
    else:
        writer = Writer(vars[0].vcf_reader, args.outfile, args.variants, args.send_placeholder, all_vcfs.vars)

    writer.write_json_out()
    writer.write_vcf_out()

if __name__ == "__main__":
    main()
