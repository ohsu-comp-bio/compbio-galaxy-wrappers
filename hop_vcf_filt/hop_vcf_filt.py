#!/usr/bin/env python

# GATK4 SelectVariants is not working in the way we need to properly filter these VCFs.  Instead of spending a bunch
# of time trying to get it to work for us, I'm writing this short script to do custom filtering.  It will only be useful
# for the purpose of post-filtering multi-sample HOP VCFs.

import argparse
import numpy
import vcf

VERSION = '0.2.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--infile', help='Input VCF')
    parser.add_argument('--infile_genes', help='Input HGNC Gene List')
    parser.add_argument('--outfile', help='Output VCF')
    parser.add_argument('--outfile_bad', help='Filtered Sites Output')
    parser.add_argument('--chrom', help='Output filtering details for coordinate.')
    parser.add_argument('--coord', help='Output filtering details for coordinate.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

class VcfRec(object):
    def __init__(self, rec):

        self.rec = rec
        self.chrom = str(rec.CHROM)
        self.coord = str(rec.POS)

        try:
            self.clnsig = rec.INFO['clinvar.CLNSIG'][0]
        except:
            self.clnsig = None

        try:
            self.clnsigconf = ','.join(rec.INFO['clinvar.CLNSIGCONF'])
        except:
            self.clnsigconf = ''

        try:
            self.gnomad = float(rec.INFO['gnomad.AF'][0])
        except:
            self.gnomad = None

        try:
            self.hgmd = rec.INFO['hgmd.CLASS'][0]
        except:
            self.hgmd = None

        self.snpeff_filt = ["3_prime_UTR_variant", "5_prime_UTR_variant", "downstream_gene_variant", "intron_variant", "intergenic_region", "synonymous_variant", "upstream_gene_variant", "non_coding_transcript_exon_variant"]

        self.snpeff = False
        try:
            for entry in rec.INFO['ANN']:
                if entry.split('|')[1] not in self.snpeff_filt:
                    self.snpeff = True
        except:
            self.snpeff = False

        try:
            self.snpeff_gene = rec.INFO['ANN'][0].split('|')[3]
        except:
            self.snpeff_gene = ''

        try:
            self.af = float(rec.INFO['AF'][0])
        except:
            self.af = None

        self.qual = float(rec.QUAL)

        self.is_path = (self.clnsig == 'Pathogenic') or \
                       (self.clnsig == 'Likely_pathogenic') or \
                       ('athogenic' in self.clnsigconf) or \
                       (self.hgmd == 'DM')

        self.no_info = not self.clnsig and not self.gnomad and not self.hgmd

        self.ab_avg, self.ab_std = self._calc_ab()
        if self.ab_avg:
            if self.ab_avg + (3 * self.ab_std) < 0.35:
                self.bad_ab = True
            elif self.ab_avg < 0.25 and self.af > 0.05:
                self.bad_ab = True
            else:
                self.bad_ab = False
        else:
            self.bad_ab = False

    def _calc_ab(self):
        """
        Get the average allele balance for a specific genotype.
        :return:
        """
        all_vafs = []
        for entry in self.rec.samples:
            if entry['GT'] == '0/1':
                ref_cnt = entry['AD'][0]
                alt_cnt = entry['AD'][1]
                try:
                    vaf = alt_cnt / (alt_cnt + ref_cnt + 0.0)
                except:
                    vaf = 0.0
                all_vafs.append(vaf)

        if len(all_vafs) > 1:
            try:
                return numpy.mean(all_vafs), numpy.std(all_vafs)
            except:
                return None, None
        return None, None

    def var_req(self, schrom, scoord):
        """
        Print out the current value of each variable.
        :param self:
        :return:
        """
        if schrom == self.chrom and scoord == self.coord:
            print("Chromosome: {0}".format(self.chrom))
            print("Position: {0}".format(self.coord))
            print("CLNSIG: {0}".format(self.clnsig))
            print("CLNSIGCONF: {0}".format(self.clnsigconf))
            print("GNOMAD: {0}".format(self.gnomad))
            print("HGMD: {0}".format(self.hgmd))
            print("SNPEFF STATUS: {0}".format(self.snpeff))
            print("SNPEFF GENE: {0}".format(self.snpeff_gene))
            print("LOCAL AF: {0}".format(self.af))
            print("QUAL: {0}".format(self.qual))
            print("IS PATH?: {0}".format(self.is_path))
            print("NO INFO?: {0}".format(self.no_info))
            print("AVG AB: {0}".format(self.ab_avg))
            print("STDEV AB: {0}".format(self.ab_std))
            print("BAD AB?: {0}".format(self.bad_ab))
            for entry in self.rec.INFO['ANN']:
                print(entry)


def process_genes(filename):
    """
    Put genes in a list from file.
    :return:
    """
    genes = []
    with open(filename, 'rU') as myfile:
        for line in myfile:
            genes.append(line.rstrip('\n'))
    return genes

def main():
    args = supply_args()
    vcf_reader = vcf.Reader(open(args.infile, 'r'))
    vcf_writer = vcf.Writer(open(args.outfile, 'w'), vcf_reader)
    vcf_writer_bad = vcf.Writer(open(args.outfile_bad, 'w'), vcf_reader)
    ingenes = process_genes(args.infile_genes)

    for record in vcf_reader:
        entry = VcfRec(record)
        if args.chrom and args.coord:
            entry.var_req(args.chrom, args.coord)
        if (entry.snpeff_gene) not in ingenes:
            vcf_writer_bad.write_record(record)
        elif entry.bad_ab:
            vcf_writer_bad.write_record(record)
        elif entry.is_path:
            vcf_writer.write_record(record)
        elif (entry.clnsig == 'Benign') or \
                (entry.clnsig == 'Likely_benign') or \
                (entry.clnsig == 'Benign/Likely_benign'):
            vcf_writer_bad.write_record(record)
        elif entry.clnsig == 'Conflicting_interpretations_of_pathogenicity':
            vcf_writer_bad.write_record(record)
        elif (not entry.snpeff) and (entry.clnsig == 'Uncertain_significance'):
            vcf_writer_bad.write_record(record)
        elif (not entry.snpeff) and (not entry.clnsig):
            vcf_writer_bad.write_record(record)
        elif (not entry.snpeff) and (entry.clnsig == 'not_provided'):
            vcf_writer_bad.write_record(record)
        elif entry.gnomad:
            if entry.gnomad < 0.02:
                vcf_writer.write_record(record)
            else:
                vcf_writer_bad.write_record(record)
        elif entry.no_info and entry.qual < 100.0:
            vcf_writer_bad.write_record(record)
        elif entry.af:
            if entry.clnsig == 'not_provided' and entry.af < 0.02:
                vcf_writer.write_record(record)
            elif entry.clnsig == 'not_provided' and entry.af >= 0.02:
                vcf_writer_bad.write_record(record)
            elif entry.af < 0.05:
                vcf_writer.write_record(record)
            else:
                vcf_writer_bad.write_record(record)
        else:
            vcf_writer_bad.write_record(record)

    vcf_writer.close()
    vcf_writer_bad.close()

if __name__ == "__main__":
    main()
