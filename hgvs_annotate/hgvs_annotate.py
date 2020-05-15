#!/usr/bin/env python

# DESCRIPTION: Apply HGVS nomenclature to an input VCF, and place in INFO field.
# https://github.com/biocommons/hgvs/
# These libraries most closely represent the current HGVS recommendations, imo.
# USAGE: hgvs_test.py <input> <output>
# CODED BY: John Letaw

import argparse
import multiprocessing as mp
import vcfpy
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser

VERSION = '0.4.0'

CHROM_MAP = {'1': 'NC_000001.10', '2': 'NC_000002.11', '3': 'NC_000003.11', '4': 'NC_000004.11', '5': 'NC_000005.9',
             '6': 'NC_000006.11', '7': 'NC_000007.14', '8': 'NC_000008.10', '9': 'NC_000009.11', '10': 'NC_000010.10',
             '11': 'NC_000011.9', '12': 'NC_000012.11', '13': 'NC_000013.10', '14': 'NC_000014.8', '15': 'NC_000015.9',
             '16': 'NC_000016.9', '17': 'NC_000017.10', '18': 'NC_000018.9', '19': 'NC_000019.9', '20': 'NC_000020.10',
             '21': 'NC_000021.8', '22': 'NC_000022.10', 'X': 'NC_000023.10', 'Y': 'NC_000024.9', 'MT': 'NC_012920.1'}

AA_MAP = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Asx': 'B', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q', 'Glx': 'Z',
          'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
          'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'}


def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input_vcf', help='Input VCF.')
    parser.add_argument('--output_vcf', help='Output VCF.')

    parser.add_argument('--chrom', help='Chromosome.')
    parser.add_argument('--pos', help='Position.')
    parser.add_argument('--ref', help='Refence Allele.')
    parser.add_argument('--alt', help='Alternate Allele.')

    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    if args.input_vcf and not args.output_vcf:
        raise Exception('Must specify an output VCF if you are specifying an input VCF.')
    if args.chrom and (not args.pos or not args.ref or not args.alt):
        raise Exception('Must specify each of chrom, pos, ref, and alt.')

    return args




class VcfRec:
    """
    Information corresponding to a single VCF entry.
    """
    def __init__(self, rec, hp):
        self.hp = hp
        self.chrom = CHROM_MAP[rec.CHROM]
        self.pos = rec.POS
        self.ref = rec.REF
        self.alt = [x.value for x in rec.ALT]
        self.hgvs_g = self._prep_hgvs()
        self.var_g = self._create_var_g()

    def _create_var_g(self):
        var_g = []
        for entry in self.hgvs_g:
            var_g.append(self.hp.parse_hgvs_variant(entry))
        return var_g

    def _prep_hgvs(self):
        """

        :return:
        """
        hgvs_g = []
        for entry in self.alt:
            pos_part = self.correct_indel_coords(self.pos, self.ref, str(entry))
            new_hgvs = self.chrom + ':g.' + pos_part
            hgvs_g.append(new_hgvs)
        return hgvs_g

    @staticmethod
    def correct_indel_coords(pos, ref, alt):
        """
        Using a VCF position, create coords that are compatible with HGVS nomenclature.
        Since we are already determining at this stage whether the event is an ins or del, also
        include the ins or del strings in the result.
        substitution event -> ac:g.[pos][ref]>[alt]
        :return:
        """
        lref = len(ref)
        lalt = len(alt)
        if lref == 1 and lalt == 1:
            # Substitution case
            change = '>'.join([ref, alt])
            new_pos = str(pos) + change
            return new_pos
        elif lref == lalt:
            # Multi-nucleotide substitution case
            # NG_012232.1: g.12_13delinsTG
            new_start = str(pos)
            new_end = str(int(pos) + lref - 1)
            new_pos = '_'.join([new_start, new_end]) + 'delins' + alt
            return new_pos
        elif lref > lalt:
            # Deletion case
            shift = lref - lalt
            if shift == 1:
                new_pos = str(int(pos) + 1) + 'del'
                return new_pos
            else:
                new_start = str(int(pos) + 1)
                new_end = str(int(pos) + shift)
                new_pos = '_'.join([new_start, new_end]) + 'del'
                return new_pos
        elif lalt > lref:
            # Insertion case
            new_start = str(pos)
            new_end = str(int(pos) + 1)
            new_pos = '_'.join([new_start, new_end]) + 'ins' + alt[1:]
            return new_pos
        else:
            # OTHER case
            print("Change type not supported: " + pos + ':' + ref + '>' + alt + '\n')


class VcfReader:
    """

    """
    def __init__(self, filename, hp):
        self.hp = hp
        self.vcf_reader = vcfpy.Reader.from_path(filename)
        self.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'HGVS_G'),
                                                           ('Number', '.'),
                                                           ('Type', 'String'),
                                                           ('Description', 'HGVS genomic reference, as produced by pyhgvs.')]))
        self.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'HGVS_C'),
                                                           ('Number', '.'),
                                                           ('Type', 'String'),
                                                           ('Description', 'HGVS coding reference, as produced by pyhgvs.')]))
        self.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'HGVS_P'),
                                                           ('Number', '.'),
                                                           ('Type', 'String'),
                                                           ('Description', 'HGVS protein reference, as produced by pyhgvs.')]))
        self.my_vrnts = self._create_rec()


    def _create_rec(self):
        """
        """
        my_vrnts = []
        for vrnt in self.vcf_reader:
            my_vrnts.append(VcfRec(vrnt, self.hp))
        return my_vrnts


        # for entry in alt:
        #     var_g = hp.parse_hgvs_variant(new_hgvs)
        #     # This is how you find out which transcripts are available.  Might as well provide HGVS for all of them.
        #
        #     if 'HGVS_G' not in record.INFO:
        #         record.INFO['HGVS_G'] = [str(var_g)]
        #     else:
        #         record.INFO['HGVS_G'].append(str(var_g))
        #
        #     tx_list = hdp.get_tx_for_region(str(var_g.ac), 'splign',
        #                                     str(var_g.posedit.pos.start), str(var_g.posedit.pos.end))
        #     for tx in tx_list:
        #         try:
        #             var_c = vm.g_to_c(var_g, str(tx[0]))
        #             if 'HGVS_C' not in record.INFO:
        #                 record.INFO['HGVS_C'] = [str(var_c)]
        #             else:
        #                 record.INFO['HGVS_C'].append(str(var_c))
        #
        #             var_p = vm.c_to_p(var_c)
        #             if 'HGVS_P' not in record.INFO:
        #                 record.INFO['HGVS_P'] = [str(var_p)]
        #             else:
        #                 record.INFO['HGVS_P'].append(str(var_p))
        #
        #         except:
        #             print("Transcript probably non-coding, ignore.")
        #
        # return record

def get_txs(my_vcf, hdp):
    """
    
    :return: 
    """
    for rec in my_vcf:
        for var_g in rec.var_g:
            print(hdp.get_tx_for_region(str(var_g.ac), 'splign', str(var_g.posedit.pos.start), str(var_g.posedit.pos.end)))

def main():

    args = supply_args()

    # Initialize the HGVS package objects.
    hdp = hgvs.dataproviders.uta.connect()
    vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method='splign')
    hp = hgvs.parser.Parser()
    pool = mp.Pool(16)

    if args.input_vcf:
        # vcf_writer = vcfpy.Writer.from_path(args.output_vcf, vcf_reader.header)
        my_vcf = VcfReader(args.input_vcf, hp).my_vrnts
        results = pool.starmap(get_txs, [(x for x in my_vcf), hdp])

        # for entry in results:
        #     print(entry)
        #     raise Exception("ASHIFS")

        # my_hgvs = HgvsVcfWriter(args.input_vcf, args.output_vcf, hp, hdp, vm)
        # results = create_hgvs(my_hgvs)

    # if args.chrom:
    #     chrom = CHROM_MAP[args.chrom]
    #     pos = args.pos
    #     ref = args.ref
    #     alt = args.alt
    #     pos_part = correct_indel_coords(pos, ref, alt)
    #     new_hgvs = chrom + ':g.' + pos_part
    #     var_g = hp.parse_hgvs_variant(new_hgvs)
    #     # This is how you find out which transcripts are available.  Might as well provide HGVS for all of them.
    #
    #     tx_list = hdp.get_tx_for_region(str(var_g.ac), 'splign',
    #                                     str(var_g.posedit.pos.start), str(var_g.posedit.pos.end))
    #     for entry in tx_list:
    #         print(entry[0])
    #         print(var_g)
    #         try:
    #             var_c = vm.g_to_c(var_g, entry[0])
    #             print(str(var_c))
    #             var_p = vm.c_to_p(var_c)
    #             print(str(var_p))
    #         except:
    #             print("NON-CODING LIKELY.")


if __name__ == "__main__":
    main()
