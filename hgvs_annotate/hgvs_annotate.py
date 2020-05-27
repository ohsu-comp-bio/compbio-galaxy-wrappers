#!/usr/bin/env python

# DESCRIPTION: Apply HGVS nomenclature to an input VCF, and place in INFO field.
# https://github.com/biocommons/hgvs/
# These libraries most closely represent the current HGVS recommendations, imo.
# USAGE: hgvs_test.py <input> <output>
# CODED BY: John Letaw

import argparse
import dill
import vcfpy
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.easy
import hgvs.parser
from pathos.multiprocessing import ProcessingPool

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
    def __init__(self, rec):
        # self.hp = hp
        # self.hdp = hdp
        # self.vm = vm
        self.hp = hgvs.parser.Parser()
        self.hdp = hgvs.dataproviders.uta.connect()
        self.vm = hgvs.assemblymapper.AssemblyMapper(self.hdp, assembly_name="GRCh37", alt_aln_method='splign')
        self.chrom = CHROM_MAP[rec.CHROM]
        self.pos = rec.POS
        self.ref = rec.REF
        self.alt = [x.value for x in rec.ALT]
        self.hgvs_g = self._prep_hgvs()
        self.var_g = self._create_var_g()
        self.uniq_key = (self.chrom, self.pos, self.ref, self.alt[0])
        self.tx_list = self._get_txs()
        self.var_c = self._create_var_c()
        self.var_p = self._create_var_p()
        self.hdp.close()

    def _get_txs(self):
        """

        :return:
        """
        return hgvs.easy.am37.relevant_transcripts(self.var_g)

    def _create_var_g(self):
        """

        :return:
        """
        return self.hp.parse_hgvs_variant(self.hgvs_g)

    def _create_var_c(self):
        """

        :return:
        """
        var_c = []
        # hdp = hgvs.dataproviders.uta.connect()
        # vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method='splign')
        for tx in self.tx_list:
            try:
                var_c.append(self.vm.g_to_c(self.var_g, tx))
            except:
                var_c.append(None)
        # hdp.close()
        return var_c

    def _create_var_p(self):
        """

        :return:
        """
        var_p = []
        # hdp = hgvs.dataproviders.uta.connect()
        # vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method='splign')
        for coding in self.var_c:
            if coding:
                var_p.append(self.vm.c_to_p(coding))
            else:
                var_p.append(None)
        # hdp.close()
        return var_p

    def _prep_hgvs(self):
        """

        :return:
        """
        pos_part = self.correct_indel_coords(self.pos, self.ref, self.alt[0])
        new_hgvs = self.chrom + ':g.' + pos_part
        return new_hgvs

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
    def __init__(self, filename):
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
        self.vcf_reader.close()

    def _create_rec(self):
        """
        """
        vrnts = {}
        for rec in self.vcf_reader:
            chrom = CHROM_MAP[rec.CHROM]
            pos = rec.POS
            ref = rec.REF
            alt = [x.value for x in rec.ALT][0]
            uniq_key = (rec.CHROM, pos, ref, alt)
            pos_part = correct_indel_coords(pos, ref, alt)
            new_hgvs = chrom + ':g.' + pos_part
            if uniq_key not in vrnts:
                vrnts[uniq_key] = new_hgvs
            else:
                raise Exception("WHY!!!")
        return vrnts


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

        # i = 0
        # results = []
        # for vrnt in self.vcf_reader:
        #     if i % 100 == 0:
        #         print(i)
        #     this_vcf = VcfRec(vrnt, self.hp, self.hdp, self.vm)
        #     results.append(this_vcf)
        #     i += 1
        # return results
    #     return VcfRec(vrnt)
    #
    #
    # def get_txs(self):
    #     """
    #
    #     :return:
    #     """
    #     pool = ProcessingPool(nodes=24)
    #     results = pool.map(self._create_rec, [x for x in self.vcf_reader])
    #     return results
#
#
# def access_region_tx(uniq_key, vrnt):
#     """
#
#     :return:
#     """
#     result = {}
#     hdp = hgvs.dataproviders.uta.connect()
#     for entry in vrnt:
#         this_tx = hdp.get_tx_for_region(str(entry.ac), 'splign', str(entry.posedit.pos.start),
#                                         str(entry.posedit.pos.end))
#         if entry not in result:
#             result[entry] = this_tx
#         else:
#             result[entry].append(this_tx)
#         # if uniq_key not in result:
#         #     result[uniq_key] = this_tx
#         # else:
#         #     result[uniq_key].append(this_tx)
#     hdp.close()
#     return result

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


def hgvs_build(vrnt):
    """

    :return:
    """
    hdp = hgvs.dataproviders.uta.connect()
    vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method='splign')
    hp = hgvs.parser.Parser()
    var_g = hp.parse_hgvs_variant(vrnt)
    tx_list = hgvs.easy.am37.relevant_transcripts(var_g)

    var_c = []
    for tx in tx_list:
        try:
            var_c.append(vm.g_to_c(var_g, tx))
        except:
            pass

    var_p = []
    for coding in var_c:
        if coding:
            var_p.append(vm.c_to_p(coding))
        else:
            pass

    hdp.close()
    return str(var_g), ','.join([str(x) for x in var_c]), ','.join([str(x) for x in var_p])


class VcfWriter:
    """

    """
    def __init__(self, filename, vcf_reader):
        self.vcf_writer = vcfpy.Writer.from_path(filename, vcf_reader.header)


def main():

    args = supply_args()

    # Initialize the HGVS package objects.
    # hdp = hgvs.dataproviders.uta.connect()
    # vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method='splign')
    hp = hgvs.parser.Parser()

    # {('NC_000001.10', 11169379, 'T', 'C'): [['NM_004958.4', 'NC_000001.10', -1, 'splign', 11166591, 11322608], ['NM_004958.3', 'NC_000001.10', -1, 'splign', 11166587, 11322608]]},

    hgvs_blob = {}

    if args.input_vcf:
        # vcf_writer = vcfpy.Writer.from_path(args.output_vcf, vcf_reader.header)
        my_vcf = VcfReader(args.input_vcf).my_vrnts

    pool = ProcessingPool(nodes=24)
    results = pool.map(hgvs_build, my_vcf.values())

    # print(results)

    # for k, v in my_vcf.items():
    #
    #     if i % 100 == 0:
    #         print(i)
    #
    # hgvs_blob[k] = {'HGVS_G':str(var_g),
    #                 'HGVS_C': ','.join([str(x) for x in var_c]),
    #                 'HGVS_P': ','.join([str(x) for x in var_p])}
    #
    # hdp = hgvs.dataproviders.uta.connect()
    # vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method='splign')
    # hp = hgvs.parser.Parser()
    #
    # i = 0
    # for k, v in my_vcf.items():
    #
    #     if i % 100 == 0:
    #         print(i)
    #     var_g = hp.parse_hgvs_variant(v)
    #     tx_list = hgvs.easy.am37.relevant_transcripts(var_g)
    #     var_c = []
    #     for tx in tx_list:
    #         try:
    #             var_c.append(vm.g_to_c(var_g, tx))
    #         except:
    #             var_c.append(None)
    #     var_p = []
    #     for coding in var_c:
    #         if coding:
    #             var_p.append(vm.c_to_p(coding))
    #         else:
    #             var_p.append(None)
    #
    #     hgvs_blob[k] = {'HGVS_G': str(var_g),
    #                     'HGVS_C': ','.join([str(x) for x in var_c]),
    #                     'HGVS_P': ','.join([str(x) for x in var_p])}
    #     i += 1
    #
    # print(hgvs_blob)
    #
    # hdp.close()




    # print(hgvs_blob)

        # tx_blobs = get_txs(my_vcf.my_varg)

    # for blob in tx_blobs:
    #     for tx_list in blob.values():
    #         for tx in tx_list:
    #             hdp = hgvs.dataproviders.uta.connect()
    #             vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method='splign')
    #             var_c = vm.g_to_c(tx, str(tx[0]))
    #             print(var_c)


    #
    # for entry in my_vcf:
    #     print(str(my_vcf.var_g))

        # new_var_g = []
        # for i in range(5000):
        #     new_var_g.append(my_vcf.var_g[0])

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
