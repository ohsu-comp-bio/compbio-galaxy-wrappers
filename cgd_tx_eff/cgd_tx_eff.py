#!/usr/bin/env python

from collections import defaultdict
from itertools import chain
from operator import methodcaller

import argparse
import hgvs
import hgvs.dataproviders.uta
import hgvs.parser
import hgvs.assemblymapper
import hgvs.exceptions
import vcfpy

VERSION = '0.1.1'
# Going to work on a better way to do this, but we definitely can't spin this parser up over and over again.
hp = hgvs.parser.Parser()
# hdp = hgvs.dataproviders.uta.connect()
# vm = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign')


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('infile', help='Input VCF to apply Annovar annotations to.')
    parser.add_argument('outfile', help='Output VCF')
    parser.add_argument('--evf', help='Input exonic_variant_function Annovar file.')
    parser.add_argument('--vf', help='Input variant_function Annovar file.')
    parser.add_argument('--ccds_evf', help='Input CCDS exonic_variant_function Annovar file.')
    parser.add_argument('--ccds_vf', help='Input CCDS variant_function Annovar file.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    if not args.evf and not args.vf:
        raise SyntaxError("Must specify either vf or evf, or both.")

    return args


class AnnovarRec:
    """
    Two files contain very similar fields, with exception of first one where it denotes line number.

    splicing	NM_001122819.3,NM_001287212.2,NM_020816.4	1	20992819	20992819
    G	C	1	240.78	8	1	20992819	.	G	C	240.78.
    AC=2;AF=1.00;AN=2;DP=8;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=30.10;SOR=0.693	GT:AD:DP:GQ:PL
    1/1:0,8:8:24:269,24,0

    line96	synonymous SNV	CDK11A:NM_001313982.2:exon4:c.321A>G:p.E107E,CDK11B:
    NM_001787.3:exon4:c.321A>G:p.E107E,CDK11A:NM_024011.4:exon4:c.321A>G:p.E107E,CDK11A:
    NM_001313896.2:exon4:c.321A>G:p.E107E,CDK11B:NM_033489.3:exon5:c.219A>G:p.E73E,CDK11A:
    NM_033529.4:exon4:c.321A>G:p.E107E,CDK11B:NM_033486.3:exon4:c.321A>G:p.E107E,CDK11B:
    NM_001291345.2:exon4:c.321A>G:p.E107E,	1	1650801	1650801	T	C	0.5	176.77	19	1	1650801	.
    T	C	176.77	.	AC=1;AF=0.500;AN=2;BaseQRankSum=0.684;ClippingRankSum=0.000;DP=19;ExcessHet=3.0103;
    FS=1.808;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=9.30;ReadPosRankSum=-0.730;SOR=1.140
    GT:AD:DP:GQ:PL	0/1:13,6:19:99:205,0,602
    """

    def __init__(self, rec, tfx_type="ANNOVAR"):
        self.tfx_type = tfx_type
        self.rec = rec

        self.chrom = self.rec[10]
        self.pos = self.rec[11]
        self.ref = self.rec[13]
        self.alt = self.rec[14]
        self.uniq_key = (self.chrom, self.pos, self.ref, self.alt)

        self.metrics = {}


class AnnovarRecVrntFunc(AnnovarRec):
    """
    Different possibilities in field 1:

    Contains c.
    NM_001291367.2(NM_001291367.2:c.-3T>C),NM_001369897.1(NM_001369897.1:c.-140T>C),
    NM_001369898.1(NM_001369898.1:c.-140T>C)

    Might contain exon numbers:
    NM_001305275.2(NM_001305275.2:exon34:c.5662+5C>T),NM_001364727.2(NM_001364727.2:exon33:c.5347+5C>T),
    NM_198576.4(NM_198576.4:exon33:c.5651+5C>T)

    Contains distances, since this is intergenic:
    NR_026818.1(dist=29825),NM_001005484.1(dist=3185)

    Contains only tx names, likely intronic
    NR_047519.1,NR_047521.1,NR_047523.1,NR_047524.1,NR_047525.1

    Depending on tx, might be intergenic (dist not in all entries):
    NM_001256456.1,NM_001256460.1,NM_001256462.1,NM_001256463.2,NM_017871.6(dist=961)
    """
    def __init__(self, rec):
        self.rec = rec.rstrip('\n').split('\t')
        super().__init__(self.rec)

        # vtype only in this file, and is in the 1st column.
        # each field will have either a splicing or vtype designation, but not both
        if self.rec[0] == 'splicing' or self.rec[0] == 'ncRNA_splicing':
            self.vtype = None
            self.splicing = self.rec[0]
        else:
            self.splicing = None
            self.vtype = self.rec[0]

        self.info = self.rec[1].split(',')
        for metric in self.info:
            new_met = self._make_metric_dict(metric)
            if new_met['TXC']:
                tx = new_met['TXC']
                self.metrics[tx] = new_met

    def _make_metric_dict(self, metrics):
        """
        Pull a metric dict out of the info string based on index.
        :return:
        """
        exon = None
        hgvs_c = None

        slen = len(metrics.split(':'))
        if slen == 1:
            if "(dist=" not in metrics:
                tx = metrics
            else:
                tx = None
        elif slen == 2:
            tx = metrics.split('(')[0]
            hgvs_c = metrics.split(':')[1].rstrip(')')
        elif slen == 3:
            tx = metrics.split('(')[0]
            exon = metrics.split(':')[1][4:]
            hgvs_c = metrics.split(':')[2].rstrip(')')
        else:
            raise UserWarning("There are an incorrect number of fields listed in: " + str(self.info))

        if hgvs_c:
            full_c = ':'.join([tx, hgvs_c])
            hgvs_parser = Hgvs(full_c)
            hgvs_basep = hgvs_parser.start
        else:
            hgvs_basep = None

        return {
            'AAP': None,
            'BASEP': hgvs_basep,
            'EXON': exon,
            'HGNC': None,
            'HGVSC': hgvs_c,
            'HGVSP1': None,
            'HGVSP3': None,
            'SOURCE': self.tfx_type,
            'SPLICE': self.splicing,
            'TXC': tx,
            'PVT': self.vtype,
            'VFX': None
        }


class AnnovarRecExonVrntFunc(AnnovarRec):
    """
    need to pull from 1st field:
    DVL1:NM_001330311.2:exon14:c.1687G>A:p.G563S,DVL1:NM_004421.3:exon14:c.1612G>A:p.G538S,
    """
    def __init__(self, rec):
        self.rec = rec.rstrip('\n').split('\t')[1:]
        super().__init__(self.rec)

        self.veff = self.rec[0]
        self.info = self.rec[1].rstrip(',').split(',')
        if self.rec[1] != 'UNKNOWN':
            for metric in self.info:
                new_met = self._make_metric_dict(metric)
                if new_met['TXC']:
                    tx = new_met['TXC']
                    self.metrics[tx] = new_met

    def _make_metric_dict(self, metrics):
        """
        Pull a metric dict out of the info string based on index.
        :return:
        """
        exon = None
        hgvs_c = None
        hgvs_p = None
        hgvs_three = None
        hgvs_basep = None
        hgvs_aap = None

        tx = metrics.split(':')[1]
        hgnc = metrics.split(':')[0]
        raw_exon = metrics.split(':')[2]
        if raw_exon != 'wholegene':
            exon = metrics.split(':')[2][4:]
            hgvs_c = metrics.split(':')[3]
            full_c = ':'.join([tx, hgvs_c])
            hgvs_parser = Hgvs(full_c)

            hgvs_basep = hgvs_parser.start
            try:
                hgvs_p = metrics.split(':')[4]
                full_p = ':'.join([tx, hgvs_p])
                hgvs_parser = Hgvs(full_p)
                hgvs_three = hgvs_parser.threep
                hgvs_aap = hgvs_parser.start[3:]
            except IndexError:
                hgvs_p = None
                full_p = None
                hgvs_three = None
                hgvs_aap = None

        return {
                'AAP': hgvs_aap,
                'BASEP': hgvs_basep,
                'EXON': exon,
                'HGNC': hgnc,
                'HGVSC': hgvs_c,
                'HGVSP1': hgvs_p,
                'HGVSP3': hgvs_three,
                'SOURCE': self.tfx_type,
                'SPLICE': None,
                'TXC': tx,
                'VFX': self.veff,
                'PVT': None
            }


class AnnovarReader:
    def __init__(self, filename):
        self.filename = filename


class AnnovarVrntFunc(AnnovarReader):
    def __init__(self, filename):
        super().__init__(filename)

    def read_annovar(self):
        anno = {}
        with open(self.filename, 'r') as myfile:
            for rec in myfile:
                this_rec = AnnovarRecVrntFunc(rec)
                if this_rec.uniq_key not in anno:
                    anno[this_rec.uniq_key] = [this_rec]
                else:
                    anno[this_rec.uniq_key].append(this_rec)

        return anno


class AnnovarExonVrntFunc(AnnovarReader):
    def __init__(self, filename):
        super().__init__(filename)

    def read_annovar(self):
        """
        Produces structure that looks like:
        (chrom, coord, ref, alt): {tx1: {field: val, ...}, tx2: {field: val, ...}}
        :return:
        """
        anno = {}
        with open(self.filename, 'r') as myfile:
            for rec in myfile:
                this_rec = AnnovarRecExonVrntFunc(rec)
                if this_rec.uniq_key not in anno:
                    anno[this_rec.uniq_key] = [this_rec]
                else:
                    anno[this_rec.uniq_key].append(this_rec)
            return anno


class CollectMetrics:
    """
    Get the metrics from evf and vf, then combine them in to one structure.
    """
    def __init__(self):
        self.metrics = {}

    @staticmethod
    def _metrics_update(orig, new):
        """
        When we are getting metrics from multiple sources, figure out what the combined on should look like.
        :return:
        """
        revsd = orig
        for k, v in orig.items():
            if not v:
                revsd[k] = new[k]
            elif new[k]:
                if new[k] != v:
                    # LOG ME, do something
                    print(orig)
                    print(new)
            elif not new[k] or v:
                pass
            else:
                raise Exception("_metrics_update error")
        return revsd

    def metrics_push(self, annovar):
        """
        We want metrics to looks like:
        {coord: {tx: {field: val, etc.}}}
        :param annovar:
        :return:
        """
        for coord, metric in annovar.items():
            for met in metric:
                if coord not in self.metrics:
                    self.metrics[coord] = met.metrics
                else:
                    for tx, fields in met.metrics.items():
                        if tx not in self.metrics[coord]:
                            self.metrics[coord][tx] = fields
                        # If we have already seen this transcript at this coordinate, will need to perform an update.
                        else:
                            new_met = self._metrics_update(self.metrics[coord][tx], fields)
                            self.metrics[coord][tx] = new_met


class VcfReader:
    """
    Handles reading the VCF and placing variants in to a dictionary.
    """
    def __init__(self, filename):
        self.vcf_reader = vcfpy.Reader.from_path(filename)

    def get_vrnts(self):
        vrnts = {}
        for entry in self.vcf_reader:
            vrnt = VcfRec(entry)
            if vrnt.uniq_key not in vrnts:
                vrnts[vrnt.uniq_key] = vrnt
            else:
                # LOG ME
                raise Exception("Duplicate variant in VcfReader.get_vrnts.")
        self.vcf_reader.close()
        return vrnts


class VcfRec:
    """
    Need fields from INFO:

    Func.refGeneWithVer=exonic;
    Gene.refGeneWithVer=NOC2L;
    GeneDetail.refGeneWithVer=.;
    ExonicFunc.refGeneWithVer=synonymous_SNV;
    AAChange.refGeneWithVer=NOC2L:NM_015658.4:exon9:c.A918G:p.E306E;

    Func.ccdsGene=exonic;
    Gene.ccdsGene=CCDS3.1;
    GeneDetail.ccdsGene=.;
    ExonicFunc.ccdsGene=synonymous_SNV;
    AAChange.ccdsGene=CCDS3.1:CCDS3.1:exon9:c.A918G:p.E306E

    """
    def __init__(self, rec):
        self.rec = rec
        self.chrom = rec.CHROM
        self.coord = str(rec.POS)
        self.ref = rec.REF
        self.alt = rec.ALT[0].serialize()
        self.uniq_key = (self.chrom, self.coord, self.ref, self.alt)


class Hgvs:
    """
    Keep HGVS-package-related stuff here.  For the parser to function accordingly, you need to provide something
    as a sequence identifier, even if it is gibberish.  Such as:
    "GIBBERISH:c.1792G>A"
    We aren't mapping between c. and p., so it doesn't matter.
    """
    def __init__(self, term):
        self.term = term
        self.hp = hp
        self.start = str(self.hgvs_start())
        self.threep = self.hgvs_three()

    def hgvs_start(self):
        """
        Take a notation, and use HGVS package to provide start position.  For c., can take as is.  For p., will
        need to later strip off AA start to isolate codon number.
        :return:
        """
        try:
            return self.hp.parse(self.term).posedit.pos.start
        except hgvs.exceptions.HGVSParseError:
            # Log me
            # print(self.term)
            return None

    def hgvs_three(self):
        """
        This package automatically provides back the three character representations of amino acid names, so we can
        simply return it after stringifying.
        :return:
        """
        try:
            return 'p.' + str(self.hp.parse(self.term).posedit)
        except hgvs.exceptions.HGVSParseError:
            # Log me
            # print(self.term)
            return None


class VcfWriter:
    def __init__(self, filename, reader):
        self.reader = reader
        self.writer = vcfpy.Writer.from_path(filename, reader.header)

    def write_metrics(self, vrnts, metrics, fld_prefix='TFX'):
        for coord, vrnt in vrnts.items():
            if coord in metrics:
                for tx, mets in metrics[coord].items():
                    tx_key = '_'.join([fld_prefix, tx])
                    if tx_key not in vrnt.rec.INFO:
                        vrnt.rec.INFO[tx_key] = mets
                a = len(set([len(v) for k, v in vrnt.rec.INFO.items() if k.startswith('TFX_')]))
                if a > 1:
                    print(vrnt.rec.INFO)
                    raise Exception("TFX field lengths are not the same, please check.")
                self.writer.write_record(vrnt.rec)
            else:
                # LOG ME
                print("NOT IN METRICS: " + str(coord))
        self.writer.close()


def main():
    args = supply_args()
    myvcf = VcfReader(args.infile)
    vrnts = myvcf.get_vrnts()
    coll = CollectMetrics()

    # Set VCF header.
    # 'AAP': hgvs_aap,
    # 'BASEP': hgvs_basep,
    # 'EXON': exon,
    # 'HGNC': hgnc,
    # 'HGVSC': hgvs_c,
    # 'HGVSP1': hgvs_p,
    # 'HGVSP3': hgvs_three,
    # 'SOURCE': self.tfx_type,
    # 'SPLICE': None,
    # 'TXC': tx,
    # 'VFX': self.veff,
    # 'PVT': None
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_AAP'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'Amino acid start position.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_BASEP'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'Coding sequence start position.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_EXON'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'Exon number associated with given '
                                                                            'transcript.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_HGNC'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'HGNC gene symbol.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_HGVSC'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'HGVS cdot nomenclature.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_HGVSP1'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'HGVS pdot nomenclature, single letter '
                                                                            'amino acids.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_HGVSP3'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'HGVS pdot nomenclature, three letter '
                                                                            'amino acids.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_SOURCE'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'Annotation source.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_SPLICE'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'Splice site annotation.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_TXC'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'Transcript identifier.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_VFX'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'Variant effect annotation.')]))
    myvcf.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_PVT'),
                                                            ('Number', '.'),
                                                            ('Type', 'String'),
                                                            ('Description', 'Variant type or location annotation.')]))

    if args.evf:
        my_evf = AnnovarExonVrntFunc(args.evf).read_annovar()
        coll.metrics_push(my_evf)
    if args.vf:
        my_vf = AnnovarVrntFunc(args.vf).read_annovar()
        coll.metrics_push(my_vf)
    if args.ccds_evf:
        my_ccds_evf = AnnovarExonVrntFunc(args.ccds_evf).read_annovar()
        coll.metrics_push(my_ccds_evf)
    if args.ccds_vf:
        my_ccds_vf = AnnovarVrntFunc(args.ccds_vf).read_annovar()
        coll.metrics_push(my_ccds_vf)

    to_write = {}
    for coord, txs in coll.metrics.items():
        comb = defaultdict(list)
        dict_items = map(methodcaller('items'), (txs.values()))
        for k, v in chain.from_iterable(dict_items):
            comb[k].append(v)
        to_write[coord] = comb
        # Look for any instances where the number of TFX fields do not match.
        num_vals = len(set([len(x) for x in comb.values()]))
        assert num_vals == 1 or not num_vals, "%s" % comb
    writer = VcfWriter(args.outfile, myvcf.vcf_reader)
    writer.write_metrics(vrnts, to_write)


if __name__ == "__main__":
    main()
