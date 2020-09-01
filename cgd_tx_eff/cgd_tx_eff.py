#!/usr/bin/env python

import argparse
import hgvs
import hgvs.parser
import vcfpy

VERSION = '0.1.0'
# Going to work on a better way to do this, but we definitely can't spin this parser up over and over again.
hp = hgvs.parser.Parser()


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
    G	C	1	240.78	8	1	20992819	.	G	C	240.78
	.	AC=2;AF=1.00;AN=2;DP=8;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;QD=30.10;SOR=0.693	GT:AD:DP:GQ:PL
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
        self.vtype = []
        self.veff = []

        self.chrom = rec[2]
        self.pos = rec[3]
        self.ref = rec[13]
        self.alt = rec[14]

        self.uniq_key = (self.chrom, self.pos, self.ref, self.alt)
        self.metrics = {self.uniq_key: {}}


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
        self.vtype.append(self.rec[0])
        self.info = self.rec[1].split(',')
        for metric in self.info:
            new_met = self._make_metric_dict(metric)
            if new_met:
                self.metrics[self.uniq_key].update(new_met)

    def _make_metric_dict(self, metrics):
        """
        Pull a metric dict out of the info string based on index.
        :return:
        """
        hgnc = []
        exon = []
        hgvs_c = []
        hgvs_p = []

        slen = len(metrics.split(':'))
        if slen == 1:
            if "(dist=" not in metrics:
                tx = metrics
            else:
                tx = None
        elif slen == 2:
            tx = metrics.split('(')[0]
        elif slen == 3:
            tx = metrics.split('(')[0]
            exon.append(metrics.split(':')[1][4:])
        else:
            raise UserWarning("There are an incorrect number of fields listed in: " + str(self.info))

        if tx:
            return {tx: {
                'AAP': [],
                'BASEP': [],
                'EXON': exon,
                'HGNC': hgnc,
                'HGVSC1': hgvs_c,
                'HGVSP1': hgvs_p,
                'HGVSP3': [],
                'SOURCE': self.tfx_type,
                'TXP': [],
                'VFX': self.veff,
                'PVT': self.vtype
            }}

        return None


class AnnovarRecExonVrntFunc(AnnovarRec):
    """
    need to pull from 1st field:
    DVL1:NM_001330311.2:exon14:c.1687G>A:p.G563S,DVL1:NM_004421.3:exon14:c.1612G>A:p.G538S,
    """
    def __init__(self, rec):
        self.rec = rec.rstrip('\n').split('\t')[1:]
        super().__init__(self.rec)
        self.veff.append(self.rec[0])
        self.info = self.rec[1].rstrip(',').split(',')
        if self.rec[1] != 'UNKNOWN':
            for metric in self.info:
                self.metrics[self.uniq_key].update(self._make_metric_dict(metric))

    def _make_metric_dict(self, metrics):
        """
        Pull a metric dict out of the info string based on index.
        :return:
        """
        hgnc = []
        exon = []
        hgvs_c = []
        hgvs_p = []
        hgvs_three = []
        hgvs_basep = []
        hgvs_aap = []

        tx = metrics.split(':')[1]
        hgnc.append(metrics.split(':')[0])
        raw_exon = metrics.split(':')[2]
        if raw_exon != 'wholegene':
            exon.append(metrics.split(':')[2][4:])

            this_c = metrics.split(':')[3]
            full_c = ':'.join([tx, this_c])
            hgvs_c.append(this_c)
            hgvs_parser = Hgvs(full_c)
            if hgvs_parser.start:
                hgvs_basep.append(hgvs_parser.start)

            this_p = metrics.split(':')[4]
            full_p = ':'.join([tx, this_p])
            hgvs_p.append(this_p)
            hgvs_parser = Hgvs(full_p)

            if hgvs_parser.threep:
                hgvs_three.append(hgvs_parser.threep)
            if hgvs_parser.start:
                hgvs_aap.append(hgvs_parser.start[3:])

        return {tx: {
                'AAP': hgvs_aap,
                'BASEP': hgvs_basep,
                'EXON': exon,
                'HGNC': hgnc,
                'HGVSC1': hgvs_c,
                'HGVSP1': hgvs_p,
                'HGVSP3': hgvs_three,
                'SOURCE': self.tfx_type,
                'TXP': [],
                'VFX': self.veff,
                'PVT': self.vtype
            }}


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

    def metrics_push(self, annovar):
        """
        Just an awful function, so bad.  The idea here is to combine all of the metrics from some number of annovar
        inputs in to one metrics blob.
        :param annovar:
        :return:
        """
        for coord, metric in annovar.items():
            for met in metric:
                if coord not in self.metrics:
                    self.metrics.update(met.metrics)
                else:
                    for tx in met.metrics.values():
                        for tx_name, fields in tx.items():
                            if tx_name not in self.metrics[coord]:
                                self.metrics[coord].update(tx)
                            else:
                                for field, val in fields.items():
                                    if val:
                                        if not self.metrics[coord][tx_name][field]:
                                            self.metrics[coord][tx_name][field] = val
                                        else:
                                            for v in val:
                                                if v not in self.metrics[coord][tx_name][field]:
                                                    self.metrics[coord][tx_name][field].append(v)


class VcfReader:
    """

    """
    def __init__(self, filename):
        self.vcf_reader = vcfpy.Reader.from_path(filename)
        self.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_HGVSBC_G'),
                                                           ('Number', '.'),
                                                           ('Type', 'String'),
                                                           ('Description',
                                                            'HGVS genomic reference, as produced by pyhgvs.')]))
        self.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_HGVSBC_C'),
                                                           ('Number', '.'),
                                                           ('Type', 'String'),
                                                           ('Description',
                                                            'HGVS coding reference, as produced by pyhgvs.')]))
        self.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_HGVSBC_P1'),
                                                           ('Number', '.'),
                                                           ('Type', 'String'),
                                                           ('Description',
                                                            'HGVS 1-letter protein reference, as produced '
                                                            'by pyhgvs.')]))
        self.vcf_reader.header.add_info_line(vcfpy.OrderedDict([('ID', 'TFX_HGVSBC_P3'),
                                                           ('Number', '.'),
                                                           ('Type', 'String'),
                                                           ('Description',
                                                            'HGVS 3-letter protein reference, as produced '
                                                            'by pyhgvs.')]))


    def get_vrnts(self):
        vrnts = {}
        for entry in self.vcf_reader:
            vrnt = VcfRec(entry)
            if vrnt.uniq_key not in vrnts:
                vrnts[vrnt.uniq_key] = vrnt
            else:
                raise Exception("sffsa")
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

    Func is

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
            return None


def main():
    args = supply_args()
    myvcf = VcfReader(args.infile).get_vrnts()
    coll = CollectMetrics()

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


    key = ('X', '153035798', 'G', 'A')
    print(coll.metrics[key])
    print(myvcf[key].rec)

    for entry in coll.metrics:
        if entry in myvcf:

            print(coll.metrics[entry])
            print(myvcf[entry].rec)
            exit(1)

    # for key, val in coll.metrics.items():
    #     for tx in val.values():
    #         if len(tx['EXON']) > 1:
    #             print(key)
    #             print(val)



if __name__ == "__main__":
    main()
