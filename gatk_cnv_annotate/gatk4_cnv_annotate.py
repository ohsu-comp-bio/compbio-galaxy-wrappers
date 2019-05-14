import argparse
import gffutils

VERSION = '0.1.0'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='Galaxy BioBlend workflow generation script.')
    parser.add_argument('--gffutils_db', default='gffutils_db', help='GffUtils DB')
    parser.add_argument('--cnv_infile', help='CNV TSV Infile')
    parser.add_argument('--genes', help='Gene list from which to produce output.')
    parser.add_argument('--outfile', help='CNV enhanced output.')
    parser.add_argument('--outfile_genes', help='CNV enhanced output, gene-based.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    args = parser.parse_args()
    return args

def coords_to_codons(start, end):
    total_bp = end - start + 1
    return ((total_bp) // 3, (total_bp) % 3)

def genomic_to_codon(cds, coord):
    for entry in cds:
        if coord < entry[0]:
            pass
        elif coord >= entry[0] and coord <= entry[1]:
            pass
        else:
            pass

def codon_map(my_cds):
    """

    :param my_cds:
    :return:
    """
    codon_map = {}
    j = 1
    k = 1
    for entry in my_cds:
        for i in range(entry[0], entry[1]+1):
            codon_map[i] = j
            if k == 3:
                k = 1
                j += 1
            else:
                k += 1
    return codon_map

class GffUtilUtil(object):
    """
    Take an interval, supply info from gffutil db.
    """
    def __init__(self, dbname, chrom, start, end):
        self.db = gffutils.FeatureDB(dbname)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.my_region = (chrom, start, end)
        self.genes = self._get_genes()
        self.start_exon, self.start_gene = self._find_start_exon()
        self.end_exon, self.end_gene = self._find_end_exon()
        self.gffutil_rec = self._coll_gffutil_metrics()

    def _coll_gffutil_metrics(self):
        this_rec = {'START_GENE': self.start_gene,
                    'END_GENE': self.end_gene,
                    'START_EXON': self.start_exon,
                    'END_EXON': self.end_exon,
                    'GENES': self.genes}
        return this_rec

    def _find_start_exon(self):
        # Find exon number of start coord
        for entry in self.db.region(self.my_region, featuretype=['transcript', 'mRNA', 'primary_transcript']):
            if entry.strand == '+':
                for exon in self.db.children(entry.id, featuretype='exon'):
                    if self.start < exon.end:
                        exon_no = list(self.db.children(entry.id, featuretype='exon')).index(exon) + 1
                        return exon_no, entry.attributes['gene'][0]
            elif entry.strand == '-':
                for exon in list(self.db.children(entry.id, featuretype='exon'))[::-1]:
                    if self.start < exon.end:
                        exon_cnt = len(list(self.db.children(entry.id, featuretype='exon')))
                        exon_no = exon_cnt - (list(self.db.children(entry.id, featuretype='exon'))[::-1].index(exon))
                        return exon_no, entry.attributes['gene'][0]
            else:
                print(self.my_region)
                raise Exception('No strand information given, please check.')

        return 'intergenic', self.genes.split(',')[0]


    def _find_end_exon(self):
        # Find exon number of end coord>
        for entry in self.db.region(self.my_region, featuretype=['transcript', 'mRNA', 'primary_transcript']):
            if entry.strand == '+':
                for exon in self.db.children(entry.id, featuretype='exon'):
                    if self.end < exon.start:
                        exon_no = list(self.db.children(entry.id, featuretype='exon')).index(exon)
                        return exon_no, entry.attributes['gene'][0]
                    elif self.end < exon.end:
                        exon_no = list(self.db.children(entry.id, featuretype='exon')).index(exon) + 1
                        return exon_no, entry.attributes['gene'][0]
            elif entry.strand == '-':
                for exon in list(self.db.children(entry.id, featuretype='exon'))[::-1]:
                    if self.end < exon.start:
                        exon_cnt = len(list(self.db.children(entry.id, featuretype='exon')))
                        exon_no = exon_cnt - (list(self.db.children(entry.id, featuretype='exon'))[::-1].index(exon)) + 1
                        return exon_no, entry.attributes['gene'][0]
                    elif self.end < exon.end:
                        exon_cnt = len(list(self.db.children(entry.id, featuretype='exon')))
                        exon_no = exon_cnt - (list(self.db.children(entry.id, featuretype='exon'))[::-1].index(exon))
                        return exon_no, entry.attributes['gene'][0]
            else:
                print(self.my_region)
                raise Exception('No strand information given, please check.')

        return 'intergenic', self.genes.split(',')[-1]

    def _get_genes(self):
        """

        :return:
        """
        my_genes = []
        for entry in self.db.region(self.my_region, featuretype='gene'):
            my_genes.extend(entry.attributes['gene'])
        return ','.join(my_genes)


class PicardInterval(object):
    def __init__(self, line, titles):
        self.line = line
        self.titles = titles
        self.ival = self._make_ival()

    def _make_ival(self):
        this_ival = {}
        for i in range(len(self.titles)):
            this_ival[self.titles[i]] = self.line[i]
        return this_ival

    def add_gffutil_field(self, dbname, common_to_refseq):
        chrom = common_to_refseq[self.ival['CONTIG']]
        start = int(self.ival['START'])
        end = int(self.ival['END'])
        self.ival.update(GffUtilUtil(dbname, chrom, start, end).gffutil_rec)


class PicardIntervals(object):
    def __init__(self, filename):
        self.myfile = open(filename, 'r')
        self._break_it()
        self.myfile.close()

    def _break_it(self):
        self.header = []
        self.regions = []
        for line in self.myfile:
            if line.startswith('@'):
                self.header.append(line)
            elif line.startswith('CONTIG'):
                self.titles = line.rstrip('\n').split('\t')
            else:
                ival = PicardInterval(line.rstrip('\n').split('\t'), self.titles)
                self.regions.append(ival)


def get_common_chroms(dbname):
    common_to_refseq = {}
    db = gffutils.FeatureDB(dbname)
    for entry in db.features_of_type('region'):
        if 'Name' in entry.attributes:
            if entry.chrom.startswith('NC'):
                common_to_refseq[entry.attributes['Name'][0]] = entry.chrom
        else:
            common_to_refseq[entry.chrom] = None
    return common_to_refseq


class CnvWriter(object):

    def write_me(self, filename, intervals):
        write_order = ['CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO', 'CALL',
                        'START_GENE', 'START_EXON', 'END_GENE', 'END_EXON', 'GENES']
        myfile = open(filename, 'w')
        for entry in intervals.header:
            myfile.write(entry)
        myfile.write('\t'.join(write_order))
        myfile.write('\n')
        for entry in intervals.regions:
            to_write = []
            for title in write_order:
                if title in entry.ival:
                    to_write.append(str(entry.ival[title]))
                else:
                    to_write.append('')
            myfile.write('\t'.join(to_write))
            myfile.write('\n')
        myfile.close()

    def gene_output_write(self, filename, my_genes, my_ints):
        gene_output = []
        handle_out = open(filename, 'w')
        header = ['GENE', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO', 'CALL', 'START_EXON', 'END_EXON']
        handle_out.write('\t'.join(header))
        handle_out.write('\n')
        for entry in my_ints:
            for gene in my_genes:
                if gene in entry.ival['GENES'].split(','):
                    to_write = [gene, str(entry.ival['NUM_POINTS_COPY_RATIO']), str(entry.ival['MEAN_LOG2_COPY_RATIO']),
                                entry.ival['CALL']]
                    if gene == entry.ival['START_GENE']:
                        to_write.append(str(entry.ival['START_EXON']))
                    else:
                        to_write.append('NA')
                    if gene == entry.ival['END_GENE']:
                        to_write.append(str(entry.ival['END_EXON']))
                    else:
                        to_write.append('NA')
                    gene_output.append(to_write)

        for entry in sorted(gene_output):
            handle_out.write('\t'.join(entry))
            handle_out.write('\n')

        handle_out.close()

    def write_amp_codon(self):
        pass


class GeneInfo(object):
    """
    From a list of genes (HGNC), find coordinates, any additional information we want that can be obtained from
    gffutils.
    """
    def __init__(self, filename, dbname, common_to_refseq):
        self.genes = self._file_to_list(filename)
        self.db = gffutils.FeatureDB(dbname)
        self.common_to_refseq = common_to_refseq
        self.gene_info = self._get_gene_info()

    def _file_to_list(self, infile):
        genes = []
        with open(infile, 'rU') as myfile:
            for line in myfile:
                if len(line.split('\t')) > 1:
                    raise Exception("Line has more than one column.")
                line = line.rstrip('\n')
                genes.append(line)
        return genes

    def _get_gene_info(self):
        gene_info = {}
        for entry in self.db.features_of_type('gene', order_by='start'):
            if entry[0] in self.common_to_refseq.values():
                gene = entry['Name'][0]
                if gene in self.genes:
                    chrom = entry[0]
                    start = entry[3]
                    stop = entry[4]
                    gene_info[gene] = {'GENE_CHROM': chrom,
                                       'GENE_START': start,
                                       'GENE_STOP': stop}
        return gene_info

def main():

    args = supply_args()
    common_to_refseq = get_common_chroms(args.gffutils_db)
    my_genes = GeneInfo(args.genes, args.gffutils_db, common_to_refseq).genes

    pic_ints = PicardIntervals(args.cnv_infile)
    my_ints = pic_ints.regions

    for entry in my_ints:
        # if entry.ival['CALL'] != '0':
        entry.add_gffutil_field(args.gffutils_db, common_to_refseq)

    CnvWriter().write_me(args.outfile, pic_ints)
    CnvWriter().gene_output_write(args.outfile_genes, my_genes, my_ints)



main()