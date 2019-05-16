import argparse
import gffutils

VERSION = '0.2.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='Galaxy BioBlend workflow generation script.')
    parser.add_argument('--gffutils_db', default='gffutils_db', help='GffUtils DB')
    parser.add_argument('--cnv_infile', help='CNV TSV Infile')
    parser.add_argument('--cnv_denoised_infile', help='CNV DenoisedReadCounts Standardized TSV Infile')
    parser.add_argument('--genes', help='Gene list from which to produce output.')
    parser.add_argument('--tumor_pct', help='Tumor percentage.')
    parser.add_argument('--outfile', help='CNV enhanced output.')
    parser.add_argument('--outfile_genes', help='CNV enhanced output, gene-based.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

    args = parser.parse_args()
    return args

def coords_to_codons(start, end):
    """
    Given a start and end coordinate, convert to codons.
    :param start:
    :param end:
    :return:
    """
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
    """
    Picard Interval == coordinates, and any additional information that can be derived per interval
      whatver additional information is available in other columns
    generic, only include raw information in data structure
    """
    def __init__(self, line, titles):
        self.line = line
        self.titles = titles
        self.ival = self._make_ival()

    def _make_ival(self):
        this_ival = {}
        for i in range(len(self.titles)):
            this_ival[self.titles[i]] = self.line[i]
        return this_ival



class KdlCopyInterval(object):
    """
    Apply additional fields to PicardInterval ival's, to get list of fields specific to KDL copy number
    analysis.

    From gffutils:
        Gene Start
        Gene Stop
        Start Exon
        Stop Exon
        Gene Name
        Gene Chrom
    """
    def __init__(self, pic_ival, dbname, common_to_refseq, my_genes):
        # PicardInterval.__init__(self, ival, titles)
        # Placing these fields here, not consistent to all PicardIntervals.
        self.pic_ival = pic_ival
        self.ival = self.pic_ival.ival
        self.common_to_refseq = common_to_refseq
        self.chrom = common_to_refseq[self.pic_ival.ival['CONTIG']]
        self.start = int(self.pic_ival.ival['START'])
        self.end = int(self.pic_ival.ival['END'])
        self.my_genes = my_genes
        self.dbname = dbname
        # Modify CALL value
        if 'CALL' in self.pic_ival.ival:
            self.ival['CALL'] = self._call_field_convert(self.pic_ival.ival['CALL'])
        self._add_gffutil_field()

    def _add_gffutil_field(self):
        self.ival.update(GffUtilUtil(self.dbname, self.chrom, self.start, self.end).gffutil_rec)

    def _call_field_convert(self, call):
        """
        Convert the 0/+/- values to human readable values, listed in call_conv.
        :return:
        """
        call_conv = {'0': 'Neutral',
                     '+': 'Gain',
                     '-': 'Loss'}
        try:
            return call_conv[call]
        except:
            return None


class GeneImport(object):
    """

    """
    def __init__(self, filename):
        self.filename = filename
        self.genes = self._file_to_list()

    def _file_to_list(self):
        genes = []
        with open(self.filename, 'rU') as myfile:
            for line in myfile:
                if len(line.split('\t')) > 1:
                    raise Exception("Line has more than one column.")
                line = line.rstrip('\n')
                genes.append(line)
        return genes



class PicardIntervals(object):
    """
    Contains header information for the file. (header, titles)
    Contains PicardInterval regions.
    Takes a filename, and an interval object type (KdlCopyInterval, PicardInterval)
    """
    def __init__(self, filename, header_start='CONTIG'):
        self.filename = filename
        self.header_start = header_start
        self.header = self._assign_header()
        self.titles = self._assign_titles()
        self.regions = self._assign_regions()

    def _assign_header(self):
        """
        Place header information in a list.
        :return:
        """
        header = []
        with open(self.filename, 'rU') as myfile:
            for line in myfile:
                if line.startswith('@'):
                    header.append(line)
        return header

    def _assign_titles(self):
        """
        Place titles in a list.
        :return:
        """
        with open(self.filename, 'rU') as myfile:
            for line in myfile:
                if line.startswith(self.header_start):
                    return line.rstrip('\n').split('\t')

    def _assign_regions(self):
        """
        Assign PicardInterval regions.
        :return:
        """
        regions = []
        with open(self.filename, 'rU') as myfile:
            for line in myfile:
                if not line.startswith('@') and not line.startswith(self.header_start):
                    ival = PicardInterval(line.rstrip('\n').split('\t'), self.titles)
                    regions.append(ival)
        return regions

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

def flip_common_chroms(common_to_refseq):
    """
    Flip the keys/values, so we can easily match off the value instead.
    :param common_to_refseq:
    :return:
    """
    return dict((v,k) for k,v in common_to_refseq.items())

class PicIntWriter(object):

    """
    Pass filename, header (list of header lines), title(list of column headers), and intervals (PicardInterval format).
    """

    def __init__(self, filename, header, title, intervals):
        self.filename = filename
        self.header = header
        self.title = title
        self.intervals = intervals

    def write_me(self, write_order=None, write_header=True):
        """
        Signal to start the writing operation.  Pass write_order if wanted.
        :param write_order:
        :return:
        """
        with open(self.filename, 'w') as myfile:
            # Write header.
            if write_header:
                for entry in self.header:
                    myfile.write(entry)
                    if not entry.endswith('\n'):
                        myfile.write('\n')
            # If write_order is specified, write the column headings just as specified.  Otherwise
            # use the title as passed to PicIntWriter.
            if write_order:
                myfile.write('\t'.join(write_order))
            else:
                myfile.write('\t'.join(self.title))
            myfile.write('\n')

            all_write = []
            for entry in self.intervals:
                to_write = []
                for title in write_order:
                    if title in entry:
                        to_write.append(str(entry[title]))
                    else:
                        to_write.append('')
                all_write.append(to_write)

            for write in sorted(all_write, key=lambda x: (x[0], int(x[1]))):
                myfile.write('\t'.join(write))
                myfile.write('\n')

    def write_amp_codon(self):
        pass

class GeneCentricCnv(object):
    """
    Need new format, based on intervals.
    """
    def __init__(self, my_ints, den_ints, genes, dbname, common_to_refseq, tumor_pct):
        self.my_ints = my_ints
        self.den_ints = den_ints
        self.genes = genes
        self.dbname = dbname
        self.common_to_refseq = common_to_refseq
        self.refseq_to_common = flip_common_chroms(common_to_refseq)
        self.tumor_pct = tumor_pct
        self.my_gene_info = GeneInfo(self.dbname, self.common_to_refseq, self.genes)
        self.gene_info = self.my_gene_info.gene_info
        self.den_genes = self._create_den_genes()
        self.gene_ints = self._gene_based_create()

    def _create_den_genes(self):
        """
        Create structure with genes as keys, and (chrom,start,stop) pairs as values.
        :return:
        """
        den_genes = {}
        for entry in self.den_ints.regions:
            key = (self.common_to_refseq[entry.ival['CONTIG']], entry.ival['START'], entry.ival['END'])
            genes = self.my_gene_info.gene_from_region(key)
            if genes:
                for gene in genes:
                    if gene not in den_genes:
                        den_genes[gene] = [key]
                    else:
                        den_genes[gene].append(key)
        return den_genes

    def _set_num_points(self, gene, default, start_exon, end_exon, int_start, int_end):
        """
        Number of points, for the gene-based output, will be based on number of regions per gene listed
        in den_genes.  If the call is broken in the middle of a gene, will need to look at coordinates to determine
        number of points.
        :return:
        """
        if gene in self.den_genes:
            if start_exon == 'NA' and end_exon == 'NA':
                return len(self.den_genes[gene])
            elif start_exon == 'NA':
                for ival in sorted(self.den_genes[gene]):
                    if ival[2] == int_end:
                        return self.den_genes[gene].index(ival) + 1
            elif end_exon == 'NA':
                for ival in sorted(self.den_genes[gene]):
                    if ival[1] == int_start:
                        return len(self.den_genes[gene]) - self.den_genes[gene].index(ival)
            else:
                for ival in sorted(self.den_genes[gene]):
                    if ival[1] == int_start:
                        start_ind = self.den_genes[gene].index(ival)
                    if ival[2] == int_end:
                        end_ind = self.den_genes[gene].index(ival)
                return end_ind - start_ind + 1
        else:
            try:
                print(self.den_genes[gene])
                print(len(self.den_genes[gene]))
            except:
                pass
            print(start_exon)
            print(end_exon)
            print(int_start)
            print(int_end)
            return default

    def _calc_raw_copy_number(self, val):
        """
        Calculate raw copy number, from MEAN_LOG2_COPY_RATIO.
        2*(2^MEAN_LOG2_COPY_RATIO)
        :return:
        """
        return 2*(2**float(val))

    def _calc_tumor_copy(self, tumor_pct, val):
        """
        Calculate tumor corrected copy number.

        :return:
        """
        val = self._calc_raw_copy_number(val)
        t = int(tumor_pct)/100.0
        n = 1.0 - t
        return (val-(2*n))/t

    def _calc_tumor_copy_ratio(self, tval, nval):
        """
        Calculate ratio of tval over nval.
        :param tval:
        :param nval:
        :return:
        """
        return float(tval)/float(nval)

    def _gene_based_create(self):
        gene_output = []
        for entry in self.my_ints:
            for gene in self.genes:
                if gene in entry['GENES'].split(','):
                    to_write = {'GENE_CHROM': str(self.refseq_to_common[self.gene_info[gene]['GENE_CHROM']]),
                                'GENE_START': str(self.gene_info[gene]['GENE_START']),
                                'GENE_STOP': str(self.gene_info[gene]['GENE_STOP']),
                                'GENE': gene,
                                'NUM_POINTS_SEGMENT': str(entry['NUM_POINTS_COPY_RATIO']),
                                'MEAN_LOG2_COPY_RATIO': str(entry['MEAN_LOG2_COPY_RATIO']),
                                'RAW_COPY_NUMBER': self._calc_raw_copy_number(str(entry['MEAN_LOG2_COPY_RATIO'])),
                                'TUMOR_COPY_NUMBER': self._calc_tumor_copy(self.tumor_pct, str(entry['MEAN_LOG2_COPY_RATIO'])),
                                'CALL': entry['CALL']}
                    if gene == entry['START_GENE']:
                        to_write['START_EXON'] = str(entry['START_EXON'])
                    else:
                        to_write['START_EXON'] = 'NA'
                    if gene == entry['END_GENE']:
                        to_write['END_EXON'] = str(entry['END_EXON'])
                    else:
                        to_write['END_EXON'] = 'NA'

                    to_write['NUM_POINTS_GENE'] = str(self._set_num_points(gene,
                                                                          entry['NUM_POINTS_COPY_RATIO'],
                                                                          to_write['START_EXON'],
                                                                          to_write['END_EXON'],
                                                                          entry['START'],
                                                                          entry['END']))
                    to_write['TUMOR_COPY_RATIO'] = str(self._calc_tumor_copy_ratio(to_write['TUMOR_COPY_NUMBER'],
                                                                                   to_write['RAW_COPY_NUMBER']))
                    gene_output.append(to_write)
        return gene_output


class GeneInfo(object):
    """
    From a list of genes (HGNC), find coordinates, any additional information we want that can be obtained from
    gffutils.
    """
    def __init__(self, dbname, common_to_refseq, genes):
        self.genes = genes
        self.db = gffutils.FeatureDB(dbname)
        self.common_to_refseq = common_to_refseq
        self.gene_info = self._get_gene_info()

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

    def gene_from_region(self, region):
        """
        Given a region (chrom, start, stop), return the gene name(s) overlapping that region.
        :return:
        """
        genes = []
        for entry in self.db.region(region, featuretype='gene'):
            genes.extend(entry.attributes['gene'])
        return genes


def main():

    args = supply_args()
    # Produce a mapping of RefSeq IDs to common chromosome names.
    common_to_refseq = get_common_chroms(args.gffutils_db)
    # Initialize gffutils db.
    dbname = args.gffutils_db
    # Produce list of reportable genes, or just any genes we care about.
    genes = GeneImport(args.genes).genes
    # Split input Picard Intervals formatted file in to headers, titles, and interval regions.
    pic_ints = PicardIntervals(args.cnv_infile)
    den_ints = PicardIntervals(args.cnv_denoised_infile)

    write_order = ['CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO', 'CALL',
                   'START_GENE', 'START_EXON', 'END_GENE', 'END_EXON']
    gene_write_order = ['GENE_CHROM', 'GENE_START', 'GENE_STOP', 'GENE', 'NUM_POINTS_SEGMENT',
                        'NUM_POINTS_GENE', 'MEAN_LOG2_COPY_RATIO', 'RAW_COPY_NUMBER', 'TUMOR_COPY_NUMBER', 'CALL']
    to_write = []
    for ival in pic_ints.regions:
        to_write.append(KdlCopyInterval(ival, dbname, common_to_refseq, args.genes).ival)
    gene_to_write = GeneCentricCnv(to_write, den_ints, genes, dbname, common_to_refseq, args.tumor_pct).gene_ints
    PicIntWriter(args.outfile, pic_ints.header, pic_ints.titles, to_write).write_me(write_order)
    PicIntWriter(args.outfile_genes, pic_ints.header, pic_ints.titles, gene_to_write).write_me(gene_write_order, write_header=False)


main()