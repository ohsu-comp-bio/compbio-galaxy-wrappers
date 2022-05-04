from collections import OrderedDict


class VcfReader(object):
    """
    Read in a VCF, provide basic header and variant structures.
    """

    def __init__(self, infile):
        self.infile = open(infile, 'rU')
        self.header, self.vcf = self._parse_infile()

    def _parse_infile(self):
        """
        Turn the input VCF in to a dictinoary of variants and list of header
        lines.
        :return:
        """
        header = []
        vcf = OrderedDict()
        for line in self.infile:
            if line.startswith('#'):
                header.append(line)
            else:
                line = line.rstrip('\n').split('\t')
                chrom = line[0]
                coord = line[1]
                vcf[(chrom, coord)] = line

        return header, vcf

    def _find_header_index(self, hfield):
        """
        Find the index of where a new header record should be inserted.
        :return:
        """
        found = False
        for entry in self.header:
            if entry.split('=')[0][2:] == hfield and found is False:
                found = True
            elif entry.split('=')[0][2:] != hfield and found is True:
                return self.header.index(entry)

    def add_header(self, hfield, hid, hnum, htype, hdesc):
        """
        Add the new header line to the VcfReader structure.
        ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of
        the event in the tumor">\n
        :return:
        """
        to_write = '##{}=<ID={},Number={},Type={},Description="{}">\n'.format(
            hfield, hid, hnum, htype, hdesc)
        my_index = self._find_header_index(hfield)
        self.header.insert(my_index, to_write)

    def add_format_field(self, chrom, coord, fname, fvalue):
        """
        Add an entry to the VCF FORMAT and Sample fields.
        :return:
        """
        self.vcf[(chrom, coord)][8] = ':'.join([self.vcf[(chrom, coord)][8],
                                                fname])
        self.vcf[(chrom, coord)][9] = ':'.join([self.vcf[(chrom, coord)][9],
                                                str(fvalue)])


class VcfWriter(object):
    """
    Write a VCF.
    """

    def __init__(self, outfile, vcffile):
        self.outfile = open(outfile, 'w')
        self.vcffile = vcffile
        self._write_me()

    def _write_me(self):
        for entry in self.vcffile.header:
            self.outfile.write(entry)
        for entry in self.vcffile.vcf.values():
            self.outfile.write('\t'.join(entry))
            self.outfile.write('\n')

        self.outfile.close()
