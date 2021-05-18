from string import Template

class VcfWriter:
    """
    Write a VCF.  A header needs to be supplied, which is simply a list of all header lines.
    Also, the body needs to be supplied, which is additionally just a list of ordered lines.
    Creating these formatted lines is done by VcfRecBase.
    """
    def __init__(self, filename, vcf, header):
        self.vcf = vcf
        self.filename = filename
        self.header = header

    def write_me(self):
        """
        Perform writing of MyVcf structure.
        :return:
        """
        with open(self.filename, 'w') as outvcf:
            for entry in self.header:
                outvcf.write(entry)
                outvcf.write('\n')
            for entry in self.vcf:
                to_write = '\t'.join(entry)
                outvcf.write(to_write)
                outvcf.write('\n')


class VcfHeader:
    """
    Represents the VCF header.
    """
    def __init__(self, raw_header):
        self.raw_header = raw_header
        self.samp_idx = self._find_sample_header()

    def _add_entry(self):
        """
        Add an entry to the header.
        FILTER entry:
        ##FILTER=<ID=LowQual,Description="Low quality">
        :return:
        """
        pass

    def _find_sample_header(self):
        """
        Figure out the index for the CHROM header.
        :return:
        """
        for entry in self.raw_header:
            if entry.startswith('#CHROM'):
                return self.raw_header.index(entry)
        return None

    def _find_entry_idx(self, htype):
        idx = None
        for entry in self.raw_header:
            search_for = '##' + htype
            if entry.startswith(search_for):
                idx = self.raw_header.index(entry)
        return idx

    def add_header_line(self, htype, header_dict):
        """
        Add a line to the header.
        :return:
        """
        htype_idx = self._find_entry_idx(htype)
        if htype_idx:
            self.raw_header.insert(htype_idx, self._filter_create(htype, header_dict))
        else:
            self.raw_header.insert(self.samp_idx, self._filter_create(htype, header_dict))


    def _filter_create(self, htype, header_dict):
        """
        Create FILTER lines to be added for each new label.
        :return:
        """
        filt_tmpls = {'FILTER': Template('##FILTER=<ID=${id},Description=\"${desc}\">')}
        return filt_tmpls[htype].substitute(header_dict)
