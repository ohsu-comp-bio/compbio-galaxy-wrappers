class VcfWriter(object):
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
