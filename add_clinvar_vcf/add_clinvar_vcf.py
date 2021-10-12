#!/usr/bin/env python

# Currently doing things this way since I'm having trouble getting VariantAnnotator from GATK to do what I want.
# In situation where ClinVar has multiple entries at same coordinate, GATK tool applies only annotation from first
# entry.  If your VCF contains multiple alternates for a record, this will cause mis-annotation based on first entry
# found in ClinVar file.  If you split your alternate alleles out in to multiple VCF entries, you will end up with
# no annotation instead of the incorrect one, which is also incorrect.

import argparse
import gzip

VERSION = '0.1.1'


def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input', help='Input VCF to annotate with ClinVar.')
    parser.add_argument('--output', help='Output VCF')
    parser.add_argument('--info', nargs='*', help='INFO annotation to pull.')
    parser.add_argument('--resource', help='Resource VCF')
    parser.add_argument('--resource_lbl', help='Resource VCF Label')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

class ClinVar(object):
    def __init__(self, filename, annot=None, prefix=None):
        try:
            self.clinvar = gzip.open(filename, 'rt')
        except:
            self.clinvar = open(filename, 'rU')
        self.annot = annot
        self.prefix = prefix
        self.headers = []
        self.clinvar_vars = self.parse_clinvar()

    def parse_clinvar(self):
        """
        Grab necessary information from ClinVar VCF.
        :return:
        """
        clinvar = {}
        with self.clinvar as myfile:
            for line in myfile:
                if not line.startswith('#'):
                    line = line.rstrip('\n').split('\t')
                    chrom = line[0]
                    coord = line[1]
                    ref = line[3]
                    alt = line[4]
                    info = line[7]
                    uniq_key = (chrom, coord, ref, alt)
                    clinvar[uniq_key] = []
                    if self.annot:
                        for ann in self.annot:
                            clinvar[uniq_key].append(self._find_info(ann, info))
                else:
                    self._header_grab(line.rstrip('\n'))
        return clinvar

    def _find_info(self, label, info):
        """
        Find something from the VCF INFO section.
        :return:
        """
        for entry in info.split(';'):
            if label == entry.split('=')[0]:
                return entry.split('=')[1]
        return ''

    def _create_new_info(self, values):
        """
        Create the new entry that will be tacked on to the existing INFO column.
        :return:
        """
        for i in range(len(self.annot)):
            new_label = ';' + self.prefix + '.' + self.annot[i]
            new_value = values[i]
            if new_value:
                yield '='.join([new_label, new_value])

    def _header_modify(self, header):
        """
        Change ID from 'AF' to '<resource_lbl>.AF'.
        :return:
        """
        new_label = 'ID=' + self.prefix + '.'
        return header.replace('ID=', new_label)

    def _header_grab(self, header):
        """
        Find and retain the header that matches the annotation we want.
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency among genotypes, for each ALT allele, in the same order as listed">
        :return:
        """
        if header.startswith('##INFO'):
            new_header = header.split('<')
            for entry in new_header[1].split(','):
                label = entry.split('=')[0]
                value = entry.split('=')[1]
                if value in self.annot and label == 'ID':
                    self.headers.append(header)
                    return None
                elif label == 'ID':
                    return None
        return None


    def vcf_rewrite(self, infile, outfile):
        """

        :return:
        """
        output = open(outfile, 'w')
        add_headers = True
        with open(infile, 'rU') as myfile:
            for line in myfile:
                if not line.startswith('#') and not line == '\n':
                    line = line.rstrip('\n').split('\t')
                    chrom = line[0]
                    coord = line[1]
                    ref = line[3]
                    alt = line[4]
                    info = line[7]
                    uniq_key = (chrom, coord, ref, alt)
                    if uniq_key in self.clinvar_vars:
                        for entry in self._create_new_info(self.clinvar_vars[uniq_key]):
                            info += entry
                    to_write = line[:7]
                    to_write.append(info)
                    to_write.extend(line[8:])
                    output.write('\t'.join(to_write))
                    output.write('\n')
                else:
                    if line.startswith('##INFO') and add_headers:
                        for header in self.headers:
                            output.write(self._header_modify(header))
                            output.write('\n')
                        add_headers = False
                    output.write(line)
        output.close()

def main():
    args = supply_args()
    to_find = args.info
    prefix = args.resource_lbl
    clinvar = ClinVar(args.resource, to_find, prefix)
    clinvar.vcf_rewrite(args.input, args.output)


if __name__ == "__main__":
    main()


