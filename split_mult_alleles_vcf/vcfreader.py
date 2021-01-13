from collections import OrderedDict
import re


class VcfAllele(object):
    """
    Define what a VCF allele looks like.
    ACTGN
    """

    def __init__(self, alleles, is_alt=False):
        self.is_alt = is_alt
        self.alleles = [self._check_allele(allele) for allele in alleles.split(',')]
        self.pr_alleles = ','.join(self.alleles)

    def _set_valid_regex(self):
        """
        Choose the set of legal regexes that will be checked against.
        :return:
        """
        if self.is_alt:
            return [r'^[ATCGNatcgn]+$', r'^[\*]$', r'^<[A-Za-z0-9]+>$']
        else:
            return [r'^[ATCGNatcgn]+$']

    def _check_allele(self, allele):
        """
        Check to see whether the given allele string follows spec, otherwise return None.
        allele with recognized bases: ^[ATCGN]+$
        asterisk allele: ^[*]$
        bracketed allele: ^<[A-Za-z0-9]+>$ (More characters are allowed, though definitely not angle brackets inside the string.
        :return:
        """
        valid = self._set_valid_regex()
        if any(re.match(regex, allele) for regex in valid):
            if allele.startswith('<'):
                return allele
            else:
                return allele.upper()
        else:
            if allele == '.':
                return None
            else:
                raise Exception('The allele \"' + str(allele) + '\" does not comply with the specification.')


class VcfFilt(object):
    def __init__(self, filt):
        self.filt = (''.join(filt.split())).split(';')
        self.pr_filt = ';'.join(self.filt)


class VcfInfo(object):
    def __init__(self, info):
        self.pr_info = "".join(info.split())
        self.info = self._create_info(self.pr_info)

    @staticmethod
    def _create_info(info):
        """
        Take the INFO string and provide it in dictionary format.
        :return:
        """
        new_info = info.split(';')
        info_dict = OrderedDict()
        for pair in new_info:
            key = pair.split('=')[0]
            try:
                val = pair.split('=')[1]
            except IndexError:
                val = None
            info_dict[key] = val

        return info_dict


class VcfSamples(object):
    def __init__(self, frmt, samps):
        self.frmt = (''.join(frmt.split())).split(':')
        self.my_samps = self._prep_samps(samps, self.frmt)

    @staticmethod
    def _prep_samps(samps, frmt):
        """
        Get the samples field setup...
        :return:
        """
        my_samps = []
        for samp in samps:
            samp_dict = OrderedDict()
            smp = ''.join(samp.split()).split(':')
            idx = 0
            for val in smp:
                samp_dict[frmt[idx]] = val
                idx += 1
            my_samps.append(samp_dict)

        return my_samps


class VcfRecBase(object):
    """
    Basic parsing of VCF records.
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	RDM-0C19G1-1
    1	16200729	.	A	AT	.	m2	DP=1180;ECNT=1;POP_AF=5e-08;RPA=9;RU=T;STR;TLOD=25.6;OLD_VARIANT
    =1:16200729:AT/ATT	GT:AD:AF:DP:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:ORIGINAL_CONTIG_MISMATCH:SA_MAP_AF:SA_POST_PROB	0/1:906,38:0.031
    :1078:906,38:0,0:36,32:202,190:60,60:38,35:0:0.121,0.121,0.129:0.011,0.002304,0.986

    all fields can be '.'

    CHROM - string
    POS - int
    ID - string
    REF - string (ATCGN)
    ALT - comma sep (ATCGN*), or <> enclosed
    QUAL - float
    FILTER - ;-sep, or PASS
    INFO - key=val pairs, commas only for val delimiters, ;-sep
    FORMAT - :-sep strings
    SAMPLE - :-sep, val per FORMAT field.  Nothing is required, but GT must be first if it exists.  GT can be / or |


    """

    def __init__(self, rec):

        self.rec = rec.rstrip('\n').split('\t')

        try:
            self.CHROM = str(self.rec[0])
        except ValueError:
            self.CHROM = None

        try:
            self.POS = int(self.rec[1])
        except ValueError:
            self.POS = None

        try:
            self.ID = str(self.rec[2])
        except ValueError:
            self.ID = None

        try:
            self.REF = VcfAllele(str(self.rec[3]), is_alt=False).alleles
        except ValueError:
            self.REF = None

        try:
            self.ALT = VcfAllele(str(self.rec[4]), is_alt=True).alleles
        except ValueError:
            self.ALT = None

        # QUAL may need to be manipulated, but it should be passed here as a string.
        try:
            self.QUAL = self.rec[5]
        except ValueError:
            self.QUAL = None

        try:
            self.FILTER = VcfFilt(self.rec[6]).filt
        except ValueError:
            self.FILTER = None

        try:
            self.INFO = VcfInfo(self.rec[7]).info
        except ValueError:
            self.INFO = None

        try:
            self.my_samps = VcfSamples(self.rec[8], self.rec[9:])
        except IndexError:
            self.my_samps = None

        try:
            self.FORMAT = self.my_samps.frmt
            self.samples = self.my_samps.my_samps
        except (ValueError, AttributeError):
            self.FORMAT = None
            self.samples = None

        # Use this as a key for the vcf dict.
        self.uniq_key = (self.CHROM, self.POS, tuple(self.REF), tuple(self.ALT))

    def print_rec(self):
        """
        Create a line that conforms to VCF spec.

        Chromosome: 17
        Position: 37884037
        ID: rs1058808
        Reference Allele: ['C']
        Alternate Allele(s): ['G']
        Quality: None
        Filter: ['panel_of_normals']
        Info: OrderedDict([('DB', None), ('ECNT', '1'), ('HCNT', '11'), ('MAX_ED', '.'), ('MIN_ED', '.'), ('NLOD', '0.00'), ('TLOD', '21596.05')])
        Format String: ['GT', 'AD', 'AF', 'ALT_F1R2', 'ALT_F2R1', 'FOXOG', 'QSS', 'REF_F1R2', 'REF_F2R1']
        Samples: OrderedDict([('GT', '0/1'), ('AD', '1037,10894'), ('AF', '0.837'), ('ALT_F1R2', '5397'), ('ALT_F2R1', '5497'), ('FOXOG', '0.505'), ('QSS', '35288,371861'), ('REF_F1R2', '529'), ('REF_F2R1', '508')])

        :return:
        """
        write_line = [self.CHROM, str(self.POS), self.ID, self.REF[0], self.ALT[0]]

        if self.QUAL:
            write_line.append(str(self.QUAL))
        else:
            write_line.append(self.QUAL)

        filt = ';'.join(self.FILTER)
        write_line.append(filt)

        info_str = []
        for key, val in self.INFO.items():
            if val:
                to_add = '='.join([key, val])
            else:
                to_add = key
            info_str.append(to_add)

        write_line.append(';'.join(info_str))

        if self.my_samps:
            frmt_str = []
            for key, val in self.samples[0].items():
                frmt_str.append(key)
            write_line.append(':'.join(frmt_str))

            for samp in self.samples:
                samp_str = []
                for val in samp.values():
                    samp_str.append(val)
                write_line.append(':'.join(samp_str))

        return self._dot_for_none(write_line)

    def _dot_for_none(self, val):
        """
        If the value is None, return a dot.
        :return:
        """
        for entry in val:
            if entry == None:
                ind = val.index(entry)
                val[ind] = '.'
        return val

    def print_me(self):
        """
        Print values of each variable.
        :return:
        """
        print("Chromosome: {0}".format(self.CHROM))
        print("Position: {0}".format(self.POS))
        print("ID: {0}".format(self.ID))
        print("Reference Allele: {0}".format(self.REF))
        print("Alternate Allele(s): {0}".format(self.ALT))
        print("Quality: {0}".format(self.QUAL))
        print("Filter: {0}".format(self.FILTER))
        print("Info: {0}".format(self.INFO))
        print("Format String: {0}".format(self.FORMAT))
        print("Samples: {0}".format(self.samples))


class VcfHeader(object):
    """

    """

    def __init__(self, header):
        self.header = header
        self.info_number = self._prep_info_number(self.header, '##INFO')
        self.samples_number = self._prep_info_number(self.header, '##FORMAT')
        self.sample_list = self._get_sample_id_list(header)

    @staticmethod
    def _prep_info_number(header, start_str):
        """
        Get all of the INFO fields, and provide them in a usable form.
        :return:
        """
        info = {}
        for entry in header:
            if entry.startswith(start_str):
                for field in entry.split('<')[1].rstrip('\n>').split(','):
                    key = field.split('=')[0]
                    if key == 'ID':
                        key_val = field.split('=')[1]
                    elif key == 'Number':
                        info[key_val] = field.split('=')[1]

        return info

    @staticmethod
    def _get_sample_id_list(header):
        """
        Create a list of all samples in the VCF.
        :return:
        """

        for line in header:
            if line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[9:]

        return samples


class VcfReader(object):
    """
    Hold VCF records, manage access to them based in chrom, pos, ref, alt.
    """

    def __init__(self, filename):
        self.filename = filename
        self.raw_header = []
        self.myvcf = self._create_vcf()
        self.header = VcfHeader(self.raw_header)
        self.info_number = self.header.info_number
        self.samples_number = self.header.samples_number

    def _create_vcf(self):
        """
        Create structures to hold the VCF data.
        :return:
        """
        myvcf = OrderedDict()
        with open(self.filename, 'r') as vcf:
            for line in vcf:
                if not line.startswith('#'):
                    vrnt = VcfRecBase(line)
                    if vrnt.uniq_key not in myvcf:
                        myvcf[vrnt.uniq_key] = vrnt
                    else:
                        print("Duplicate entry found: " + str(vrnt.uniq_key) + " ... ignoring.")
                else:
                    self.raw_header.append(line.rstrip('\n'))
        return myvcf
