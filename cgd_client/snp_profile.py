from scipy.stats import binom_test
import vcfpy

VERSION = '0.2.2'


class SnpProfile(object):
    """
    Structure of json sent back and forth to/from CGD.
    {"message":"ok","items":[{"chrompposome":"22","position":700,"genotype":0},
    {"chromosome":"3","position":500,"genotype":1},{"chromosome":"7","position":600,"genotype":2}]}

    This works (GET):
    java8 -jar ~/cgd_client-1.2.4.jar -c ~/cgd_client.properties -u <url>
    {"message":"ok","items":[]}

    (POST)
    java8 -jar ~/cgd_client-1.2.4.jar -c ~/cgd_client.properties -u <url> -j blah -d
    {"message":"ok","items":[]}
    """
    def __init__(self, filename):

        self.vcf_read = vcfpy.Reader.from_path(filename)
        self.ab_dict = {}
        self.geno_items = self._vcf_parse()

    def _vcf_parse(self):
        """
        Get coordinates and genotypes from the input VCF file.
        :return:
        """
        geno_items = []
        for record in self.vcf_read:
            samp_cnt = len([call for call in record.calls])
            chrom = str(record.CHROM)
            pos = int(record.POS)
            if samp_cnt == 1:
                geno = self._get_geno([call.data.get('GT') for call in record.calls][0])
                if geno != -1:
                    ref = [call.data.get('AD') for call in record.calls][0][0]
                    alt = [call.data.get('AD') for call in record.calls][0][1]
                    ab = self._calc_ab(ref, alt, geno)
                    self.ab_dict[(chrom, pos)] = ab
                    geno = self._assess_ab(ab, geno)
            else:
                raise Exception("The input VCF should only have one sample column, this one has " +
                                str(len(record.samples)))
            geno_items.append({"chromosome": chrom, "position": pos, "genotype": geno})

        return geno_items

    @staticmethod
    def _get_geno(geno):
        """
        For genotypes, set values 0, 1, and 2.
        :return:
        """
        if geno == '0/0':
            return 0
        elif geno == '0/1':
            return 1
        elif geno == '1/1':
            return 2
        else:
            return -1

    @staticmethod
    def _assess_ab(ab, geno, hom_thresh=0.01, het_thresh=0.1):
        """
        """
        if ab is None or geno is None:
            return -1
        elif ab > hom_thresh and geno == 0:
            return -1
        elif ab > het_thresh and geno == 1:
            return -1
        elif ab > hom_thresh and geno == 2:
            return -1
        else:
            return geno

    @staticmethod
    def _calc_ab(ref, alt, geno):
        """
        Calculate the allele balance.
        :return:
        """
        if geno == 0:
            return float(alt / (ref + alt + 0.0))
        elif geno == 1:
            return abs(0.5 - (float(alt / (ref + alt + 0.0))))
        elif geno == 2:
            return float(ref / (ref + alt + 0.0))
        else:
            return None


class CompareProfiles(object):
    """
    Stuff related to comparing two SNP profiles.
    proft = template profile
    profp = patient profile
    """
    def __init__(self, proft, profp):
        self.proft = self._create_dict(proft)
        self.profp = self._create_dict(profp)
        self.mm_loc = []
        self._compare_them()
        self.for_cgd = self.dict_to_json(self.new_snps)
        self.pvalue = self._perform_binom()

    @staticmethod
    def _create_dict(prof):
        """
        Make a dictionary, to facilitate easy comparison between profiles.
        {(chrom, pos): genotype}
        :return:
        """
        this_dict = {}
        for entry in prof:
            chrom = entry["chromosome"]
            pos = entry["position"]
            geno = entry["genotype"]
            this_dict[(chrom, pos)] = geno
        return this_dict

    @staticmethod
    def dict_to_json(to_json):
        """
        Put the dictionary that will be sent back in a format that the CGD will accept.
        The format we want is: [{"chromosome": <str>, "position": <int>, "genotype": <int>}, {}, ...]
        Needs to be dumped to json.
        :return:
        """
        new_json = []
        for pos, geno in to_json.items():
            chrom = str(pos[0])
            coord = int(pos[1])
            to_add = {"chromosome": chrom, "position": coord, "genotype": geno}
            new_json.append(to_add)
        return new_json

    def _count_matches(self):
        """
        Look for the valid genotypes contained within the list of SNPs, and count them.
        :return:
        """
        return len([x for x in (list(set(self.proft) & set(self.profp))) if x != '-1'])

    def _perform_binom(self, prob=0.001, pfail=0.05):
        """
        Look at total and mismatch values from compare_them, and find a p-value.
        :return:
        """
        pvalue = binom_test(self.mismatch, self.total, prob, alternative='greater')
        if pvalue <= pfail:
            print(f"Null hypothesis(p={str(pfail)}) of sample profiles being the same is rejected: {str(pvalue)}")
        else:
            print(f"Null hypothesis(p={str(pfail)}) of sample profiles being the same is confirmed: {str(pvalue)}")
        return pvalue

    def _compare_them(self):
        """
        Compare the two profiles.  Different things can happen:
        1. entry is in both, and proft == profp
            A: count one match for this
        2. entry is in both, but proft != profp
            A: count one mismatch for this
        3. entry is in proft, but not in profp
            A: Not relevant, do not count this event.
        4. entry is in profp, but not in proft
            A: This needs to be added to the profile, but it is not to be counted.

        Any instance of a '-1' value as a genotype should be ignored, and not counted.
        NOTE: Problem with -1, if this is set up front in CGD profile, can't be overwritten.  No bueno.

        :return:
        """
        self.total = 0
        self.mismatch = 0
        self.new_snps = {}

        for pos, geno in self.profp.items():
            if geno != -1:
                if pos in self.proft:
                    # Situation 2
                    if geno != self.proft[pos] and self.proft[pos] != -1:
                        self.mismatch += 1
                        self.mm_loc.append(pos)
                    self.total += 1
                else:
                    # Situation 4
                    self.new_snps[pos] = geno
