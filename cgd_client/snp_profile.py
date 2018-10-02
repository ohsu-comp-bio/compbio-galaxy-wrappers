from scipy.stats import binom_test

import json
import vcf

class SnpProfile(object):
    """
    Structure of json sent back and forth to/from CGD.
    {"message":"ok","items":[{"chromosome":"22","position":700,"genotype":0},
    {"chromosome":"3","position":500,"genotype":1},{"chromosome":"7","position":600,"genotype":2}]}

    This works (GET):
    letaw@exaclinical:~/clinical/users/letaw/compbio-galaxy-wrappers/snp_profile$ java8 -jar /home/exacloud/clinical/installedTest/cgd_client-1.2.4.jar -c /home/exacloud/clinical/installedTest/cgd_client.properties -u https://kdlwebuser01.ohsu.edu/cgd_next/service/run/180408_NS500390_0222_AHMMV5BGX5/barcodeId/A01/snpProfile
    {"message":"ok","items":[]}

    (POST)
    letaw@exaclinical:~/clinical/users/letaw/compbio-galaxy-wrappers/snp_profile$ java8 -jar /home/exacloud/clinical/installedTest/cgd_client-1.2.4.jar -c /home/exacloud/clinical/installedTest/cgd_client.properties -u https://kdlwebuser01.ohsu.edu/cgd_next/service/run/180408_NS500390_0222_AHMMV5BGX5/barcodeId/A01/snpProfile -j blah -d
    URL endpoint: https://kdlwebuser01.ohsu.edu/cgd_next/service/run/180408_NS500390_0222_AHMMV5BGX5/barcodeId/A01/snpProfile
    {"message":"ok","items":[]}


    """
    def __init__(self, filename):
        handle = open(filename, 'r')
        self.vcf_read = vcf.Reader(handle)
        self.geno_items = self._vcf_parse()

    def _vcf_parse(self):
        """
        Get coordinates and genotypes from the input VCF file.
        :return:
        """
        geno_items = []
        for record in self.vcf_read:
            chrom = str(record.CHROM)
            pos = int(record.POS)
            if len(record.samples) == 1:
                geno = record.samples[0].gt_type
                if geno is None:
                    geno = -1
            else:
                raise Exception("The input VCF should only have one sample column, this one has " + str(len(record.samples)))
            geno_items.append({"chromosome": chrom, "position": pos, "genotype": geno})

        return geno_items


class CompareProfiles(object):
    """
    Stuff related to comparing two SNP profiles.
    proft = template profile
    profp = patient profile
    """
    def __init__(self, proft, profp):
        # self.proft = self._create_dict(json.loads(proft))
        # self.profp = self._create_dict(json.loads(profp))
        self.proft = self._create_dict(proft)
        self.profp = self._create_dict(profp)
        self._compare_them()
        self.for_cgd = self.dict_to_json(self.new_snps)
        self._perform_binom()

    def _create_dict(self, prof):
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


    def dict_to_json(self, to_json):
        """
        Put the dictionary that will be sent back in a format that the CGD will accept.
        The format we want is: [{"chromosome": <str>, "position": <int>, "genotype": <int>}, {}, ...]
        Needs to be dumped to json.
        :return:
        """
        new_json = []
        for pos, geno in to_json.iteritems():
            chrom = str(pos[0])
            coord = int(pos[1])
            to_add = {"chromosome": chrom, "position": coord, "genotype": geno}
            new_json.append(to_add)
        return json.dumps(new_json)


    def _count_matches(self):
        """
        Look for the valid genotypes contained within the list of SNPs, and count them.
        :return:
        """
        return len([x for x in (list(set(self.proft) & set(self.profp))) if x != '-1'])


    def _perform_binom(self, prob=0.01, pfail=0.05):
        """
        Look at total and mismatch values from compare_them, and find a p-value.
        :return:
        """
        pvalue = binom_test(self.mismatch, self.total, prob)
        if pvalue < pfail:
            raise Exception("Null hypothesis(p=" + str(pfail) + ") of patient profiles being the same rejected with p-value: " + str(pvalue))
        else:
            print("Null hypothesis(p=" + str(pfail) + ") of patient profiles being the same accepted with p-value: " + str(pvalue))


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

        for pos, geno in self.profp.iteritems():
            if geno != -1:
                if pos in self.proft:
                    # Situation 2
                    if geno != self.proft[pos]:
                        self.mismatch += 1
                    self.total += 1
                else:
                    # Situation 4
                    self.new_snps[pos] = geno
