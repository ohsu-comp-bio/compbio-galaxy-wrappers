import amplicon_assess
import pysam
import unittest
import vcfpy


class VrntMetricsTestCase(unittest.TestCase):

    def setUp(self):
        self.samfile = pysam.AlignmentFile("test/umi_dedup.bam", "rb")
        self.primers = amplicon_assess.QiagenPrimers("test/bed_test.bed").my_primers
        print(self.primers)
        self.reader = vcfpy.Reader.from_path("test/vcf_test.vcf.gz")

    def test_check_all_primers(self):
        vrnt = ("3", 178952085)
        primers = {('3', 178952219): 1}
        metrics = amplicon_assess.VrntMetrics(self.samfile, vrnt, primers)
        print(metrics.raw_metrics)
        print(metrics.primers)
        coord = ("3", 178952200)
        self.assertEqual(178952219, metrics._check_all_primers(coord))

    def tearDown(self):
        self.samfile.close()
        self.reader.close()


if __name__ == '__main__':
    unittest.main()
