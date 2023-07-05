#!/usr/bin/env python

import unittest
import os
import tempfile
from filecmp import cmp


class PonFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        base_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../pon_filter'))
        cls.test_data_dir = os.path.join(base_dir, 'tests', 'test_data')

    def test_pon_filter(self):
        tempdir = tempfile.TemporaryDirectory()
        os.symlink(os.path.join(self.test_data_dir, 'input.vcf'), os.path.join(tempdir.name, 'input.vcf'))
        os.symlink(os.path.join(self.test_data_dir, 'pon.vcf'), os.path.join(tempdir.name, 'pon.vcf'))
        cmd = [
            os.path.join(tempdir.name, 'input.vcf'),
            os.path.join(tempdir.name, 'pon.vcf'),
            "output.vcf"
            "--min_cnt" "40"
            "--pon_flag" "PON"
            "--pon_flag_filtered" "PON_FILT"
            "--bkgd_avg" "0.2"
            "--bkgd_std" "0.06"
            "--bkgd_min_cnt" "4"
            "--pon_flag_above_bkgd" "PON_OV"
        ]
        pon_filt.pon_filter.main(cmd)
        self.assertTrue(cmp(os.path.join(self.test_data_dir, 'output.vcf'), os.path.join(tempdir.name, 'output.vcf')))
        tempdir.cleanup()
