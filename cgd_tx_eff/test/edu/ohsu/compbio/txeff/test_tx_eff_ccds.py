import unittest
from edu.ohsu.compbio.txeff.tx_eff_ccds import TxEffCcds
from edu.ohsu.compbio.txeff.variant_transcript import VariantTranscript

class TxEffCcdsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tx_eff_ccds = TxEffCcds(None)

    def test__get_preferred_refseq_transcript_different_accesion(self):
        vt0 = VariantTranscript('1', 123, 'C', 'G')
        vt0.refseq_transcript = 'NM_000051.4'

        vt1 = VariantTranscript('1', 123, 'C', 'G')
        vt1.refseq_transcript = 'NM_001351834.2'

        preferred = self.tx_eff_ccds._get_preferred_refseq_transcript(vt0, vt1)
        self.assertEqual(preferred.refseq_transcript, 'NM_000051.4', "Earliest accession should be preferred")

    def test__get_preferred_refseq_transcript_different_version(self):
        vt0 = VariantTranscript('1', 123, 'C', 'G')
        vt0.refseq_transcript = 'NM_000051.4'

        vt1 = VariantTranscript('1', 123, 'C', 'G')
        vt1.refseq_transcript = 'NM_000051.5'
        preferred = self.tx_eff_ccds._get_preferred_refseq_transcript(vt0, vt1)

        self.assertEqual(preferred.refseq_transcript, 'NM_000051.5', "Latest version should be preferred")
    
    def test__get_preferred_refseq_transcript_non_refseq(self):
        vt0 = VariantTranscript('1', 123, 'C', 'G')
        vt0.refseq_transcript = 'NP_000000.0'

        vt1 = VariantTranscript('1', 123, 'C', 'G')
        vt1.refseq_transcript = 'NM_000051.5'

        self.assertRaises(ValueError, self.tx_eff_ccds._get_preferred_refseq_transcript, vt0, vt1)
