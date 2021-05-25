#!/usr/bin/env python

# DESCRIPTION: Create sample level metrics to be passed to the CGD for RNA.  Metrics
#  are passed as a json dump.
# usage: sample_metrics.py -h


import argparse
import json
import numpy


VERSION = '1.0.0'


def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('counts', help='Kallisto counts matrix.')
    parser.add_argument('hk_file', help='List of housekeeping genes to output')
    parser.add_argument('outjson', help='json output for CGD')
    parser.add_argument('out', help='Output file in human readable text format.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class KallistoCounts:
    """
    header:
    target_id	length	eff_length	est_counts	tpm
    """
    def __init__(self, filename):
        self.filename = filename
        self.header = self._get_header()
        self.recs = self._gather_recs()

    def _get_header(self):
        """
        Grab the header from the first line of the file.
        :return:
        """
        with open(self.filename, 'r') as myfile:
            return myfile.readline().rstrip().split('\t')

    def _gather_recs(self):
        """
        Put all the records in a structure.
        :return:
        """
        recs = {}
        with open(self.filename, 'r') as myfile:
            next(myfile)
            for line in myfile:
                sline = line.rstrip().split('\t')
                target_id = sline[0]
                hgnc_tx = target_id.split('|')[4]
                if target_id not in recs:
                    recs[hgnc_tx] = KallistoRec(sline)
                else:
                    raise Exception("DUP")
        return recs


class KallistoMetrics(KallistoCounts):
    """
    Metrics based on Kallisto counts structures.
    """
    def __init__(self, filename, hk):
        super(KallistoMetrics, self).__init__(filename)
        self.tpms = [x.tpm for x in self.recs.values()]
        self.raws = [x.raw for x in self.recs.values()]
        self.summ_tpms = self._summ_tpm()
        self.summ_raw = self._summ_raw()
        # For housekeeping gene stats.
        self.hk = hk
        self.hk_counts = self._hk_filt()

    def _hk_filt(self):
        """
        Return Kallisto records that intersect with hk_genes.
        :return:
        """
        return {x.hgnc_tx: [x.hgnc_tx, x.enst, x.raw, x.tpm] for x in self.recs.values() if x.hgnc_tx in self.hk}

    def _summ_tpm(self):
        """
        Put TPM values in to bins, as such:
        {'TPM_0': 198909, 'TPM_0.01': 100329, 'TPM_0.1': 67443, 'TPM_1': 25875, 'TPM_10': 3840, 'TPM_100': 468,
        'TPM_1000': 61}
        :return:
        """
        sumrz_tpm = {}
        np_tpm = numpy.asarray(self.tpms, dtype=numpy.float64)
        sumrz_tpm["TPM_0"] = len(np_tpm[np_tpm >= 0])
        for i in range(-2, 4):
            idx = pow(10, i)
            sumrz_tpm["TPM_" + str(idx)] = len(np_tpm[np_tpm >= idx])
        return sumrz_tpm

    def _summ_raw(self):
        """
        Similar to _summ_tpm but for raw COUNTS instead.
        :return:
        """
        sumrz_counts = {}
        np_counts = numpy.asarray(self.raws, dtype=numpy.float64)
        sumrz_counts["COUNT_0"] = len(np_counts[np_counts >= 0])
        for i in range(0, 6):
            idx = pow(10, i)
            sumrz_counts["COUNT_" + str(idx)] = len(np_counts[np_counts >= idx])
        return sumrz_counts


class KallistoRec:
    """
    Gather information related to a single Kallisto record.
    ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-002|DDX11L1|
    1657|processed_transcript|	1657	1389.41	105.16	0.143597
    """
    def __init__(self, rec):
        self.rec = rec
        self.enst = self.rec[0].split('|')[0]
        self.hgnc_tx = self.rec[0].split('|')[4]
        self.hgnc = self.rec[0].split('|')[5]
        self.tpm = self.rec[4]
        self.raw = self.rec[3]


class HousekeepingGenes:
    """
    Load housekeeping genes from a file.
    """
    def __init__(self, filename):
        self.filename = filename
        self.hk_genes = self._get_hkgenes()

    def _get_hkgenes(self):
        """
        Housekeeping genes, filter for those genes are in the target file
        """
        hk_genes = []
        with open(self.filename, 'r') as fh_hk:
            for lines in fh_hk:
                geneid = lines.rstrip('\n').split('\t')
                hk_genes.append(geneid[0])
        return hk_genes


class MetricsWriter:
    """
    Write the KallistoMetrics to json/tsv outputs.
    """
    @staticmethod
    def write_json(filename, *args):
        """
        Write metrics to json.
        :param filename:
        :return:
        """
        with open(filename, 'w') as jsonout:
            to_dump = {k: v for d in args for k, v in d.items()}
            jsonout.write(json.dumps(to_dump))

    @staticmethod
    def write_txt(outfile, sumrz_cnts, sumrz_tpm, hk_counts):
        """
        Put in to human-readable format.
        :param outfile:
        :param sumrz_cnts:
        :param sumrz_tpm:
        :param hk_counts:
        :return:
        """
        with open(outfile, 'w') as out:
            # COUNTS section
            out.write('COUNTS_greaterorequalto_Limit\tnumber_of_transcripts\n')
            for k in sorted(sumrz_cnts):
                out.write('%s\t%s\n' % (k, sumrz_cnts[k]))
            out.write('\n')

            # TPM section
            out.write('TPM_greaterorequalto_Limit\tnumber_of_transcripts\n')
            for k in sorted(sumrz_tpm):
                out.write('%s\t%s\n' % (k, sumrz_tpm[k]))
            out.write('\n')

            # HOUSEKEEPING section
            out.write('HOUSEKEEPING GENES\n')
            out.write('hgvs_tx_id\tensembl_id\tcounts\ttpm\n')
            # sort by HUGO gene name, because life is hard
            for val in sorted(hk_counts.values()):
                out.write('%s\t%s\t%s\t%s\n' % (val[0], val[1], val[2], val[3]))
            out.write('\n')


def main():
    args = supply_args()
    hk = HousekeepingGenes(args.hk_file).hk_genes
    metrics = KallistoMetrics(args.counts, hk)
    writer = MetricsWriter()
    writer.write_json(args.outjson, metrics.summ_tpms, metrics.summ_raw)
    writer.write_txt(args.out, metrics.summ_raw, metrics.summ_tpms, metrics.hk_counts)


if __name__ == "__main__":
    main()
