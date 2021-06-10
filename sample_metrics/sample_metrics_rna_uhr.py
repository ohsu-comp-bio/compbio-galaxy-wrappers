#!/usr/bin/env python

# DESCRIPTION: Create sample level metrics to be passed to the CGD for RNA.  Metrics
#  are passed as a json dump.
# usage: sample_metrics.py -h


import argparse
import json
import numpy


VERSION = '1.0.2'


def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('counts', help='Kallisto counts matrix.')
    parser.add_argument('hk_file', help='List of housekeeping genes to output')
    parser.add_argument('outjson', help='json output for CGD')
    parser.add_argument('out', help='Output file in human readable text format.')
    parser.add_argument('--gene_filt', required=False, help='Further filter the output Kall metrics by a gene list.')
    parser.add_argument('--by_gene', action="store_true", help='Collect metrics by gene.')
    parser.add_argument('--uhr_zero', default=None, help='List of genes that should be zero in a UHR sample')
    parser.add_argument('--uhr_low', default=None, help='List of genes that should be low in a UHR sample')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class KallistoCounts:
    """
    header:
    target_id	length	eff_length	est_counts	tpm
    """
    def __init__(self, filename, gene_filt=None):
        self.filename = filename
        self.gene_filt = gene_filt
        self.header = self._get_header()
        self.recs = self._gather_recs()
        self.agg_recs = self._agg_recs()
        if gene_filt:
            self.recs = self._gene_filt(self.recs, gene_filt)
            self.agg_recs = self._gene_filt(self.agg_recs, gene_filt)

    @staticmethod
    def _gene_filt(recs, genes):
        """
        Based on an input list of gene ids, filter the counts down further.
        :return:
        """
        new_recs = {}
        for entry in recs:
            if entry in genes:
                new_recs[entry] = recs[entry]
        return new_recs

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

    def _agg_recs(self):
        """
        Collect TPM and raw metrics together per-gene.
        :return:
        """
        agg_recs = {}
        for val in self.recs.values():
            if val.hgnc not in agg_recs:
                agg_recs[val.hgnc] = KallistoGeneCounts(val.hgnc)
            agg_recs[val.hgnc].incr(float(val.raw), float(val.tpm))
        return agg_recs


class KallistoGeneCounts:
    def __init__(self, hgnc):
        self.hgnc = hgnc
        self.tpm = 0.0
        self.raw = 0.0

    def incr(self, raw, tpm):
        """
        Increment tpm and raw when values come in.
        This can probably use __add__ or something like that.
        :return:
        """
        self.raw += raw
        self.tpm += tpm


class KallistoMetrics(KallistoCounts):
    """
    Metrics based on Kallisto counts structures.
    """
    def __init__(self, filename, hk, gene_filt=None, agg=False):
        super(KallistoMetrics, self).__init__(filename, gene_filt)
        self.agg = agg
        if agg:
            self.recs = self.agg_recs
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
        if self.agg:
            return {x.hgnc: [x.hgnc, x.raw, x.tpm] for x in self.recs.values() if x.hgnc in self.hk}
        else:
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


class Genes:
    """
    Load genes from a single-column file.
    """
    def __init__(self, filename):
        self.filename = filename
        self.genes = self._get_genes()

    def _get_genes(self):
        """
        Put the list of genes from file in to a list.
        """
        genes = []
        with open(self.filename, 'r') as fh:
            for lines in fh:
                geneid = lines.rstrip('\n').split('\t')
                genes.append(geneid[0])
        return genes

class UHRGenes:
    """
    Load list of genes and thresholds from a coverage file
    """

    def __init__(self, filename):
        self.filename = filename
        self.uhr_genes = self._get_uhr_genes()
        self.thresholds = self._get_uhr_thresholds()
        
    def _get_uhr_genes(self):
        """
        Get the gene names from the coverage file
        """
        uhr_genes = []
        with open(self.filename, 'r') as uhr_file:
            for line in uhr_file:
                geneid = line.rstrip('\n').split('|')
                uhr_genes.append(geneid[4])
        return uhr_genes
    
    def _get_uhr_thresholds(self):
        """
        Get the thresholds for acceptable uhr coverage
        """
        
        uhr_thresholds = {}
        with open(self.filename, 'r') as uhr_file:
            for line in uhr_file:
                gene_line = line.rstrip('\n').split('\t')
                geneid = gene_line[0].split('|')
                if len(gene_line) > 1:
                    uhr_thresholds[geneid[4]] = float(gene_line[1])
                else:
                    uhr_thresholds[geneid[4]] = 0
        return uhr_thresholds


class MetricsWriter:
    """
    Write the KallistoMetrics to json/tsv outputs.
    """
    def __init__(self, metrics):
        self.metrics = metrics

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

    def write_txt(self, outfile):
        """
        Put in to human-readable format.
        :param outfile:
        :return:
        """
        with open(outfile, 'w') as out:
            # COUNTS section
            out.write('COUNTS_greaterorequalto_Limit\tnumber_of_transcripts\n')
            for k in sorted(self.metrics.summ_raw):
                out.write('%s\t%s\n' % (k, self.metrics.summ_raw[k]))
            out.write('\n')

            # TPM section
            out.write('TPM_greaterorequalto_Limit\tnumber_of_transcripts\n')
            for k in sorted(self.metrics.summ_tpms):
                out.write('%s\t%s\n' % (k, self.metrics.summ_tpms[k]))
            out.write('\n')

            # HOUSEKEEPING section
            out.write('HOUSEKEEPING GENES\n')
            if self.metrics.agg:
                out.write('hgvs_id\tcounts\ttpm\n')
                # sort by HUGO gene name, because life is hard
                for val in sorted(self.metrics.hk_counts.values()):
                    out.write('%s\t%s\t%s\n' % (val[0], val[1], val[2]))
                out.write('\n')
            else:
                out.write('hgvs_tx_id\tensembl_id\tcounts\ttpm\n')
                # sort by HUGO gene name, because life is hard
                for val in sorted(self.metrics.hk_counts.values()):
                    out.write('%s\t%s\t%s\t%s\n' % (val[0], val[1], val[2], val[3]))
                out.write('\n')


def main():
    args = supply_args()
    hk = Genes(args.hk_file).genes
    if args.gene_filt:
        fgenes = Genes(args.gene_filt).genes
    else:
        fgenes = None
    if args.by_gene:
        metrics = KallistoMetrics(args.counts, hk, fgenes, agg=True)
    else:
        metrics = KallistoMetrics(args.counts, hk, fgenes, agg=False)
        
    if args.uhr_zero != None:
        uhr_zero = UHRGenes(args.uhr_zero)
        uhr_zero_genes = uhr_zero.uhr_genes
        zero_metrics = KallistoMetrics(args.counts, uhr_zero_genes)
        zero_counts = zero_metrics.hk_counts
        print(zero_counts)
        for geneid in zero_counts:
            if zero_counts[geneid][2] != '0':
                raise ValueError(geneid + ' has a raw count of ' + zero_counts[geneid][2] + '. For a UHR sample all genes in the uhr_zero file should have a raw count of 0. If this is not a UHR sample do not provide a uhr_zero file.')
    if args.uhr_low != None:
        uhr_low = UHRGenes(args.uhr_low)
        uhr_low_genes = uhr_low.uhr_genes
        low_metrics = KallistoMetrics(args.counts, uhr_low_genes)
        low_counts = low_metrics.hk_counts
        print(low_counts)
        for geneid in low_counts:
            if float(low_counts[geneid][2]) > uhr_low.thresholds[geneid]:
                raise ValueError(geneid + ' has a raw count of ' + low_counts[geneid][2] + ' but a UHR threshold of ' + str(uhr_low.thresholds[geneid]) + '. For a UHR sample all genes in the uhr_low file should have a raw count less than their threshold. If this is not a UHR sample do not provide a uhr_low file.')
    writer = MetricsWriter(metrics)
    writer.write_json(args.outjson, metrics.summ_tpms, metrics.summ_raw)
    writer.write_txt(args.out)


if __name__ == "__main__":
    main()