class Gtf:
    def __init__(self, filename, txs=None):
        self.filename = filename
        self.txs = txs
        self.raw = self._gather_txs()

    @staticmethod
    def _info_to_dict(line):
        """
        Provide the information in the final GTF column as a dictionary.

        gene_id "ENSG00000223972.4"; transcript_id "ENSG00000223972.4"; gene_type "pseudogene"; gene_status "KNOWN";
        gene_name "DDX11L1"; transcript_type "pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1";
        level 2; havana_gene "OTTHUMG00000000961.2";
        :return:
        """
        info_dict = {}
        info = line.rstrip('\n').split('\t')[8].split(';')
        for entry in info:
            if entry:
                entry = entry.strip()
                key = entry.split(' ')[0].strip()
                val = entry.split(' ')[1].strip().strip('\"')
                info_dict[key] = val
        return info_dict

    @staticmethod
    def filt_feature(choices, feat='exon'):
        """
        Given a specific GTF feature, return entries that match.
        :param choices:
        :param feat:
        :return:
        """
        output = []
        for choice in choices:
            if choice['feature'] == feat:
                output.append(choice)
        return output

    @staticmethod
    def filt_tx(choices, tx):
        """
        Given a "transcript_id", return matching entries.
        :return:
        """
        output = []
        for choice in choices:
            tx_id = choice['info']['transcript_id']
            if tx_id == tx:
                output.append(choice)
        return output

    def _gather_txs(self):
        """
        Find records corresponding to the transcripts we care about.
        {TX_ID: {FEAT: [rec, rec, ...]}}
        :return:
        """
        gtf_recs = {}
        with open(self.filename) as myfile:
            for rec in myfile:
                if not rec.startswith('#'):
                    info = self._info_to_dict(rec)
                    tx_id = info['transcript_id']
                    if self.txs:
                        if tx_id in self.txs:
                            rec = rec.rstrip('\n').split('\t')
                            feat = rec[2]
                            this_rec = {'seqname': rec[0],
                                        'source': rec[1],
                                        'feature': feat,
                                        'start': rec[3],
                                        'end': rec[4],
                                        'score': rec[5],
                                        'strand': rec[6],
                                        'frame': rec[7],
                                        'info': info}

                            if tx_id not in gtf_recs:
                                gtf_recs[tx_id] = {}
                            if feat not in gtf_recs[tx_id]:
                                gtf_recs[tx_id][feat] = []

                            gtf_recs[tx_id][feat].append(this_rec)

                    else:
                        raise NotImplementedError("Functionality to provide records for complete GTF not currently "
                                                  "implemented.")
        return gtf_recs

    @staticmethod
    def exon_to_ccds(txs):
        """
        Get a CCDS ID, if it exists.
        :param txs:
        :return:
        """
        for tx in txs:
            if 'ccdsid' in tx['info']:
                return tx['info']['ccdsid']

    def exon_to_coord(self, chrom, coord, tx_id, trm='exon'):
        if trm in self.raw[tx_id]:
            for exon in self.raw[tx_id][trm]:
                exon_chrom = exon['seqname']
                start = int(exon['start'])
                end = int(exon['end'])
                exon_no = exon['info']['exon_number']

                if chrom == exon_chrom:
                    if int(coord) >= start and int(coord) <= end:
                        return exon_no
