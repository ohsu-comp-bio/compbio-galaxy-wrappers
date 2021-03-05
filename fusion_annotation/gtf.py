class Gtf:
    def __init__(self, filename, txs=None):
        self.filename = filename
        self.txs = txs
        if self.txs:
            self.raw = self._gather_txs()
        else:
            self.raw = None

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
        :return:
        """
        gather = []
        with open(self.filename) as myfile:
            for rec in myfile:
                if not rec.startswith('#'):
                    info = self._info_to_dict(rec)
                    tx_id = info['transcript_id']
                    rec = rec.rstrip('\n').split('\t')
                    this_rec = {'seqname': rec[0],
                                'source': rec[1],
                                'feature': rec[2],
                                'start': rec[3],
                                'end': rec[4],
                                'score': rec[5],
                                'strand': rec[6],
                                'frame': rec[7],
                                'info': info}
                    if tx_id in self.txs:
                        gather.append(this_rec)
        return gather

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

    @staticmethod
    def exon_to_coord(exons, first=True):
        """
        Match coordinates to exon numbers.  Rules:

        Only pass exon entries for a particular tx.
        Must be an ONLY_REF_SPLICE entry.
        If this is the first gene (left):
            if plus:
                get end coord
            if minus:
                get start coord+1
        If it is the second gene (right):
            if plus:
                get start coord+1
            if minus:
                get end coord
        :return:
        """
        exon_coord = {}
        for exon in exons:
            start = exon['start']
            end = exon['end']
            strand = exon['strand']
            exon_no = exon['info']['exon_number']

            if first:
                if strand == '+':
                    exon_coord[end] = exon_no
                elif strand == '-':
                    exon_coord[start] = exon_no
                else:
                    raise ValueError("Strand is not - or + in _exon_to_coord.")
            else:
                if strand == '+':
                    exon_coord[start] = exon_no
                elif strand == '-':
                    exon_coord[end] = exon_no
                else:
                    raise ValueError("Strand is not - or + in _exon_to_coord.")
        return exon_coord
