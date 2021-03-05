class EnsemblDbImport:
    """
    We are going to run the following command to retrieve ENST to RefSeq mappings:

    echo "use homo_sapiens_core_75_37; SELECT transcript.stable_id, transcript.version, xref.display_label
    FROM transcript, object_xref, xref,external_db WHERE transcript.transcript_id = object_xref.ensembl_id
    AND object_xref.ensembl_object_type = 'Transcript' AND object_xref.xref_id = xref.xref_id AND
    xref.external_db_id = external_db.external_db_id AND external_db.db_name = 'RefSeq_mRNA'" |
    mysql -u anonymous -h ensembldb.ensembl.org > enst_refseq_map.tsv

    I was not able to get this to work using INTO OUTFILE due to anonymous user login.

    Tab-delimited output file looks like:

    stable_id	version	display_label
    ENST00000517143	1	NR_046932.1
    ENST00000362897	1	NR_046944.1
    """
    def __init__(self, filename):
        self.filename = filename
        self.enst_refseq = self._ens_parse()

    @staticmethod
    def _find_smaller_enst(val1, val2):
        """
        Grab the smaller integer NM value.
        :return:
        """
        a = int(val1.split('_')[1].split('.')[0])
        b = int(val2.split('_')[1].split('.')[0])

        if a < b:
            return val1
        elif b < a:
            return val2
        elif a == b:
            return val1
        else:
            raise Exception(f"Please check the mapping file for {val1} and {val2}.")

    def _ens_parse(self):
        """
        This file contains multiple enst ids related to multiple refseq (many-to-many).
        If there is a conflict, we will simply pick the lowest NM refseq id.
        mapping = {<ENST>: <RefSeq>}
        :return:
        """
        mapping = {}
        with open(self.filename, 'r') as myfile:
            # Skip header.
            next(myfile)
            for line in myfile:
                line = line.rstrip('\n').split('\t')
                enst = '.'.join([line[0], line[1]])
                refseq = line[2]
                if enst not in mapping:
                    mapping[enst] = refseq
                else:
                    mapping[enst] = self._find_smaller_enst(mapping[enst], refseq)
        return mapping
