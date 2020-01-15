#!/usr/bin/env python

# Process FASTQ files to obtain potential ITD calls from given genomic regions.
# Should be used on your assay-specific data.  Generally, you will pre-process your FASTQ files in any way
# you'd like, then run PEAR to merge reads (if using paired-end).  Then, run BWA and output a coordinate sorted
# BAM file.  This is an appropriate input file for the detection algorithm.
# Paired mode is not working very well especially when estimating VAFs, use at own risk.
# VAFs are also off in paired-end mode because I am only looking for exact matches during overlap determination.

import argparse
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.parser
import pysam
from collections import defaultdict
from copy import deepcopy
from itertools import groupby
from itertools import imap
from operator import itemgetter
from string import Template

VERSION = '0.4.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sam', required=True, help='Input SAM File')
    parser.add_argument('--outfile', required=True, help='Output File')
    parser.add_argument('--ref', required=True, help='Input Reference Sequence')
    parser.add_argument('--ref_build', choices=['hg19'], default='hg19', help='Which reference build to utilize')
    parser.add_argument('--target', choices=['flt3', 'flt3_e13', 'flt3_e14', 'flt3_e15'], default='flt3_e14', help='Region to target')
    parser.add_argument('--coords', help='Coordinate range, in the format [chrom:start-stop], 1-based.')
    parser.add_argument('--paired', action='store_true', help='Data is paired-end data.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    if args.target and args.coords:
        raise Exception("You can't specify both a known target region and a custom coordinate at the same time.")
    if not args.target and not args.coords:
        raise Exception("You must specify either the target or coords option.")

    return args


class SuffixArray(object):
    """Analyze all common strings in the text.

    Short substrings of the length _step a are first pre-sorted. The are the
    results repeatedly merged so that the garanteed number of compared
    characters bytes is doubled in every iteration until all substrings are
    sorted exactly.

    Arguments:
        text:  The text to be analyzed.
        _step: Is only for optimization and testing. It is the optimal length
               of substrings used for initial pre-sorting. The bigger value is
               faster if there is enough memory. Memory requirements are
               approximately (estimate for 32 bit Python 3.3):
                   len(text) * (29 + (_size + 20 if _size > 2 else 0)) + 1MB

    Return value:      (tuple)
      (sa, rsa, lcp)
        sa:  Suffix array                  for i in range(1, size):
               assert text[sa[i-1]:] < text[sa[i]:]
        rsa: Reverse suffix array          for i in range(size):
               assert rsa[sa[i]] == i
        lcp: Longest common prefix         for i in range(1, size):
               assert text[sa[i-1]:sa[i-1]+lcp[i]] == text[sa[i]:sa[i]+lcp[i]]
               if sa[i-1] + lcp[i] < len(text):
                   assert text[sa[i-1] + lcp[i]] < text[sa[i] + lcp[i]]
    >>> suffix_array(text='banana')
    ([5, 3, 1, 0, 4, 2], [3, 2, 5, 1, 4, 0], [0, 1, 3, 0, 0, 2])

    Explanation: 'a' < 'ana' < 'anana' < 'banana' < 'na' < 'nana'
    The Longest Common String is 'ana': lcp[2] == 3 == len('ana')
    It is between  tx[sa[1]:] == 'ana' < 'anana' == tx[sa[2]:]
    """

    def __init__(self, text, _step=16):
        self.tx = text
        self.step = min(max(_step, 1), len(self.tx))
        self.size = len(self.tx)
        self.sa, self.rsa = self._create_sa()
        self.lcp = self._create_lcp()

    def _create_sa(self):
        sa = list(range(len(self.tx)))
        sa.sort(key=lambda i: self.tx[i:i + self.step])
        # a boolean map for iteration speedup.
        grpstart = self.size * [False] + [True]
        # It helps to skip yet resolved values. The last value True is a sentinel.
        rsa = self.size * [None]
        stgrp, igrp = '', 0

        for i, pos in enumerate(sa):
            st = self.tx[pos:pos + self.step]
            if st != stgrp:
                grpstart[igrp] = (igrp < i - 1)
                stgrp = st
                igrp = i
            rsa[pos] = igrp
            sa[i] = pos
        grpstart[igrp] = (igrp < self.size - 1 or self.size == 0)

        while grpstart.index(True) < self.size:
            # assert step <= size
            nextgr = grpstart.index(True)
            while nextgr < self.size:
                igrp = nextgr
                nextgr = grpstart.index(True, igrp + 1)
                glist = []
                for ig in range(igrp, nextgr):
                    pos = sa[ig]
                    if rsa[pos] != igrp:
                        break
                    newgr = rsa[pos + self.step] if pos + self.step < self.size else -1
                    glist.append((newgr, pos))
                glist.sort()
                for ig, g in groupby(glist, key=itemgetter(0)):
                    g = [x[1] for x in g]
                    sa[igrp:igrp + len(g)] = g
                    grpstart[igrp] = (len(g) > 1)
                    for pos in g:
                        rsa[pos] = igrp
                    igrp += len(g)
            self.step *= 2
        del grpstart

        return sa, rsa

    def _create_lcp(self):
        lcp = self.size * [None]
        h = 0
        for i in range(self.size):
            if self.rsa[i] > 0:
                j = self.sa[self.rsa[i] - 1]
                while i != self.size - h and j != self.size - h and self.tx[i + h] == self.tx[j + h]:
                    h += 1
                lcp[self.rsa[i]] = h
                if h > 0:
                    h -= 1
        if self.size > 0:
            lcp[0] = 0

        return lcp


    def longest_common_substring(self):
        """
        """
        result = {}
        for i in range(1, len(self.tx)):
            if self.lcp[i] >= 8:
                j1, j2, h = self.sa[i - 1], self.sa[i], self.lcp[i]
                assert self.tx[j1:j1 + h] == self.tx[j2:j2 + h]
                substring = self.tx[j1:j1 + h]
                if not substring in result:
                    result[substring] = [j1]
                result[substring].append(j2)
        return dict((k, sorted(v)) for k, v in result.items())


    def longest_common_substring_sample(self):
        """
        """
        result = {}
        maxlen = max(self.lcp)
        for i in range(1, len(self.tx)):
            if self.lcp[i] == maxlen:
                j1, j2, h = self.sa[i - 1], self.sa[i], self.lcp[i]
                assert self.tx[j1:j1 + h] == self.tx[j2:j2 + h]
                substring = self.tx[j1:j1 + h]
                if not substring in result:
                    result[substring] = [j1]
                result[substring].append(j2)
        return dict((k, sorted(v)) for k, v in result.items())


class GetSeq(object):
    """
    Retrieve the sequence we need to search through, utilizing the pysam package.
    """
    def __init__(self, filename, region):
        self.coords = self._known_targets(region)
        self.coords = self._buffer_targets()
        self.fasta = pysam.FastaFile(filename)
        self.my_seq = self.fasta.fetch(reference=self.coords[0], start=self.coords[1], end=self.coords[2])
        self.my_arr = SuffixArray(self.my_seq)
        self.curr_long = self.my_arr.longest_common_substring()
        self.flt3_dups = self._get_ref_dups()

    def _buffer_targets(self, buffer=50):
        """
        Provide an extra genomic coordinate buffer to the regions of interest.
        :return:
        """
        return (self.coords[0], self.coords[1]-buffer, self.coords[2]+buffer)

    def _known_targets(self, target):
        """
        Hold the coordinates for commonly interrogated targets.
        These have been pulled from the RefSeq GFF, and are 1-based.
        :return:
        """
        known_targ = {'flt3': ('13', 28577411, 28682904),
                    'flt3_e13': ('13', 28608438, 28608544),
                    'flt3_e14': ('13', 28608219, 28608351),
                    'flt3_e15': ('13', 28608024, 28608128)}
        return known_targ[target]

    def _get_ref_dups(self):
        flt3_dups = {}
        if max(imap(len, self.curr_long)) > 0:
            for k, v in self.curr_long.items():
                diff = v[1] - v[0]
                if k not in flt3_dups:
                    flt3_dups[k] = diff
                else:
                    raise Exception(k + " already in flt3_dups, please check.")
        return flt3_dups


class Sequence(object):
    """
    Dissecting pysam objects for necessary info, but will also include additional necessary structures on the sequence level.
    """
    def __init__(self, seq):
        self.gnmic_seq = seq.seq
        self.chrom = seq.rname
        self.cigar = seq.cigar
        self.soft_clip = self._soft_clip_amt()
        self.pos = seq.pos
        self.itd_arr = SuffixArray(self.gnmic_seq)
        self.curr_long = self.itd_arr.longest_common_substring_sample()
        self.dup_idx = self._dup_idx()

    def _dup_idx(self):
        dup_idx = {}
        for k, v in self.curr_long.items():
            diff = v[1] - v[0]
            dup_idx[k] = [self.pos + v[0], self.pos + v[1]]
        return dup_idx

    def _soft_clip_amt(self):
        """
        Determine total amount covered by soft clips, so that we can properly deduce the true start site.
        :return:
        """
        if self.cigar:
            if self.cigar[0][0] == 4:
                return self.cigar[0][1]
        return 0


class SequenceCollection(object):
    """
    Structures resulting from operations on multiple Sequence objects.
    """
    def __init__(self, filename, outfile, refseq, paired=False):
        self.filename = filename
        self.samfile = pysam.AlignmentFile(filename, 'rb')
        self.itd_list = {}
        self.max_pos = {}
        self.refseq = refseq
        self.flt3 = refseq.my_seq

        if not paired:
            for line in self.samfile.fetch(reference=refseq.coords[0], start=refseq.coords[1], end=refseq.coords[2]):
                self.this_seq = Sequence(line)
                self._seq_diffs(self.this_seq)
        else:
            for read1, read2 in self.read_pair_generator():
                if read1 and read2:
                    combine = self._combine_overlap_reads(read1.seq, read2.seq)
                    if combine:
                        if combine[1]:
                            # use r2 as template
                            new_read = self._record_modify(combine[0], read2)
                        else:
                            # use r1 as template
                            new_read = self._record_modify(combine[0], read1)
                        self.this_seq = Sequence(new_read)
                        self._seq_diffs(self.this_seq)
                    else:
                        pass
                        # self.this_seq = Sequence(read1)
                        # self._seq_diffs(self.this_seq)

        self.samfile.close()

        self.handle_out = open(outfile, 'w')
        self.handle_out.write('\t'.join(['CHROM', 'START', 'STOP', 'HGVS_G', 'HGVS_C', 'HGVS_P', 'ITD LENGTH', 'ITD COUNT', 'REF COUNT', 'VAF ESTIMATE', 'SEQUENCE']))
        self.handle_out.write('\n')

        if not paired:
            self._spray_coords(False)
        else:
            self._spray_coords(True)

        self.handle_out.close()

    def _record_modify(self, new_seq, rec):
        """
        Record should be modified to then be passed to a Sequence object.
        :return:
        """
        new_rec = deepcopy(rec)
        new_rec.seq = new_seq
        return new_rec

    def _combine_overlap_reads(self, r1, r2, o=20):
        """
        Look for the similar portion, then combine.
        :return:
        """
        if r1[:o] in r2:
            mtch_pnt = r2.index(r1[:o])
            return (r2[:mtch_pnt] + r1, True)
        elif r2[:o] in r1:
            mtch_pnt = r1.index(r2[:o])
            return (r1[:mtch_pnt] + r2, False)
        return None

    def read_pair_generator(self):
        """
        Generate read pairs in a BAM file or within a region string.
        Reads are added to read_dict until a pair is found.
        """
        read_dict = defaultdict(lambda: [None, None])
        for read in self.samfile.fetch(reference=self.refseq.coords[0], start=self.refseq.coords[1], end=self.refseq.coords[2]):
            if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
                continue
            qname = read.query_name
            if qname not in read_dict:
                if read.is_read1:
                    read_dict[qname][0] = read
                else:
                    read_dict[qname][1] = read
            else:
                if read.is_read1:
                    yield read, read_dict[qname][1]
                else:
                    yield read_dict[qname][0], read
                del read_dict[qname]

    def _sam_seq_count(self, seq_to_find, seq_to_excl):
        """
        Move through the BAM file and count records when they contain a sequence we are looking for.
        :return:
        """
        samfile = pysam.AlignmentFile(self.filename, 'rb')
        match_cnt = 0
        for line in samfile.fetch(reference=self.refseq.coords[0], start=self.refseq.coords[1], end=self.refseq.coords[2]):
            if seq_to_find in line.seq and seq_to_excl not in line.seq:
                match_cnt += 1
        samfile.close()
        return match_cnt

    def _ref_seq_create(self, pos, buffer=10):
        """
        Create sequence that represents both sides of a duplication breakpoint.
        :return:
        """
        return self.refseq.fasta.fetch(reference=self.refseq.coords[0], start=int(pos-buffer), end=int(pos+buffer))

    def _dup_junc_seq_create(self, seq):
        """
        Create the sequence that you would only find in a read containing the duplication.
        :return:
        """
        return seq[-10:] + seq[:10]

    def _calc_vaf(self, ref, alt):
        """
        Do a VAF calculation.
        :return:
        """
        return alt / (ref + alt + 0.0)

    def _spray_coords(self, paired):
        """
        Go through the max_pos dict and provide coords for everything based off of what's in here.
        max_pos looks like:
         5: {'TTTTATTTTA': <__main__.Sequence object at 0x2b6bde468b90>}
         itd_list looks like:
         35: ['GGAATGGAATGGAATGGAATGGAATGGAA', 'GGAATGGAATGGAATGGAATGGAATGGAA', 'GGAATGGAATGGAATGGAATGGAATGGAA', 'GGAATGGAATGGAATGGAATGGAATGGAA']
        TODO: This function is a disaster.
        :return:
        """

        for k, v in self.itd_list.items():
            in_ref_cnt = 0.0
            total_cnt = len(v)
            for seq in sorted(v, key=len, reverse=True):
                if seq in self.flt3:
                    in_ref_cnt += 1
            real = in_ref_cnt / total_cnt
            # TODO: parameter
            if real > 0.75:
                long_seq = max(v, key=len)
                maxlen = len(long_seq)
                itd_len_diff = k - maxlen
                flt3_pos = self.flt3.find(long_seq)
                itd_start = flt3_pos - itd_len_diff + self.refseq.coords[1]
                itd_stop = itd_start + k

            maxseq = self.max_pos[k].keys()[0]
            entry = self.max_pos[k][maxseq]
            diff = entry.curr_long[maxseq][1] - entry.curr_long[maxseq][0]
            if len(entry.gnmic_seq) - entry.curr_long[maxseq][1] >= len(maxseq):
                chrom = self.refseq.coords[0]
                start_coord = entry.pos + entry.curr_long[maxseq][0] + 1 - entry.soft_clip
                stop_coord = start_coord + diff - 1
                this_seq = self.refseq.fasta.fetch(reference=chrom, start=start_coord - 1, end=stop_coord)
                itd_cnt = len(self.itd_list[diff])
                if paired:
                    ref_cnt = (self._sam_seq_count(self._ref_seq_create(stop_coord), self._dup_junc_seq_create(this_seq))) / 2
                else:
                    ref_cnt = self._sam_seq_count(self._ref_seq_create(stop_coord), self._dup_junc_seq_create(this_seq))
                # TODO: Parameter
                if len(set(self.itd_list[diff])) > 10 or max([len(x) for x in self.itd_list[diff]]) > 30:
                    #and self._calc_vaf(ref_cnt, itd_cnt) > 0.005:
                    my_hgvs = HgvsVars(chrom, start_coord, stop_coord)
                    to_write = [chrom, str(start_coord), str(stop_coord),
                                my_hgvs.var_g, my_hgvs.hgvs_print(my_hgvs.var_c),
                                my_hgvs.hgvs_print(my_hgvs.var_p), str(diff),
                                str(itd_cnt), str(ref_cnt), "{0:0.3f}".format(self._calc_vaf(ref_cnt, itd_cnt)),
                                this_seq]
                    self.handle_out.write('\t'.join(to_write))
                    self.handle_out.write('\n')


    def _seq_diffs(self, seq):
        """
        Create the structure containing number of times distances between duplicate sequences occur.
        :return:
        """
        if max(imap(len, seq.curr_long)) > 8:
            for k, v in seq.curr_long.items():
                diff = v[1] - v[0]
                # TODO: parameter
                if len(k) > 8:
#                    v = [v[0] + seq.pos - seq.soft_clip, v[1] + seq.pos - seq.soft_clip]
#                    v = [v[0] + seq.pos, v[1] + seq.pos]
                    if diff not in self.itd_list:
                        self.itd_list[diff] = [k]
                    else:
                        self.itd_list[diff].append(k)

                    if diff not in self.max_pos:
                        self.max_pos[diff] = {k: seq}
                    else:
                        if len(k) > len(self.max_pos[diff].keys()[0]):
                            self.max_pos[diff] = {k: seq}


class HgvsVars(object):
    """
    Provide HGVS c./p./g.
    #'NC_000013.10:g.28608250_28608300dup'
    :param object:
    :return:
    """
    def __init__(self, chrom, start, stop, genome='GRCh37'):
        self.hdp = hgvs.dataproviders.uta.connect()
        self.vm = hgvs.assemblymapper.AssemblyMapper(self.hdp, assembly_name=genome)
        self.hp = hgvs.parser.Parser()
        self.chrom = self._chrom_map(chrom)
        self.start = start
        self.stop = stop
        self.var_g = self._create_hgvs_g()
        self.var_c = self._create_hgvs_c()
        self.var_p = self._create_hgvs_p()

    def _create_hgvs_g(self):
        """
        Provide the actual HGVS nomenclature.  To start, we are just providing the g. value.  To provide c./p., we will need
        to bring in the python hgvs package, which is probably worth doing.
        return cht_tmpl.substitute(self.xml_out)
        13:g.28608251_28608301dup
        :return:
        """
        hgvs_g = Template('${chrom}:g.${start}_${stop}dup')
        return hgvs_g.substitute(chrom=self.chrom, start=self.start, stop=self.stop)

    def _create_hgvs_c(self):
        """
        Create a list of the HGVS c. values.
        :return:
        """
        var_c = []
        var_g = self.hp.parse_hgvs_variant(self.var_g)
        self.tx_list = self.hdp.get_tx_for_region(str(var_g.ac), 'splign', str(var_g.posedit.pos.start),
                                        str(var_g.posedit.pos.end))
        for entry in self.tx_list:
            try:
                var_c.append(self.vm.g_to_c(var_g, str(entry[0])))
            except:
                pass
        return var_c

    def _create_hgvs_p(self):
        """
        Create a list of the HGVS p. values.
        :return:
        """
        var_p = []
        for entry in self.var_c:
            try:
                var_p.append(self.vm.c_to_p(entry))
            except:
                pass
        return var_p

    def hgvs_print(self, hgvs):
        """

        :return:
        """
        return ', '.join([str(x) for x in hgvs])

    def _chrom_map(self, chrom):
        """
        Mapping of current common to RefSeq chromsome ids.
        :return:
        """
        chrom_map = {'1': 'NC_000001.10', '2': 'NC_000002.11', '3': 'NC_000003.11', '4': 'NC_000004.11', '5': 'NC_000005.9',
             '6': 'NC_000006.11', '7': 'NC_000007.14', '8': 'NC_000008.10', '9': 'NC_000009.11', '10': 'NC_000010.10',
             '11': 'NC_000011.9', '12': 'NC_000012.11', '13': 'NC_000013.10', '14': 'NC_000014.8', '15': 'NC_000015.9',
             '16': 'NC_000016.9', '17': 'NC_000017.10', '18': 'NC_000018.9', '19': 'NC_000019.9', '20': 'NC_000020.10',
             '21': 'NC_000021.8', '22': 'NC_000022.10', 'X': 'NC_000023.10', 'Y': 'NC_000024.9', 'MT': 'NC_012920.1'}
        # In case someone is trying to utilize chr prefixed chromosomes.
        if chrom.startswith('chr'):
            chrom = chrom[3:]
        try:
            return chrom_map[chrom]
        except:
            return None


def main():

    # http://stackoverflow.com/questions/13560037/effcient-way-to-find-longest-duplicate-string-for-python-from-programming-pearl
    # Original implementation horribly slow.  If I can find my original sa code, will replace with this.

    args = supply_args()
    my_seq = GetSeq(args.ref, args.target)
    if args.paired:
        SequenceCollection(args.sam, args.outfile, my_seq, True)
    else:
        SequenceCollection(args.sam, args.outfile, my_seq, False)

if __name__ == "__main__":
    main()
