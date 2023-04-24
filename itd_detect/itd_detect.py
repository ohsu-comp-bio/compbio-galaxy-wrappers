#!/usr/bin/env python

# Process FASTQ files to obtain potential ITD calls from given genomic regions.
# Should be used on your assay-specific data.  Generally, you will pre-process your FASTQ files in any way
# you'd like, then run PEAR to merge reads (if using paired-end).  Then, run BWA and output a coordinate sorted
# BAM file.  This is an appropriate input file for the detection algorithm.
# Paired mode is not working very well especially when estimating VAFs, use at own risk.
# VAFs are also off in paired-end mode because I am only looking for exact matches during overlap determination.
# 0.4.2 - If we can't get HGVS results, just return an output file where they are blank.
# 1.0.0 - Remove hgvs from script, create vcf as output.
# 1.0.2 - Added bcor and fgfr targets
# 1.1.0 - Fixed error when read sequences are empty due to amplicon clipping.  Provide parameter to limit output size
# 1.1.1 - Added BCOR and FGFR regions to support STP4 workflows
# of potential ITDs.

import argparse
import pysam
from collections import defaultdict
from copy import deepcopy
from itertools import groupby
from operator import itemgetter

VERSION = '1.1.1'

def supply_args():
    """
    Populate args.
    https://docs.python.org/2.7/library/argparse.html
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--sam', required=True, help='Input SAM File')
    parser.add_argument('--outfile', required=True, help='Output TSV')
    parser.add_argument('--outfile_vcf', required=True, help='Output VCF')
    parser.add_argument('--ref', required=True, help='Input Reference Sequence')
    parser.add_argument('--sample_id', required=True, help='Sample ID string to be added to the VCF.')
    parser.add_argument('--ref_build', choices=['hg19'], default='hg19', help='Which reference build to utilize')
    parser.add_argument('--target', default='flt3_e14', help='Region to target')
    parser.add_argument('--coords', help='Coordinate range, in the format [chrom:start-stop], 1-based.')
    parser.add_argument('--paired', action='store_true', help='Data is paired-end data.')
    parser.add_argument('--chr_prefix', action='store_true', help='Add chr prefix to all chromosome identifiers.')
    parser.add_argument('--min_size', type=int, default=12, help='Minimum size of ITD call to write to VCF.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()

    if args.target and args.coords:
        raise Exception("You can't specify both a known target region and a custom coordinate at the same time.")
    if not args.target and not args.coords:
        raise Exception("You must specify either the target or coords option.")

    return args


class SuffixArray:
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


class GetSeq:
    """
    Retrieve the sequence we need to search through, utilizing the pysam package.
    """
    def __init__(self, filename, region, chr_prefix):
        self.chr_prefix = chr_prefix
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
        if self.chr_prefix:
            known_targ = {'flt3': ('chr13', 28577411, 28682904),
                          'flt3_e13': ('chr13', 28608438, 28608544),
                          'flt3_e14': ('chr13', 28608219, 28608351),
                          'flt3_e15': ('chr13', 28608024, 28608128),
                          'bcor_e15': ('chrX', 39910499, 39911653),
                          'fgfr1_e10': ('chr8', 38275746, 38275891),
                          'fgfr1_e18': ('chr8', 38268656, 38271322),
                          'kmt2a': ('chr11', 118307205, 118397539)}
        else:
            known_targ = {'flt3': ('13', 28577411, 28682904),
                          'flt3_e13': ('13', 28608438, 28608544),
                          'flt3_e14': ('13', 28608219, 28608351),
                          'flt3_e15': ('13', 28608024, 28608128),
                          'bcor_e15': ('X', 39910499, 39911653),
                          'fgfr1_e10': ('8', 38275746, 38275891),
                          'fgfr1_e18': ('8', 38268656, 38271322),
                          'kmt2a': ('11', 118307205, 118397539)}
        return known_targ[target]

    def _get_ref_dups(self):
        flt3_dups = {}
        if max(list(map(len, self.curr_long))) > 0:
            for k, v in self.curr_long.items():
                diff = v[1] - v[0]
                if k not in flt3_dups:
                    flt3_dups[k] = diff
                else:
                    raise Exception(k + " already in flt3_dups, please check.")
        return flt3_dups


class Sequence:
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


class SequenceCollection:
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

        self.samfile.close()

        self.handle_out = open(outfile, 'w')
        self.handle_out.write('\t'.join(['CHROM', 'START', 'STOP', 'ITD LENGTH', 'ITD COUNT', 'REF COUNT', 'VAF ESTIMATE', 'SEQUENCE']))
        self.handle_out.write('\n')

        if not paired:
            self.sample_data = self._spray_coords(False)
        else:
            self.sample_data = self._spray_coords(True)

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
            if not read.is_proper_pair or read.is_secondary or read.is_supplementary or not read.seq:
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
            if line.seq:
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
        sample_data = []
        for k, v in self.itd_list.items():
            in_ref_cnt = 0
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

            maxseq = list(self.max_pos[k].keys())[0]
            entry = self.max_pos[k][maxseq]
            diff = entry.curr_long[maxseq][1] - entry.curr_long[maxseq][0]
            if len(entry.gnmic_seq) - entry.curr_long[maxseq][1] >= len(maxseq):
                chrom = self.refseq.coords[0]
                start_coord = entry.pos + entry.curr_long[maxseq][0] + 1 - entry.soft_clip
                stop_coord = start_coord + diff - 1
                this_seq = self.refseq.fasta.fetch(reference=chrom, start=start_coord - 1, end=stop_coord)
                # for the VCF, need to get the base that precedes the insertion
                padding_base = self.refseq.fasta.fetch(reference=chrom, start=start_coord - 2, end=start_coord - 1)
                itd_cnt = len(self.itd_list[diff])
                if paired:
                    ref_cnt = (self._sam_seq_count(self._ref_seq_create(stop_coord), self._dup_junc_seq_create(this_seq))) / 2
                else:
                    ref_cnt = self._sam_seq_count(self._ref_seq_create(stop_coord), self._dup_junc_seq_create(this_seq))
                # TODO: Parameter
                if len(set(self.itd_list[diff])) > 10 or max([len(x) for x in self.itd_list[diff]]) > 30:
                    to_write = [chrom, str(start_coord), str(stop_coord), str(diff), str(itd_cnt), str(ref_cnt),
                                    "{0:0.3f}".format(self._calc_vaf(ref_cnt, itd_cnt)), this_seq]
                    self.handle_out.write('\t'.join(to_write))
                    self.handle_out.write('\n')
                    sample_data.append(ItdCallVcf(chrom, str(start_coord), str(stop_coord), str(diff), itd_cnt,
                                               ref_cnt, "{0:0.3f}".format(self._calc_vaf(ref_cnt, itd_cnt)),
                                               this_seq, padding_base))
        return sample_data

    def _seq_diffs(self, seq):
        """
        Create the structure containing number of times distances between duplicate sequences occur.
        :return:
        """
        if max(list(map(len, seq.curr_long))) > 8:
            for k, v in seq.curr_long.items():
                diff = v[1] - v[0]
                # TODO: parameter
                if len(k) > 8:
                    if diff not in self.itd_list:
                        self.itd_list[diff] = [k]
                    else:
                        self.itd_list[diff].append(k)

                    if diff not in self.max_pos:
                        self.max_pos[diff] = {k: seq}
                    else:
                        if len(k) > len(list(self.max_pos[diff].keys())[0]):
                            self.max_pos[diff] = {k: seq}


class ItdCall:
    def __init__(self, chrom, start, stop, diff, itd_cnt, ref_cnt, vaf, seq, pad_seq):
        """
        We need fields for the VCF.
        CHROM: chrom
        POS: self.start - 1
        ID: .
        REF: pad_seq
        ALT: pad_seq+seq
        QUAL: .
        FILTER: .
        INFO: diff can go here (ITD_LENGTH=)
        FORMAT:
            GT: pretty much always 0/1
            DP: itd_cnt + ref_cnt
            AD: ref_cnt, itd_cnt
            AF: vaf
        samples
        :param chrom:
        :param start:
        :param stop:
        :param diff:
        :param itd_cnt:
        :param ref_cnt:
        :param vaf:
        :param seq:
        :param pad_seq:
        """
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.diff = diff
        self.itd_cnt = itd_cnt
        self.ref_cnt = ref_cnt
        self.vaf = vaf
        self.seq = seq
        self.pad_seq = pad_seq


class ItdCallVcf(ItdCall):
    def __init__(self, chrom, start, stop, diff, itd_cnt, ref_cnt, vaf, seq, pad_seq):
        super(ItdCallVcf, self).__init__(chrom, start, stop, diff, itd_cnt, ref_cnt, vaf, seq, pad_seq)
        self.pos = str(int(self.start) - 1)
        self.id = '.'
        self.ref = self.pad_seq
        self.alt = self.pad_seq + self.seq
        self.qual = '.'
        self.filter = '.'
        self.info = 'ITD_LENGTH=' + diff
        self.frmt = 'GT:DP:AD:AF'
        self.gt = '0/1'
        self.dp = str(int(self.itd_cnt + self.ref_cnt))
        self.ad = ','.join([str(int(self.ref_cnt)), str(int(self.itd_cnt))])
        self.sample = ':'.join([self.gt, self.dp, self.ad, self.vaf])
        self.to_write = '\t'.join([self.chrom, self.pos, self.id, self.ref, self.alt,
                         self.qual, self.filter, self.info, self.frmt, self.sample])


class VcfWrite:
    def __init__(self, filename, sample_id, sample_data, min_size):
        self.filename = filename
        self.sample_id = sample_id
        self.sample_data = sample_data
        self.min_size = min_size

    def write_me(self):
        with open(self.filename, 'w') as self.myfile:
            self._write_header()
            for vrnt in self.sample_data:
                if len(vrnt.alt) > self.min_size:
                    self.myfile.write(vrnt.to_write)
                    self.myfile.write('\n')
        self.myfile.close()

    def _write_header(self):
        """
        Write the header to the VCF.
        :return:
        """
        self.myfile.write("##fileformat=VCFv4.2\n")
        self.myfile.write("##INFO=<ID=ITD_LENGTH,Number=1,Type=Integer,Description=\"ITD estimated length\">\n")
        self.myfile.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        self.myfile.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"k-mer Depth\">\n")
        self.myfile.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"k-mer depth supporting reference/indel at the site\">\n")
        self.myfile.write("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Variant allele frequency\">\n")
        self.myfile.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', self.sample_id]))
        self.myfile.write('\n')


def main():

    # http://stackoverflow.com/questions/13560037/effcient-way-to-find-longest-duplicate-string-for-python-from-programming-pearl
    # Original implementation horribly slow.  If I can find my original sa code, will replace with this.

    args = supply_args()
    my_seq = GetSeq(args.ref, args.target, args.chr_prefix)
    if args.paired:
        coll = SequenceCollection(args.sam, args.outfile, my_seq, True)
    else:
        coll = SequenceCollection(args.sam, args.outfile, my_seq, False)
    sample_data = coll.sample_data

    VcfWrite(args.outfile_vcf, args.sample_id, sample_data, args.min_size).write_me()


if __name__ == "__main__":
    main()
