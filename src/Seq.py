#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @author: Tian Chen
# @times: 2020/10/13  15:37

import TmDeltaG
import re
import os
import sys
import itertools
import logging
from functools import reduce
import numpy as np
import globalVar as gl


class Sequence:
    """ the attributes of sequence."""

    def __init__(self, chr, start, seq, Prange, Flag, snp_dic):
        # self.chr = chr
        self.chr = chr.split('_')[0].lower()   # for methylation, 20220319

        self.start = start
        self.strand = Flag
        self.len = len(seq)
        self.end = start + self.len - 1
        if self.strand == '+':
            self.seq = seq
        else:
            self.seq = rev_com(seq)
        self.gc = gc(self.seq)
        self.poly = Poly(self.seq)
        self.new_seq, self.snp, self.Maf = self.SNP_new(snp_dic)  # return new_seq, snp info list
        self.cg, self.new_seq, self.tm = self._ct_cg_cnt() if "_CT" in chr else self._ga_cg_cnt()  # for methylation, 20220324

    def _ct_cg_cnt(self):
        cg_cnt, tm, new_seq = 0, 0, ''
        seqs = self.seq.upper()
        # print(f"before merge base: {self.strand}, {seqs}")
        for i in range(1, len(self.seq)):
            if seqs[i-1: i+1] == 'CG':
                cg_cnt += 1
        if self.strand == '+' and seqs[-1] == 'C':
            cg_cnt += 1
            seqs = seqs[:-1] + 'Y'
        elif self.strand == '-' and seqs[0] == 'G':
            cg_cnt += 1
            seqs = 'R' + seqs[1:]
        if self.new_seq != '':
            new_seq = seqs.replace('CG', 'YG') if self.strand == '+' else seqs.replace('CG', 'CR')
            tm_seq = new_seq.replace('Y', 'T') if self.strand == '+' else new_seq.replace('R', 'A')
            tm = TmDeltaG.Cal(tm_seq.upper(), complement(tm_seq.upper())).Tm
        if new_seq != '' and ('Y' in new_seq[-3:] or 'R' in new_seq[-3:]):
            # if new_seq != '' and (new_seq[-1] == 'Y' or new_seq[-1] == 'R'):
            new_seq = ''
        # print(f"after merge base: {self.strand}, {seqs}")
        return cg_cnt, new_seq, tm

    def _ga_cg_cnt(self):
        cg_cnt, tm, new_seq = 0, 0, ''
        seqs = self.seq.upper()
        # print(f"before merge base: {self.strand}, {seqs}")
        for i in range(1, len(self.seq)):
            if seqs[i - 1: i + 1] == 'CG':
                cg_cnt += 1
        if self.strand == '+' and seqs[0] == 'G':
            cg_cnt += 1
            seqs = 'R' + seqs[1:]
        elif self.strand == '-' and seqs[-1] == 'C':
            cg_cnt += 1
            seqs = seqs[:-1] + 'Y'
        if self.new_seq != '':  # 2022-08-05
            new_seq = seqs.replace('CG', 'CR') if self.strand == '+' else seqs.replace('CG', 'YG')
            tm_seq = new_seq.replace('R', 'A') if self.strand == '+' else new_seq.replace('Y', 'T')
            tm = TmDeltaG.Cal(tm_seq.upper(), complement(tm_seq.upper())).Tm
        if new_seq != '' and ('Y' in new_seq[-3:] or 'R' in new_seq[-3:]):
            # if new_seq != '' and (new_seq[-1] == 'Y' or new_seq[-1] == 'R'):
            new_seq = ''
        # print(f"after merge base: {self.strand}, {seqs}")
        return cg_cnt, new_seq, tm

    def SNP(self, snp_dic):
        """
        Function: after find the snp postion, use new base replace alternative base.
        """
        snp = []
        if self.strand == '-':
            new_seq = list(rev_com(self.seq))
        else:
            new_seq = list(self.seq)
        lens = 0
        # consider 3' can't appear snp position.
        for i in range(self.start, self.start + self.len):

            if self.chr + ':' + str(i) in snp_dic.keys():
                if '-' in snp_dic[self.chr + ':' + str(i)]:  # avoid indel position
                    return '', snp
                snp.append(snp_dic[self.chr + ':' + str(i)])
                base = ("".join(snp_dic[self.chr + ':' + str(i)][3:5])).replace(',', '')
                try:
                    new_seq[lens] = merge_base[base]
                except:
                    # print('appear indel')
                    pass
                # 3' 5bp can't appear snp.

                if (self.strand == '+' and self.len - lens <= 5) or (self.strand == '-' and lens <= 5):
                    return '', snp
            lens += 1
        return ''.join(new_seq), snp

    def SNP_new(self, snp_dic):
        """
        Caculate SNP Maf function: no degeneracy base just accumulate MAF.
        :param snp_dic:
        :return: new_seq, snp_info, Maf
        """
        snp = []
        maf = []
        Maf = 0
        leng = 0
        # consider 3' can't appear snp position.
        for i in range(self.start, self.start + self.len):
            if self.chr.lower() + ':' + str(i) in snp_dic.keys():
                # avoid indel position, 3' 10-mer can't contain snp position.
                snp.append(snp_dic[self.chr + ':' + str(i)])
                maf.append(max([float(x) for x in snp_dic[self.chr + ':' + str(i)][-1].split(',')]))
                if ('-' in snp_dic[self.chr + ':' + str(i)]) or (self.strand == '+' and self.len - leng <= 10) or (
                        self.strand == '-' and leng <= 10):
                    return '', snp, Maf
            leng += 1
        if len(maf) >= 1:
            if max(maf) > 0.1:
                return '', snp, Maf
            else:
                for i in range(2, len(maf) + 1):
                    for arr in list(itertools.combinations(maf, i)):
                        Maf += reduce(lambda x, y: x * y, arr)
        return self.seq, snp, Maf


def complement(seq):
    '''Return the complement sequence'''
    return ''.join([ambiguous_dna_complement[base] for base in seq])


def reverse(seq):
    """Reverses a string given to it."""
    return seq[::-1]


def rev_com(seq):
    """Return the reverse_complement sequence"""
    com = complement(seq)
    rev_com = reverse(com)
    return rev_com


def gc(seq):
    g_count = (seq.upper()).count('G')
    c_count = (seq.upper()).count('C')
    gc = (float(g_count + c_count) / len(seq)) * 100

    return round(gc, 2)


def Poly(seq):
    seq = seq.upper()
    seq_3 = seq[-10:]
    # add 1bp base for Methylation, 20220324
    if re.search('AAAAAAA', seq, flags=0) or re.search('AAAAAA', seq_3, flags=0):
        return 1
    if re.search('TTTTTTT', seq, flags=0) or re.search('TTTTTT', seq_3, flags=0):
        return 1
    if re.search('GGGGGG', seq, flags=0) or re.search('GGGGG', seq_3, flags=0):
        return 1
    if re.search('CCCCCC', seq, flags=0) or re.search('CCCCC', seq_3, flags=0):
        return 1
    return 0


def hits(script_dir, path, chr_name, start, lens, strand):
    """
    Function Hits: find Aln hits from hits file
        chr1:1234:   18| 19|  (1234 is the right pos)
        chr2:-1234:  18| 19|  (1234 is the left pos)
    """
    if lens < 18:
        lens = 18
    if strand == '+':
        name = chr_name + ':' + str(start + lens - 1)
    else:
        name = chr_name + ':-' + str(start)
    try:
        cmd = f"{script_dir}/catPos {path}/Hits_result.ln.gz {gl.hit_idx[name][0]} " \
            f"{gl.hit_idx[name][1]}  {lens} {gl.hit_idx[name][lens - 17]} {gl.hit_idx[name][lens - 16]}"
    except:
        logging.info(f"can't find {chr_name}:{start}:{lens} specificity value, default to 19999.")
        return 19999
    read_obj = os.popen(cmd)
    aln_info = read_obj.readlines()
    read_obj.close()
    if not aln_info:
        # sys.stderr.write(f"can't find {chr_name}:{start}:{lens} specificity value, default to -1\n")
        logging.info(f"{cmd}")
        logging.info(f"can't find {chr_name}:{start}:{lens} specificity value, default to -1")
        # sys.exit(1)
        return -1

    return int(aln_info[0].split('|')[0])


def self_hairpin(seq):
    """
    检测引物自身是否存在成发夹结构
    :param seq: 5'-3’ primer seq
    :return: True | False
    """
    string2 = complement(seq[::-1])
    len1 = len(seq)
    res = np.zeros((len1 + 1, len1 + 1))
    result = 0
    for i in range(1, len1 + 1):
        for j in range(1, len1 + 1):
            if string2[i - 1] == seq[j - 1]:
                res[i, j] = res[i - 1, j - 1] + 1
    result = set([int(j) for i in res for j in i if j >= 4])  # stem length >=4
    for i in list(result)[::-1]:
        idx = np.argwhere(res == i)
        if len(idx) == 1:  # avoid a part of palindromic(回文序列) in 3'
            if i > 5 and idx[0, -1] == len1:
                # print("exists self dimer in 3'")
                return True
            else:
                continue
        stem = []  # seq
        pos = []  # seq start position
        for x, y in idx:
            Rckmer = string2[x - i:x]  # get the rev_comp seq
            stem.append(Rckmer)
            pos.append(x)
            if rev_com(Rckmer) in stem:
                pos_value = pos[stem.index(rev_com(Rckmer))]
                dis = x - i - pos_value
                if dis in range(3, 9) and pos_value - i <= 10:  # loop length in [3,8]bp, hairpin not in 3' 10bp range
                    deltG = TmDeltaG.calDeltaG(Rckmer, complement(Rckmer), mono_conc=50, diva_conc=3, dntp_conc=0.8,
                                               deltaH=None, deltaS=None)
                    if deltG < -3:
                        # print("exists hairpin")
                        return True
    return False


'''
def Score(self):

    return self.hits
'''

ambiguous_dna_complement = {
    "a": "t",
    "c": "g",
    "g": "c",
    "t": "a",
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "m": "k",
    "R": "Y",
    "r": "y",
    "W": "W",
    "w": "w",
    "S": "S",
    "s": "s",
    "Y": "R",
    "y": "r",
    "K": "M",
    "k": "m",
    "V": "B",
    "v": "b",
    "H": "D",
    "h": "d",
    "D": "H",
    "d": "h",
    "B": "V",
    "b": "v",
    "X": "X",
    "x": "x",
    "N": "N",
    "n": "n",
    "-": "-",
}

merge_base = {
    'AG': 'R', 'GA': 'R',
    'CT': 'Y', 'TC': 'Y',
    'AC': 'M', 'CA': 'M',
    'GT': 'K', 'TG': 'K',
    'GC': 'S', 'CG': 'S',
    'AT': 'W', 'TA': 'W',
    'ATC': 'H', 'TAC': 'H',
    'TCA': 'H', 'CAT': 'H',
    'CTA': 'H', 'ACT': 'H',
    'GTC': 'B', 'GCT': 'B',
    'TGC': 'B', 'CGT': 'B',
    'TCG': 'B', 'CTG': 'B',
    'GAC': 'V', 'GCA': 'V',
    'AGC': 'V', 'CGA': 'V',
    'ACG': 'V', 'CAG': 'V',
    'GAT': 'D', 'GTA': 'D',
    'AGT': 'D', 'TGA': 'D',
    'ATG': 'D', 'TAG': 'D',
    'ATCG': 'N', 'ATGC': 'N',
    'ACTG': 'N', 'AGTC': 'N',
    'AGCT': 'N', 'ACGT': 'N',
    'TACG': 'N', 'TAGC': 'N',
    'TCAG': 'N', 'TGAC': 'N',
    'TCGA': 'N', 'TGCA': 'N',
    'CATG': 'N', 'CAGT': 'N',
    'CTAG': 'N', 'CTGA': 'N',
    'CGAT': 'N', 'CGTA': 'N',
    'GATC': 'N', 'GACT': 'N',
    'GTCA': 'N', 'GCTA': 'N',
    'GTAC': 'N', 'GCAT': 'N',
}

if __name__ == "__main__":
    # primer = Seq(chr='chr1',start=44462,seq='GAAGTTAGGTTCatcaatccgctccaa',Prange=250, Flag=0)
    # print(primer.rev_com)
    # print(primer.len)
    # print(primer.tm)
    # print(primer.poly)
    #  print(primer.hits)
    strs1 = 'TAAAAGAATCTAGAAATGATATTTCAC'
    strs2 = 'GTCACCAGGACGTCCTGG'
    seq = 'GAAAAAGGCGTTTTGCCT'
    print(self_hairpin(strs2.upper()))
    print(self_hairpin(seq.upper()))
    print(gc(seq))
    print(rev_com(seq))