#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Time  : 2019/8/21 15:42
# @Author: Chentian
# @File  : Dimer.py

from numpy.core.multiarray import ndarray

import Seq
import re
import math
import numpy as np
import ThermodynamicsParameters as TP
import re
import globalVar as gl

gl.init()


def primer_dimer(primer_f, primer_r, pd_full, bind_string, pd_extensible=gl.pd_extensible):
    F_len = len(primer_f)
    R_len = len(primer_r)

    rating_dict = []

    max_len = F_len if F_len > R_len else R_len
    score = []
    rcompr = (Seq.complement(primer_r))[::-1]
    rcomprlc = rcompr.lower()

    fprimer_r = primer_f[::-1].lower()

    rprimer_r = primer_r[::-1]

    primer_matrix = np.zeros((F_len, 6))  # type: ndarray is binary

    # build a g c t index according Fprimer reverse seq.
    for i in range(F_len):
        if fprimer_r[i] == 'a':
            primer_matrix[i, 0] = 1
        elif fprimer_r[i] == 'g':
            primer_matrix[i, 1] = 1
        elif fprimer_r[i] == 'c':
            primer_matrix[i, 2] = 1
        elif fprimer_r[i] == 't':
            primer_matrix[i, 3] = 1
        elif fprimer_r[i] == 'y':
            primer_matrix[i, 4] = 1
            # primer_matrix[i, 2], primer_matrix[i, 3] = 1, 1
        elif fprimer_r[i] == 'r':
            primer_matrix[i, 5] = 1
            # primer_matrix[i, 0], primer_matrix[i, 1] = 1, 0

    primer_comp = np.zeros((R_len, F_len))

    # build matrix for reverse Fprimer and reverse complement Rprimer.
    for i in range(R_len):
        if rcomprlc[i] == 'a':
        # if rcomprlc[i] is 'a' or rcomprlc[i] is 'r':
            primer_comp[i] = primer_matrix[:, 0]
        # elif rcomprlc[i] is 'g' or rcomprlc[i] is 'r':
        elif rcomprlc[i] == 'g':
            primer_comp[i] = primer_matrix[:, 1]
        elif rcomprlc[i] == 'c':
        # elif rcomprlc[i] is 'c' or rcomprlc[i] is 'y':
            primer_comp[i] = primer_matrix[:, 2]
        elif rcomprlc[i] == 't':
            primer_comp[i] = primer_matrix[:, 3]
        elif rcomprlc[i] == 'y':
            primer_comp[i] = primer_matrix[:, 4]
        elif rcomprlc[i] == 'r':
            primer_comp[i] = primer_matrix[:, 5]

    pd_len = F_len + R_len - 1 if pd_full else max_len - 2
    primer_comp1 = primer_comp.astype(np.int32)
    # print(primer_comp1)
    row, col = primer_comp1.shape
    for k in range(pd_len):
        # score[k] = 0
        score_p = 0

        bind = ''
        start = F_len - 1 if k > F_len - 1 else k
        end = R_len - 1 if k > R_len - 1 else k

        # print("extend and pd_full=0:")
        if pd_extensible and not pd_full:
            if abs(start - k) > col:
                continue
            if primer_comp1[0][start] == 0:
                continue
            if primer_comp1[end][start - k] == 0:
                continue
        for i in range(R_len):
            if k - i < F_len:
                if k - i >= 0:
                    bind = bind + (str(int(primer_comp[i][k - i])))

            else:
                bind = bind + "2"

        """ Single matched based surrounded by mismatches are unstable,
        so we remove them with the regexp (look ahead is needed otherwise strings of consecutive
         match/mismatches are not caught)"""

        bind = re.sub('01(?=[^1])', '00', bind, count=0)

        if not re.search('[1]', bind): continue

        # Find start and end of similarity
        pb_init = None
        pb_end = None
        # print(bind)
        for i in range(len(bind)):
            if bind[i] == '1':
                if pb_init == None:
                    pb_init = i
                pb_end = i
        if pb_init != None:
            for i in range(pb_init, pb_end):
                if bind[i:i + 2] == '00': continue
                if bind[i] == '2': continue
                seq_key = primer_f[F_len - k + i - 1:F_len - k + i + 1] + rprimer_r[i:i + 2]
                # print(seq_key)
                if 'Y' in seq_key or 'R' in seq_key:
                    if seq_key.replace('Y', 'C').replace('R', 'G') in gl.oligo_dG.keys():
                        seq_key = seq_key.replace('Y', 'C').replace('R', 'G')
                    elif seq_key.replace('Y', 'T').replace('R', 'G') in gl.oligo_dG.keys():
                        seq_key = seq_key.replace('Y', 'T').replace('R', 'G')
                    elif seq_key.replace('Y', 'C').replace('R', 'A') in gl.oligo_dG.keys():
                        seq_key = seq_key.replace('Y', 'C').replace('R', 'A')
                    elif seq_key.replace('Y', 'T').replace('R', 'A') in gl.oligo_dG.keys():
                        seq_key = seq_key.replace('Y', 'T').replace('R', 'A')
                    # print("new_seq: ", seq_key)
                # score_p += gl.oligo_dG[primer_f[F_len - k + i - 1:F_len - k + i + 1] + rprimer_r[i:i + 2]]
                score_p += gl.oligo_dG[seq_key]


            initterm = 'init' + rprimer_r[pb_init:pb_init + 1]
            initterm = initterm.replace('Y', 'C').replace('R', 'G')
            score_p += TP.oligo_dG[initterm]

            endterm = 'init' + rprimer_r[pb_end:pb_end + 1]
            endterm = endterm.replace('R', 'G').replace('Y', 'C')

            score_p += TP.oligo_dG[endterm]

            # score[k] = round(score_p,2)
            score.append(round(score_p, 2))
            bind_string.append(bind)
            rating_dict.append([round(score_p, 2), k])
            # rating_dict[round(score_p,2)] = k
            # print(round(score_p,2),k)
    # print("score len:",len(score))
    # print("ratind_dict len:",len(rating_dict))
    #     print(bind)
    if len(score) == 0:
        return None, 0, rating_dict
    else:
        return min(score), score.index(min(score)), rating_dict


def recalculate_dG(temperature, mg_conc=gl.mg_conc, dntp_conc=gl.dntp_conc,
                   monovalent_cation_conc=gl.monovalent_cation_conc):
    if mg_conc > dntp_conc:
        salt_correction = math.sqrt(mg_conc - dntp_conc)
    else:
        salt_correction = 0

    na_eq = (monovalent_cation_conc + 120 * salt_correction) / 1000.0

    # print monovalent_cation_conc, na_eq

    entropy_adjust = 0.368 * math.log(na_eq, math.e)

    for key in TP.dH_full:
        if re.match('init', key): continue
        dS = TP.dS_full[key] + entropy_adjust
        dG = TP.dH_full[key] - ((273.15 + temperature) * (dS / 1000))
        gl.oligo_dG[key] = dG


def draw_dimer(primer_f, primer_r, pos, index, bind_string):
    rprimer_r = primer_r[::-1]
    dimer_binding = ''
    pr_space = ''
    fspace = ''
    rspace = ''

    F_len = len(primer_f)
    R_len = len(rprimer_r)
    max_len = F_len if F_len > R_len else R_len
    fspace_def = max_len - F_len if max_len - F_len > 0 else 0
    fspace_def2 = pos - max_len + 1 if pos > max_len - 1 else 0
    fspace = " " * (fspace_def + fspace_def2)
    fspace2 = pos - max_len + 1 if pos > max_len - 1 else 0
    # print("max_len:",max_len)
    if pos + 1 >= F_len:
        rspace = " " * (max_len - pos - 1)
    else:
        rspace = fspace

    pr_space = " " * (F_len - pos - 1)

    # print("bind_sting:",len(bind_string),pos,bind_string[index])
    # print(len(bind_string[index]),F_len)
    for i in range(pos + 1):
        # print(i, F_len)
        if i >= R_len:
            break
        if bind_string[index][i] == '1':
            dimer_binding = dimer_binding + "|"
        elif bind_string[index][i] == '0':
            dimer_binding = dimer_binding + "."
        else:
            dimer_binding = dimer_binding + " "
        # print(i, pos, bind_string[index])

    # print (fspace+"5' "+primer_f+" 3'\n"+rspace+"   "+pr_space+dimer_binding+"\n"+rspace+pr_space+"3' "+rprimer_r+" 5'")
    return (
                fspace + "5' " + primer_f + " 3'\n" + rspace + "   " + pr_space + dimer_binding + "\n" + rspace + pr_space + "3' " + rprimer_r + " 5'\n")


if __name__ == '__main__':
    pd_extensible = 1
    mg_conc = 3
    oligo_conc = 50
    monovalent_cation_conc = 50
    dntp_conc = 0.8
    pd_temperature = 37
    oligo_dG = {}
    DG_dict = []
    recalculate_dG(37)
    bind_string = []
    seq1 = 'GCCATGCGAAATGCTGG'  # 5'-3'
    # seq2 = 'CTAACTGAAATTTCTCTCAAACCTCACAG'
    # seq1 = Seq.rev_com('GCCATGCGAAATGCTGG')
    seq2 = 'AGCTAAAATGTAATCACCAGCATT'  # 5'- 3'
    seq1 = 'GTGTAGTGGGAACGTTTGCGTT'
    seq2 = 'GTGTAGTGGGAACGTTTGCGTT'
    # DG, index, DG_dict = primer_dimer('AATCATTGCCCAGACCAAGC','ATTGGGGCTTGGAATACTCTTGGA',0,bind_string)
    # print(DG,index,DG_dict[index][1])
    # draw_dimer('AATCATTGCCCAGACCAAGC','ATTGGGGCTTGGAATACTCTTGGA',DG_dict[index][1], index, bind_string)
    DG, index, DG_dict = primer_dimer(seq1, seq2, 0, bind_string)
    print(DG, index)
    # print(bind_string)
    print(DG_dict[index])
    result = draw_dimer(seq1, seq2, DG_dict[index][1], index, bind_string)
    print(result)
