#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @author: Tian Chen
# @times: 2020/9/12  17:21
import numpy as np
import pandas as pd
# from multiprocessing import Process, Pool, Manager, Lock
from functools import lru_cache

import Seq
# from PrepareData import time_this
import globalVar as gl


def build_snGSP_tagGSP(amp_min, amp_max, f_primer, r_primer):
    """构建上下游引物集合"""
    # manager = Manager()
    # res_dic = manager.dict()
    # P = Pool(40)
    res = []
    for f_id, f_row in enumerate(f_primer):
        for r_id, r_row in enumerate(r_primer):
            k = r_row[2] - f_row[1] + 1
            if k in range(amp_min, amp_max + 1):
                res.append([f_id, r_id, f_row[0], f_row[1], f_row[2], f_row[4], r_row[1], r_row[2],
                            r_row[4], f_row[-1] + r_row[-1]])
                # res.append([f_id, r_id, f_row[0], f_row[1], f_row[2], f_row[4], r_row[1], r_row[2],
                #             r_row[4], f_row[6] + r_row[6]])
            # P.apply_async(cut_dt, args=(f_id, r_id, f_row, r_row, opt.amp_max, opt.amp_min, res_dic, ))
    # P.close()
    # P.join()
    combine = pd.DataFrame(data=res,
                           columns=['f_idx', 'r_idx', 'chr', 'Fstart', 'Fend', 'Fseq', 'Rstart', 'Rend', 'Rseq',
                                    'score'])
    # combine['score'] = combine['score'].astype(int)
    # print('before sort: ', combine.head())
    combine.sort_values(by='score', ascending=True, inplace=True)
    # combine.to_excel("snGSP_tagGSP.xlsx", index=False)
    # print('after sort: ', combine.head(100))
    return combine.values


def cut_dt(f_id, r_id, f_row, r_row, amp_max, amp_min, res_dic):
    """选择合适扩增子的引物对"""
    k = r_row[2] - f_row[1] + 1
    if k in range(amp_min, amp_max + 1):
        res_dic[f'{f_id}_{r_id}'] = [f_id, r_id, f_row[0], f_row[1], f_row[2], f_row[4], r_row[1], r_row[2],
                                     r_row[4], f_row[6] + r_row[6]]
        # return res
    # return None


def GSP_snGSP_tagGSP(dis_ud_min, dis_ud_max, amp_max, snp_flag, snGSP2tar, path,
                ref, name, count, f_primer, r_primer, snGSP_tagGSP, flag, tar_start, tar_end,
                res_lst, used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst):
    """选择合适的GSP, snGSP, tagGSP"""
    #step_size = 10
    dt_cnt = len(snGSP_tagGSP)
    step = int(dt_cnt/2000)
    step_size = step if step != 0 else 1
    print(f"step_size: {step_size}")
    for i in range(0, len(snGSP_tagGSP), step_size):
        dt = snGSP_tagGSP[i]
        if flag % 2 == 1:
            udt = select_udt(dis_ud_min, dis_ud_max, amp_max, snp_flag, snGSP2tar,
                             res_lst, used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst,
                             f_primer, r_primer, dt[0], dt[1], chr_name=dt[2], d_start=dt[3], d_end=dt[4],
                             t_start=dt[6], t_end=dt[7], t_type='TR', u_type='UF', tar_start=tar_start, tar_end=tar_end)

        else:
            udt = select_udt(dis_ud_min, dis_ud_max, amp_max, snp_flag, snGSP2tar,
                             res_lst, used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst,
                             f_primer, r_primer, dt[0], dt[1], chr_name=dt[2], d_start=dt[6], d_end=dt[7],
                             t_start=dt[3], t_end=dt[4], t_type='TF', u_type='UR', tar_start=tar_start, tar_end=tar_end)
        # print(udt)
        if not udt: continue
        new_start = udt[2][1] if flag % 2 == 1 else udt[1][1]
        amp_len = dt[7] - dt[3] + 1
        amp_gc = Seq.gc(ref[dt[3] - 1:dt[7]])
        name_pattern = lambda s: str(name) + s + str(count)
        name_id, GSP_name, snGSP_name, tagGSP_name = name_pattern('_'), name_pattern('-GSP-'),\
                                                name_pattern('-snGSP-'), name_pattern('-tagGSP-')

        primer = [name_id, dt[2], amp_gc, amp_len,
                  GSP_name, udt[0][4], udt[0][1], udt[0][2], udt[0][5], udt[0][6],
                  snGSP_name, 'GAACGACATGGCTACGATCCGACTT' + udt[1][4], udt[1][1], udt[1][2], udt[1][5], udt[1][6],
                  tagGSP_name, 'CGCTTGGCCTCCGACTT' + udt[2][4], udt[2][1], udt[2][2], udt[2][5], udt[2][6],
                  udt[0][-1], udt[1][-1], udt[2][-1]]
        # gl.lock.acquire()
        with open(path+'/standard.txt', 'a+') as f:
            f.write(','.join([str(x) for x in primer]))
            f.write('\n')
        res_lst.append(primer)
        store_used(udt, used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst, flag)
        # gl.lock.release()
        return 'ok', flag, new_start
    return 'not ok', flag, 'None'

# @time_this
def select_udt(dis_ud_min, dis_ud_max, amp_max, snp_flag, snGSP2tar,
               res_lst, used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst,
               f_primer, r_primer, f_id, r_id,
               chr_name, d_start, d_end, t_start, t_end, t_type, u_type, tar_start, tar_end):
    """选择GSP-snGSP-tagGSP"""
    # detect snGSP and tagGSP overlap in history primers.
    if over_check(amp_max, chr_name, d_start, d_end, "D",
                  used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst)\
            or over_check(amp_max, chr_name, t_start, t_end, t_type, used_d_lst, used_uf_lst, used_tr_lst,
                           used_ur_lst, used_tf_lst):
        print(f"meet conflict with history primers")
        return None
    snGSP, tagGSP = f_primer[f_id], r_primer[r_id]
    if u_type == 'UR':
        snGSP, tagGSP = tagGSP, snGSP
    # u_d_lst = [x[5] for x in res_lst] + [x[11][21:] for x in res_lst]
    u_d_lst = [x[5] for x in res_lst]
    tagGSP_dimer_flag = dimer_check(tagGSP[4], u_d_lst) if len(res_lst) > 0 else False #
    if snp_flag and snGSP2tar_check(d_start, d_end, tar_start, tar_end, u_type, snGSP2tar):
        return None

    GSP, u_flag = find_GSP(amp_max, dis_ud_min, dis_ud_max, f_primer, r_primer, d_start, d_end, u_type,
                           res_lst, used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst)
    if not u_flag or tagGSP_dimer_flag:
        return None
    return [GSP, snGSP, tagGSP]


def over_check(amp_max, chr_name, start, end, types, used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst):
    """检测序列与已存在的引物序列，位置上是否存在冲突"""
    if len(used_d_lst) == 0:
        return False
    if 'F' in types:
        ur_lst = [x for x in used_ur_lst if x[0] == chr_name]
        tr_lst = [x for x in used_tr_lst if x[0] == chr_name]
        # print("types:",types)
        if len(ur_lst) > 0 or len(tr_lst) > 0:
            if types == 'UF':
                flag_ur = pos_overlap(amp_max, start, end, ur_lst, types, diff_name_id=False)
                flag_tr = pos_overlap(amp_max, start, end, tr_lst, types, diff_name_id=True)
                print(flag_ur, flag_tr)
                return flag_ur or flag_tr
            elif types == 'TF':
                flag_tr = pos_overlap(amp_max, start, end, tr_lst, types, diff_name_id=False)
                flag_ur = pos_overlap(amp_max, start, end, ur_lst, types, diff_name_id=True)
                return flag_ur or flag_tr
    elif 'R' in types:
        uf_lst = [x for x in used_uf_lst if x[0] == chr_name]
        tf_lst = [x for x in used_tf_lst if x[0] == chr_name]
        # print("types:",types)
        if len(uf_lst) > 0 or len(tf_lst) > 0:
            if types == 'UR':
                flag_uf = pos_overlap(amp_max, start, end, uf_lst, types, diff_name_id=False)
                flag_tf = pos_overlap(amp_max, start, end, tf_lst, types, diff_name_id=True)
                return flag_uf or flag_tf
            elif types == 'TR':
                flag_tf = pos_overlap(amp_max, start, end, tf_lst, types, diff_name_id=False)
                flag_uf = pos_overlap(amp_max, start, end, uf_lst, types, diff_name_id=True)
                return flag_uf or flag_tf
    elif types == 'D':
        d_lst = [x for x in used_d_lst if x[0] == chr_name]
        tr_lst = [x for x in used_tr_lst if x[0] == chr_name]

        flag = pos_overlap(amp_max, start, end, d_lst, types, diff_name_id=False)
        flag_tr = pos_overlap(amp_max, start, end, tr_lst, types+'F', diff_name_id=True)
        return flag or flag_tr
    return False


def pos_overlap(amp_max, start, end, small_lst, types, diff_name_id=False):

    # print(types, small_lst, diff_name_id)
    for row in small_lst:
        un_amp = start - row[2] in range(-amp_max, -4) if 'F' in types else start - row[2] in range(-4, amp_max)
        #print("unexpected amp: ", row, start, end, un_amp)
        if (row[1] >= start and end >= row[1] + 4) or (row[1] <= start and start + 4 <= row[2]):
            return True
        elif diff_name_id and un_amp:
            return True
    return False


def snGSP2tar_check(d_start, d_end, tar_start, tar_end, u_type, distance):
    """检测目标区域与目标区域的距离是否满足要求"""
    # DF: tar_end - snGSP_start + 1, DR: snGSP_end - tar_start + 1
    dis = tar_end - d_start + 1 if u_type == 'UF' else d_end - tar_start + 1
    if dis <= distance:
        return False
    return True


def find_GSP(amp_max, dis_ud_min, dis_ud_max, f_primers, r_primer, start, end, u_type,
             res_lst, used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst):
    """设计GSP"""
    primers = f_primers if u_type == 'UF' else r_primer
    ut_lst = []
    if len(res_lst) > 0:
        for x in res_lst:
            ut_lst.extend([x[5], x[17][17:]])
    for row in primers:
        dis = start - row[2] + 1 if u_type == 'UF' else row[1] - end + 1
        # print(f"UD dis: {dis}")
        if dis_ud_min <= dis <= dis_ud_max\
                and not over_check(amp_max, row[0], row[1], row[2], u_type, used_d_lst, used_uf_lst, used_tr_lst,
                                   used_ur_lst, used_tf_lst)\
                and not dimer_check(row[4], ut_lst):
            return row, True
    return None, False


def dimer_check(seq, seq_lst):
    """依次遍历新添序列与已存在引物序列的dimer状况"""
    seq1 = Seq.rev_com(seq).upper()
    for seq2 in seq_lst:
        if seq_lcs_dimer(seq1, seq2.upper()):
            #print(f"Dimer: {seq}, {seq2}")
            return True
    return False


@lru_cache(maxsize=128)
def seq_lcs_dimer(string1, string2):
    """检测序列3'端的dimer"""
    len1, len2 = len(string1), len(string2)
    res = np.zeros((len2 + 1, len1 + 1))
    for i in range(1, len2 + 1):  # len2 + 1
        for j in range(1, len1 + 1):  # len1 + 1
            # print(i, j)
            if string2[i - 1] != string1[j - 1] and string2[i - 1] not in ['Y', 'R'] and string1[j - 1] not in ['Y', 'R']:
                pass
            elif string2[i - 1] == string1[j - 1] \
                    or (string2[i - 1] in ['Y', 'T', 'C'] and string1[j - 1] in ['Y', 'T', 'C']) \
                    or (string2[i - 1] in ['R', 'G', 'A'] and string1[j - 1] in ['R', 'G', 'A']):
                res[i, j] = res[i - 1, j - 1] + 1
    res = np.tril(res)
    # for a in res:
    #     print("".join([str(_) for _ in a]))

    result = int(np.max(res))
    # print(result)
    idx = np.argwhere(res == result)
    # print(idx)
    if result >= 7:
        for x, y in idx:
            #print(x, y, len1, len2)
            if (x == len2 and (0 not in [res[x - result + 1 + i, i + 1] for i in range(result)])):
                    # (y == len1 and (0 not in [res[i + 1, -result + i] for i in range(result)])):
                #print('dimer: ', string1, string2)
                return True
    return False


def store_used(udt, used_d_lst, used_uf_lst, used_tr_lst, used_ur_lst, used_tf_lst, flag):
    """存储设计好引物的位置信息"""
    used_d_lst.append([udt[0][0], udt[1][1], udt[1][2], udt[1][4]])
    if flag % 2 == 1:
        used_uf_lst.append([udt[0][0], udt[0][1], udt[0][2], udt[0][4]])
        used_tr_lst.append([udt[0][0], udt[2][1], udt[2][2], udt[2][4]])
    else:
        used_ur_lst.append([udt[0][0], udt[0][1], udt[0][2], udt[1][4]])
        used_tf_lst.append([udt[0][0], udt[2][1], udt[2][2], udt[2][4]])


if __name__ == '__main__':
    # seq1 = 'GTGTTACAAACCAGTGTGTG'
    # seq2 = ['AAATCAATAACAGAACACACACACT']
    # seq1 = 'GTTAAGTCATAAGTTAGTTAAGTCATAAGTTAGT'
    # seq2 = ['TATGACTTAACTAACTTATGACTTAACTAACTTAT']
    seq1 = 'GCTAGATAGATTGAATCAACCAGTTAAT'
    seq2 = ['CCCTTGTTGTTCATTAGTATAGTATTAACT']
    res = dimer_check(seq1, seq2)
    print(res)
