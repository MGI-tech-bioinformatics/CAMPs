#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @author: Tian Chen
# @times: 2020/9/8  18:04
import pandas as pd


def init_lock(l, lock_2):
    global lock, lock_fail
    lock = l
    lock_fail =  lock_2


def init():
    global report_dic, hit_idx, failed
    global thread
    report_dic = {}
    hit_idx = {}
    thread = 2
    failed = []

    global used_tf, used_tr, used_uf, used_ur, used_d, std
    used_tf = pd.DataFrame(columns=['chr', 'tagGSP start', 'tagGSP end', 'tagGSP seq'])
    used_tr = pd.DataFrame(columns=['chr', 'tagGSP start', 'tagGSP end', 'tagGSP seq'])
    used_uf = pd.DataFrame(columns=['chr', 'GSP start', 'GSP end', 'GSP seq'])
    used_ur = pd.DataFrame(columns=['chr', 'GSP start', 'GSP end', 'GSP seq'])
    used_d = pd.DataFrame(columns=['chr', 'snGSP start', 'snGSP end', 'snGSP seq'])
    std = pd.DataFrame(columns=['name_id', 'chr', 'Product GC percent', 'product size',
                                'GSP ID', 'GSP seq', 'GSP start', 'GSP end', 'GSP strand', 'GSP score',
                                'snGSP ID', 'snGSP seq', 'snGSP start', 'snGSP end', 'snGSP strand', 'snGSP score',
                                'tagGSP ID', 'tagGSP seq', 'tagGSP start', 'tagGSP end', 'tagGSP strand', 'tagGSP score',
                                'GSP snp', 'snGSP snp', 'tagGSP snp'])

    global pd_extensible, mg_conc, oligo_conc, monovalent_cation_conc, dntp_conc, pd_temperature
    global oligo_conc, oligo_dG
    pd_extensible = 1
    mg_conc = 3  # our mg2+ is 3nM, 1.5 default
    oligo_conc = 50
    monovalent_cation_conc = 50
    dntp_conc = 0.8  # our dntp is 0.2*4=0.8, 0.25 default
    pd_temperature = 37
    oligo_dG = {}

