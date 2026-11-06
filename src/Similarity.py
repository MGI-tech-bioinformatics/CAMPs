#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @author: Tian Chen
# @times: 2020/9/25  15:23
import re
import pandas as pd
import numpy as np


def edit_dis(line):
    mismatch = 0
    for i in re.compile('([0-9]+)[DIH]').findall(line[5]):
        mismatch += int(i)
    return mismatch+int(line[11].split(':')[-1])


def HomoSim(file, merge_dic):
    """
    通过bwa比对文件，找到其同源位点，并计算同源位点与对应目标位点的相似性(simliarity)和同源数目(dupnum)
    :param file: sam文件
    :param merge_dic: the target position which have been merge
    :return: 包含name、simliarity、dupnum的Dataframe
    """
    index = 0
    with open(file, 'r') as r:
        while True:
            line = r.readline()
            if not line.startswith('@'):
                break
        line = line.split()
        rename = line[0]
        length = len(line[9])
        dup_num, dup_region = 0, []
        # is_first = True
        sim = []
        line = r.readline().split()
        Tag = pd.DataFrame(columns=['name', 'sim', 'sim_cnt', 'sim_region'])
        while line:
            if line[0] == rename:
                # dup_num += 1
                similarity = 1 - edit_dis(line) / length
                if similarity >= 0.8:
                    sim.append(similarity)
                    dup_region.append(f'{line[2]}:{line[3]}:{line[5]}:{line[11]}')
                    dup_num += 1
            else:
                Tag.loc[index] = [rename, round(np.mean(sim), 2), dup_num, dup_region]
                # Tag.loc[index] ={'name':rname,'sim':simliarity,'simsum':dupnum}
                rename = line[0]
                length = len(line[9])
                dup_num, dup_region = 0, []
                # is_first = True
                similarity = 0
                sim = []
                index += 1
            line = r.readline().split()
    Tag.loc[index] = [rename, round(np.mean(sim), 2), dup_num, dup_region]
    Tag.fillna(0, inplace=True)
    Tag['chr'] = Tag['name'].str.split(':').str[1]
    Tag['start'] = Tag['name'].str.split(':').str[2]
    Tag['end'] = Tag['name'].str.split(':').str[3]
    Tag['length'] = Tag['name'].str.split(':').str[4]
    Tag['design_flag'] = Tag['name'].str.split(':').str[5]
    Tag['name'] = Tag['name'].str.split(':').str[0]
    Tag['TarLen'] = Tag['length']
    Tag = Tag[['name', 'chr', 'start', 'end', 'length', 'design_flag', 'sim', 'sim_cnt', 'TarLen','sim_region']]
    Tag['start'] = Tag['start'].astype(int)
    Tag['end'] = Tag['end'].astype(int)
    Tag['length'] = Tag['length'].astype(int)
    Tag['TarLen'] = Tag['TarLen'].astype(int)
    Tag['sim'] = Tag['sim'].astype('float64')
    Tag['sim_cnt'] = Tag['sim_cnt'].astype(int)
    Tag['sim_region'] = Tag['sim_region'].apply(lambda x: np.array(x))
    row, col = Tag.shape
    for idx in range(row):
        try:
            set_i = [x for x in range(Tag.iloc[idx, 2], Tag.iloc[idx, 3]+1)]
            # print(Tag.iloc[idx, 1])
            intersection = set(set_i).intersection(set(merge_dic[Tag.iloc[idx, 1]]))
            # print(set(merge_dic[Tag.iloc[idx, 1]]))
            Tag.iloc[idx, -2] = len(intersection)
        except:
            print("can't caculate the intersection")

    return Tag

#tag = HomoSim('Product_aln.sam')
#print(tag)
