# -*- coding:utf-8 -*-
"""
@File   : RecommendPrimer.py
@Time   : 2025/2/5 14:00
@Author : Tian Chen
@Email  : chentian@genomics.cn
"""
import sys
import re
import pandas as pd

def load_two_chain_result(path, cover_amp_file, primer_result):
    try:
        print(f"{path}{cover_amp_file}")
        target = pd.read_csv(f"{path}/{cover_amp_file}", sep='\t', header=None)
        target.head()
        target.columns = ['chrom', 'start', 'end', 'name_id']
    except IOError as e:
        raise e
    try:
        ct_chain = pd.read_excel(f"{path}/{primer_result}", sheet_name="chain_CT_primer")
        ga_chain = pd.read_excel(f"{path}/{primer_result}", sheet_name="chain_GA_primer")
    except IOError as e:
        raise e
    ct_chain['U_merge_cnt'] = ct_chain['GSP seq'].apply(lambda x: x.count('Y'))
    ct_chain['D_merge_cnt'] = ct_chain['snGSP seq'].apply(lambda x: x.count('Y'))
    ct_chain['T_merge_cnt'] = ct_chain['tagGSP seq'].apply(lambda x: x.count('R'))
    ct_chain['score'] = 1/((ct_chain['GSP score']+1)*(2**ct_chain['U_merge_cnt'])) + 1/((ct_chain['snGSP score']+1)*(2**ct_chain['D_merge_cnt'])) + 1/((ct_chain['tagGSP score']+1)*(2**ct_chain['T_merge_cnt']))

    ga_chain['U_merge_cnt'] = ga_chain['GSP seq'].apply(lambda x: x.count('Y'))
    ga_chain['D_merge_cnt'] = ga_chain['snGSP seq'].apply(lambda x: x.count('Y'))
    ga_chain['T_merge_cnt'] = ga_chain['tagGSP seq'].apply(lambda x: x.count('R'))
    ga_chain['score'] = 1 / ((ga_chain['GSP score'] + 1) * (2 ** ga_chain['U_merge_cnt'])) + 1 / ((ga_chain['snGSP score'] + 1) * (2 ** ga_chain['D_merge_cnt'])) + 1 / ((ga_chain['tagGSP score'] + 1) * (2 ** ga_chain['T_merge_cnt']))

    target['name'] = target['name_id'].apply(lambda x: "_".join(x.split('_')[:-2]))
    ct_chain['target_name'] = ct_chain['name_id'].apply(lambda x: "_".join(x.split('_')[:-2]))
    ga_chain['target_name'] = ga_chain['name_id'].apply(lambda x: "_".join(x.split('_')[:-2]))

    ct_cover, ga_cover = ct_chain.target_name.tolist(), ga_chain.target_name.tolist()
    target_set = set(target['name'].tolist())
    #res = pd.DataFrame(columns=ct_chain.columns)
    res = []
    print(len(target_set), target_set)
    for name in target_set:
        ct_flag, ga_flag = False, False
        if name in ct_cover:
            ct_res = ct_chain[ct_chain['target_name'] == name]
            ct_flag = True
        if name in ga_cover:
            ga_res = ga_chain[ga_chain['target_name'] == name]
            ga_flag = True
        if ct_flag and ga_flag:
            print(name)
            print(ct_res)
            print(ga_res)
            if ct_res['score'].tolist()[0] >= ga_res['score'].tolist()[0]:
                #res.concat(ct_res)
                #res = pd.concat([res, ct_res])
                res.extend(ct_res.values.tolist())
            else:
                res.extend(ga_res.values.tolist())
        #elif len(ct_res) > 0:
        elif ct_flag:
            res.extend(ct_res.values.tolist())
        #elif len(ga_res) > 0:
        elif ga_flag:
            res.extend(ga_res.values.tolist())
    #print(type(res))
    #print(res)
    result = pd.DataFrame(data=res, columns=ct_chain.columns)
    result.to_excel(f'{path}/selected_primers.xlsx', index=False)

def trans_GSP_T2U(path, primer_file, output):
    try:
        primer_df = pd.read_excel(f"{path}/{primer_file}")
    except IOError as e:
        raise e
    #primer_df["new_GSP_seq"] = primer_df['GSP seq'].apply(lambda x: T2U_at3_end(x))
    primer_df["new_GSP_seq_mid_U"] = primer_df['new_GSP_seq'].apply(lambda x: process_primer(x, verbose=False))
    primer_df.to_excel(path+f'{output}')


def T2U_at3_end(seq):
    # dege_base = {'Y', 'G'}
    for i in range(len(seq)-1, -1, -1):
        if seq[i] == 'T':
            return seq[:i] + 'U' + seq[i+1:]
    return seq

def find_closest_t_to_center(primer_seq):
    """
    找到引物序列中最接近中间位置的T碱基的位置
    
    参数:
        primer_seq: 引物序列字符串
        
    返回:
        最接近中间位置的T的索引，如果没有T则返回-1
    """
    if not primer_seq:
        return -1
    
    seq_len = len(primer_seq)
    center = seq_len / 2
    
    # 查找所有T的位置
    t_positions = [i for i, base in enumerate(primer_seq) if base.upper() == 'T']
    
    if not t_positions:
        return -1
    
    # 找到最接近中心的T
    closest_pos = min(t_positions, key=lambda pos: abs(pos - center))
    return closest_pos

def replace_t_with_u(primer_seq, position):
    """
    将指定位置的T替换为U
    
    参数:
        primer_seq: 原始引物序列
        position: 需要替换的位置索引
        
    返回:
        修改后的序列字符串
    """
    if position < 0 or position >= len(primer_seq):
        return primer_seq
    
    # 创建新序列（字符串不可变，需要重建）
    seq_list = list(primer_seq)
    if seq_list[position].upper() == 'T':
        seq_list[position] = 'U'
    else:
        print(f"警告：位置{position}不是T碱基，实际是{seq_list[position]}，跳过替换")
        return primer_seq
    
    return ''.join(seq_list)

def process_primer(primer_seq, verbose=True):
    """
    处理引物序列的主要函数
    
    参数:
        primer_seq: 输入引物序列
        verbose: 是否打印详细信息
        
    返回:
        修改后的引物序列
    """
    # 移除可能的空白字符
    primer_seq = primer_seq.strip().upper()
    
    if verbose:
        print(f"原始序列: {primer_seq}")
        print(f"序列长度: {len(primer_seq)}")
    
    # 验证序列是否只包含ATCGU
    valid_bases = set('ATCGU')
    if not set(primer_seq.upper()).issubset(valid_bases):
        invalid = set(primer_seq.upper()) - valid_bases
        print(f"警告：序列包含非标准碱基: {invalid}")
    
    # 找到最接近中间位置的T
    closest_t_pos = find_closest_t_to_center(primer_seq)
    
    if closest_t_pos == -1:
        print("错误：序列中没有找到T碱基，无法进行替换")
        return None
    
    if verbose:
        center = len(primer_seq) / 2
        print(f"序列中心位置: {center:.1f}")
        print(f"找到最接近中心的T在位置: {closest_t_pos}")
        print(f"距离中心距离: {abs(closest_t_pos - center):.1f}")
    
    # 替换T为U
    modified_seq = replace_t_with_u(primer_seq, closest_t_pos)
    
    if verbose:
        print(f"修改后序列: {modified_seq}")
        print(f"替换位置: {closest_t_pos+1} (1-based索引)")
        
        # 显示比对
        spacer = ' ' * closest_t_pos
        print(f"比对图示:    {primer_seq}")
        print(f"             {spacer}↓")
        print(f"             {modified_seq}")
    
    return modified_seq

if __name__ == '__main__':
    path = "/home/pcr/work/ATOPlexMethylation/project/2026/BGI_multi_cancer/redesign/BGI_DMR_multi_cancer_CpG_primer/"
    cover_amplicon = 'amplicon_region.bed'
    primer_file = 'BGI_DMR_multi_cancer_CpG_primer.xlsx'
    #load_two_chain_result(path, cover_amplicon, primer_file)
    trans_GSP_T2U(path, 'multi-canser_DMR_CpG_redesign_primers_20260522.xlsx', 'new_multi-canser_DMR_CpG_redesign_primers_20260522.xlsx')
