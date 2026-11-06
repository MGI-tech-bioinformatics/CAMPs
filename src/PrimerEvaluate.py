"""
@File   : PrimerEvaluate.py
@Time   : 2023/11/15 9:51
@Author : Tian Chen
@Email  : chentian@genomics.cn
"""
import pandas as pd
import itertools as it

from multiprocessing import Process, Pool, Manager, Lock
from openpyxl import load_workbook

import Dimer
import globalVar as gl
import PrimerAln as PA


merge_base = {
    'Y':'CT', 'R':'GA',
    'S':'GC', 'W':'AT',
    'K':'GT', 'M':'AC',
}

def dimers_detect(path, file, is_contain_adapter=True):
    """针对含简并碱基的引物进行扩充，并评估二聚体"""
    try:
        data = pd.read_excel(path+'/'+file)
    except IOError as e:
        print(f"No found {path}/{file}")
        raise e
    if is_contain_adapter:
        data['snGSP seq'] = data['snGSP seq'].apply(lambda x: x[21:].upper())
        data['tagGSP seq'] = data['tagGSP seq'].apply(lambda x: x[17:].upper())
    primer_u, primer_d, primer_t = list(), list(), list()
    for idx, row in data.iterrows():
        primer_u.extend([f"{row['GSP ID']}.{idx + 1}\n{seq}" for idx, seq in enumerate(PA.primer_expand(row['GSP seq']))])
        primer_d.extend([f"{row['snGSP ID']}.{idx + 1}\n{seq}" for idx, seq in enumerate(PA.primer_expand(row['snGSP seq']))])
        primer_t.extend([f"{row['tagGSP ID']}.{idx + 1}\n{seq}" for idx, seq in enumerate(PA.primer_expand(row['tagGSP seq']))])
        # primer_u.extend(expand_seqs(row.name_id, row['GSP seq'], flag='-U'))
        # primer_d.extend(expand_seqs(row.name_id, row['snGSP seq'], flag='-D'))
        # primer_t.extend(expand_seqs(row.name_id, row['tagGSP seq'], flag='-T'))
    # sub_seq_dimer(path, 'dimer_test.xlsx', primer_u, primer_u, 'GSP_GSP')
    lock = Lock()
    p1 = Process(target=sub_seq_dimer, args=(path, file, primer_u, primer_u, lock, 'GSP_GSP'))
    p2 = Process(target=sub_seq_dimer, args=(path, file, primer_u, primer_t, lock, 'GSP_tagGSP'))
    p3 = Process(target=sub_seq_dimer, args=(path, file, primer_d, primer_t, lock, 'snGSP_tagGSP'))
    p4 = Process(target=sub_seq_dimer, args=(path, file, primer_d, primer_d, lock, 'snGSP_snGSP'))
    p5 = Process(target=sub_seq_dimer, args=(path, file, primer_t, primer_t, lock, 'tagGSP_tagGSP'))
    p1.start()
    p2.start()
    p3.start()
    p4.start()
    p5.start()
    p1.join()
    p2.join()
    p3.join()
    p4.join()
    p5.join()


def sub_seq_dimer(path, file, sub1, sub2, lock, flag='GSP_GSP'):
    manager = Manager()
    dimer_dict = manager.dict()
    p = Pool(10)
    for row in sub1:
        name1, seq1 = row.split('\n')
        p.apply_async(two_seq_dimer, args=(name1, seq1, sub2, dimer_dict))
    p.close()
    p.join()

    dimer = pd.DataFrame(dimer_dict.values(), columns=['name1:name2', 'DG', 'dimer'])
    dimer.sort_values(by='DG', ascending=True, inplace=True)
    # dimer.to_excel(f"{path}/{file}", index=False)

    lock.acquire()
    book = load_workbook(f"{path}/{file}")
    writer = pd.ExcelWriter(f"{path}/{file}", engine='openpyxl')
    writer.book = book
    dimer[dimer['DG'] < -4].to_excel(writer, sheet_name=f'{flag}_dimer_sort', index=False)
    writer.save()
    lock.release()



def two_seq_dimer(name1, seq1, seq2_lst, dimer_dic):
    seq1 = seq1.upper()
    for row in seq2_lst:
        name2, seq2 = row.split('\n')
        seq2 = seq2.upper()
        # print(seq1, seq2)
        dg, dimer = detect_dimer(seq1, seq2)
        # print(name1, name2, seq1, seq2, dg, dimer)
        if dg and dg < -4:
            # print(name1, name2, seq1, seq2, dg, dimer)
        # if dg < 0:
            dimer_dic[f'{name1}:{name2}'] = [f'{name1}:{name2}', dg, dimer]


def detect_dimer(seq1, seq2):
    gl.oligo_dG = {}
    Dimer.recalculate_dG(37)
    bind_string = []
    DG, index, DG_lst = Dimer.primer_dimer(seq1, seq2, 0, bind_string)
    if DG:
        # print(DG, index, DG_lst[index][1])
        dimer = Dimer.draw_dimer(seq1, seq2, DG_lst[index][1], index, bind_string)
        if DG < 0:
            return DG, dimer

    return DG, None


def expand_seqs(name_id, primer_seq, flag='-U'):
    res = []
    if 'Y' in primer_seq:
        merge_base = 'Y'
    elif 'R'in primer_seq:
        merge_base = 'R'
    elif 'W' in primer_seq:
        merge_base = 'W'
    elif 'S' in primer_seq:
        merge_base = 'S'
    elif 'M' in primer_seq:
        merge_base = 'M'
    elif 'K' in primer_seq:
        merge_base = 'K'
    else:
        res.append(f'{name_id}{flag}\n{primer_seq}')
        return res
    for idx, prim in enumerate(replace_base_seq(primer_seq, merge_base)):
        res.append(f'{name_id}{flag}.{idx+1}\n{prim}')
    return  res


def replace_base_seq(seq, m_base='Y'):
    bases = merge_base[m_base]
    seqs_lst = seq.split(m_base)
    merge_lst = [list(x) for x in it.product(bases, repeat=len(seqs_lst)-1)]
    for b in merge_lst:
        new_seq = list("".join(x) for x in list(zip(seqs_lst[:-1], b))) + seqs_lst[-1:]
        yield "".join(new_seq)


if __name__ == '__main__':
    path = '/home/pcr/work/ATOPlexMethylation/ATOPlexMethDesigner/merge_dimer_test/'
    file = 'giant_panda_snp_STR_primer_20260624.xlsx'
    dimers_detect(path, file, True)
