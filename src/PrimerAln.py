"""
@File   : PrimerAln.py
@Time   : 2023/11/16 10:14
@Author : Tian Chen
@Email  : chentian@genomics.cn
"""
import os
import pandas as pd
import itertools as it

import globalVar as gl


merge_base = {
    'Y':'CT', 'R':'GA',
    'S':'GC', 'W':'AT',
    'K':'GT', 'M':'AC',
    'H':'ATC', 'B':'GTC',
    'V':'GAC', 'D':'GAT',
#    'N':'ATGC', 'I':'ATCG',
}

def generate_fq_file(path, primer_res, is_contain_adapter=True):
    # try:
    #     data1 = pd.read_excel(path+'/'+primer_file, sheet_name='chain_CT_primer')
    #     data2 = pd.read_excel(path+'/'+primer_file, sheet_name='chain_GA_primer')
    #     data = pd.concat([data1,data2])
    # except IOError as e:
    #     print(f"No found {path}/{primer_file}")
    #     raise e
    if is_contain_adapter:
        primer_res['snGSP seq'] = primer_res['snGSP seq'].apply(lambda x: x[21:])
        primer_res['tagGSP seq'] = primer_res['tagGSP seq'].apply(lambda x: x[17:])
    # primer_u, primer_d, primer_t = list(), list(), list()
    with open(path+'/GSP.fq', 'w') as fu, open(path+'/snGSP.fq', 'w') as fd, open(path+'/tagGSP.fq', 'w') as ft:
        for idx, row in primer_res.iterrows():
            primer_u = primer_expand(row['GSP seq'])
            primer_d = primer_expand(row['snGSP seq'])
            primer_t = primer_expand(row['tagGSP seq'])
            udt = 0
            for t_seq in primer_t:
                t_qua = 'F' * len(t_seq)
                for u_seq in primer_u:
                    u_qua = 'F' * len(u_seq)
                    for d_seq in primer_d:
                        d_qua = 'F' * len(d_seq)
                        udt += 1
                        fu.write(f"@{row.name_id}.{udt}\n{u_seq}\n+\n{u_qua}\n")
                        fd.write(f"@{row.name_id}.{udt}\n{d_seq}\n+\n{d_qua}\n")
                        ft.write(f"@{row.name_id}.{udt}\n{t_seq}\n+\n{t_qua}\n")


def primer_expand(primer_seq, primer_prefix=['']):
    ''' Returns a list containing degenerate nucleic acid sequences based on an amino acid sequence.
    Reference: http://lists.open-bio.org/pipermail/biopython/2006-September/003190.html '''
    if len(primer_seq) == 0:
        return  primer_prefix
    else:
        first = primer_seq[0]
        last_seq = primer_seq[1:]
        new_prefix = []
        for prefix in primer_prefix:
            if first.upper() in 'ATCG':
                new_prefix.append(prefix+first)
            elif first.upper() in merge_base.keys():
                for base in merge_base[first.upper()]:
                    new_prefix.append(prefix+base)
            else:
                print(f"can't recognized base '{first}' in {primer_seq}, please check it.")
        return primer_expand(last_seq, new_prefix)


def create_sam_file(ref_path, path, primer_file):
    generate_fq_file(path, primer_file, True)
    bwa_aln_sai(ref_path, path, 'GSP.fq', 'aln_GSP.sai')
    bwa_aln_sai(ref_path, path, 'snGSP.fq', 'aln_snGSP.sai')
    bwa_aln_sai(ref_path, path, 'tagGSP.fq', 'aln_tagGSP.sai')
    bwa_sampe_aln(ref_path, path, 'aln_GSP.sai', 'aln_tagGSP.sai', 'GSP.fq', 'tagGSP.fq', 'GSP_tagGSP_aln.sam')
    bwa_sampe_aln(ref_path, path, 'aln_snGSP.sai', 'aln_tagGSP.sai', 'snGSP.fq', 'tagGSP.fq', 'snGSP_tagGSP_aln.sam')
    sam2bam(path, 'GSP_tagGSP_aln.sam', 'GSP_tagGSP_aln.bam')
    sam2bam(path, 'snGSP_tagGSP_aln.sam', 'snGSP_tagGSP_aln.bam')
    bam_index(path, 'GSP_tagGSP_aln.bam')
    bam_index(path, 'snGSP_tagGSP_aln.bam')
    # os.system(f'rm {path}/*.sai {path}/*.sam {path}/*SP.fq')
    # os.system(f'rm {path}/*.sai {path}/*SP.fq')

def bwa_aln_sai(ref_path, path, input_fq, output_sai):
    """bwa 构建sai 索引"""
    os.system(f'bwa aln -t 20 {ref_path} {path}/{input_fq} 1>{path}/{output_sai} 2>>{path}/logs/bwa.log')


def bwa_sampe_aln(ref_path, path, sai_file1, sai_file2, fq1, fq2, output_sam):
    """bwa sampe 真实比对"""
    # os.system(f'bwa sampe -a 10000 -A {gl.Ref_path}  {gl.path}/{sai_file1} {gl.path}/{sai_file2} {gl.path}/
    # {fq1} {gl.path}/{fq2} >{gl.path}/{output_sam}')
    os.system(f'bwa sampe {ref_path}  {path}/{sai_file1} {path}/{sai_file2} {path}/{fq1} {path}/{fq2} '
              f'1>{path}/{output_sam} 2>>{path}/logs/bwa.log')


def sam2bam(path, input_sam_file, output_bam_file):
    """sam 格式 转bam 格式"""
    os.system(f'samtools view -bSu {path}/{input_sam_file}|samtools sort 1>{path}/{output_bam_file} '
              f'2>>{path}/logs/samtools.log')


def bam_index(path, bam_file):
    """bam 文件构建索引"""
    os.system(f'samtools index {path}/{bam_file} 2>>{path}/logs/samtools.log')


def load_head(sam_file):
    """加载头文件"""
    with open(sam_file, 'r') as f:
        heads = f.read()
    return heads


if __name__ == '__main__':
    # ref_path = '/home/pcr/D1/QuickPrimer/RunPrimer/ATOPlex_Designer/2023/DENV1-4/data/meta_DENV1/meta_NC_001477.1.fa'
    path = 'G:\\0.Work document\\2023\\ATOPlex-个性化设计\\冠状病毒\\'
    file = 'primers.xlsx'
    # create_sam_file(ref_path, path, file)
    # res = replace_base_seq(seq='AAAAGACGGSTCGSCCGT', m_base='S')
    # print([x for x in res])
    generate_fq_file(path, file, False)
