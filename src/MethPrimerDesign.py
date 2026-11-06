"""
@File   : MethPrimerDesign.py
@Time   : 2023/11/1 17:21
@Author : Tian Chen
@Email  : chentian@genomics.cn
"""
import os
import re
import sys
import math
import logging

import pandas as pd
import numpy as np
from multiprocessing import Pool, Manager, Lock

import Seq
import TmDeltaG
import CombinePrimer
import globalVar as gl
import PrepareData as PreD


def load_hit_index(path):
    """ 加载候选序列比对特异性数值 """
    index_file = path + '/Hits_result.ln.gz.idx'
    cnt = 0
    # hit_idx = {}
    try:
        with open(index_file, 'r') as idx:
            for line in idx:
                cnt += 1
                if cnt < 4: continue
                k = line.replace('\n', '').split('\t')
                gl.hit_idx[k[0]] = k[2:]  # 第２列数值没有意义
        if len(gl.hit_idx) == 0:
            sys.stderr.write(f"primersDB run failed.\n")
            logging.info(f"primersDB run failed.")
        else:
            logging.info("load hit_idx have done")
    except:
        logging.info("load hit_idx failed")
        sys.stderr.write(f"Can't find the primersDB result, please check it.\n")
        sys.exit()
    # return hit_idx

def load_unique_pos(path, file_name):
    """加载同源分析后的目标区域，并去除_CT/_GA tag, 对区域进行去重或合并
    name    chr     start   end     length  sim     sim_cnt
    """
    pos, pos_CT, pos_GA = [], [], []
    try:
        with open(path+'/'+file_name, 'r') as f:
            line = f.readline().strip().split()
            if '_CT' in line[0]:
                pos_CT.append(line)
            elif '_GA' in line[0]:
                pos_GA.append(line)
    except IOError as e:
        logging.error(f"No found {path}/{file_name}, please check it.")
        raise e


def split_seq(len_short, len_long, extend_range, tm_min, tm_max, script_dir, path, pos, snp_dic, ref_dic):
    """Load the target positions candidate primers"""
    # chr_ref = load_chr_ref(db_dir+pos[1]+'.fa')
    chr_ref = ref_dic[pos[1]]
    # manager = Manager()
    # res_dic = manager.dict()
    # P1 = Process(target=get_primers, args=(opt, pos, '+', chr_ref, snp_dic, res_dic,))
    # P2 = Process(target=get_primers, args=(opt, pos, '-', chr_ref, snp_dic, res_dic,))
    f_primers = get_primers(len_short, len_long, extend_range, tm_min, tm_max,
                            script_dir, path, pos, '+', chr_ref, snp_dic)
    r_primers = get_primers(len_short, len_long, extend_range, tm_min, tm_max,
                            script_dir, path, pos, '-', chr_ref, snp_dic)
    # P1.start()
    # P2.start()
    # P1.join()
    # P2.join()
    return f_primers, r_primers, chr_ref


def get_primers(len_short, len_long, extend_range, tm_min, tm_max, script_dir, path, pos, strand, ref_seq, snp_dic):
    """获取候选序列"""
    start = pos[2] - extend_range if strand == '+' else pos[3] + 1
    end = pos[2] - len_short if strand == '+' else pos[3] + extend_range - len_short
    # pool = Pool(processes=10)
    # manager = Manager()
    # res_dic = manager.dict()
    primers = []
    for i in range(start, end):
        # pool.apply_async(filter_primers, args=(opt.len_short, opt.len_long, opt.extend_range, opt.tm_min, opt.tm_max,
        #                                        opt.script_dir, opt.path, pos, end, i, ref_seq, strand, snp_dic, res_dic))
        primers.extend(filter_primers(len_short, len_long, extend_range, tm_min, tm_max,
                                      script_dir, path, pos, end, i, ref_seq, strand, snp_dic))

    # pool.close()
    # pool.join()
    primers_df = pd.DataFrame(data=primers, columns=['chr', 'start', 'end', 'len', 'seq', 'strand',
                                                              'hits', 'gc', 'tm', 'cg_cnt', 'poly', 'snp', 'score'])
    primers_df['score'] = primers_df['score'].astype(int)
    primers_df.sort_values(by=['score', 'cg_cnt'], ascending=[True, True], inplace=True)

    # res_dic[strand] = primers_df.values
    return primers_df.values


def filter_primers(len_short, len_long, extend_range, tm_min, tm_max, script_dir, path,
                   pos, end, i, ref_seq, strand, snp_dic):
    """根据设定参数过滤不符合要求的短序列"""
    res = []
    for j in range(len_short, len_long + 1):
        sub_seq = ref_seq[i - 1:i + j - 1]
        if 'N' in sub_seq:
            continue
        if (i + j >= end + len_short or i - 1 < 0 or len(sub_seq.replace(' ', '')) < len_short) or \
                adapter_deltaG(Seq.complement(sub_seq)) < -10 or detect_hairpin(sub_seq, strand):
            # print(f'skip ... {sub_seq}')
            continue
        primer = Seq.Sequence(chr=pos[1], start=i, seq=sub_seq, Prange=extend_range,
                              Flag=strand, snp_dic=snp_dic)
        # test print
        # 剔除不满足poly、TM值、snp要求的候选引物

        if (primer.poly > 0 or primer.tm < tm_min - 1.8 * primer.cg  or primer.tm > tm_max - 1.8 * primer.cg
                or primer.Maf > 0.005 or primer.cg > 2 or primer.new_seq == ''):  # cg < 1
            # if strand == '-':
            # print("primer skip ...", pos[1], i, primer.end, primer.len, primer.new_seq)
            continue
        hits = Seq.hits(script_dir, path, chr_name=pos[1], start=i, lens=primer.len, strand=strand)
        #hits += 10000 * primer.cg  # 2022-5-25, increasing the effection from cg count
        if hits > 15000:
            # if strand == '-':
            #     print('hit skip ...')
            continue
        # if strand == '-':
        #     print(pos[1], i, primer.end, primer.len, primer.new_seq, primer.strand,
        #        hits, primer.gc, primer.tm, hits, primer.poly, primer.snp)
        score = 1 / ((hits + 1) * (2 ** primer.cg))
        #stm = primer_score(primer.tm, 10, 56, 58, 52, 62, 50, 64)
        #sgc = primer_score(primer.gc, 10, 40, 60, 30, 70, 25, 75)
        #score_hits = round(math.log10(1/((hits+1)*(2 ** primer.cg))),2)
        #score = 2 * stm + sgc - 2 * primer.cg + 2 *  score_hits
        res.append([pos[1], i, primer.end, primer.len, primer.new_seq, primer.strand,
               hits, primer.gc, primer.tm, primer.cg, primer.poly, primer.snp, score])

        # res_dic[f"{pos[1]}_{i}_{j}"] = res
    return res


def primer_score(value, score, min_b, max_b, min_v, max_v, min_l, max_l ):
    dis_l, dis_r = min_b - min_v, max_v - max_b
    if dis_l < 0 or dis_r < 0:
        print(f"Wrong Score set: {min_b}, {max_b}, {min_v}, {max_v}")
        exit(1)
    K, a, b, e = 1.0, 1.0, 1.0, 2.718
    if value < min_b:
        x = (value - min_l)/((min_b - min_l)/10) - 5
        y1 = K / round((1 + b * e ** (-1 * a * x)),2)
        x = (min_v - min_l)/((min_b - min_l)/10) - 5
        y10 = K / round((1 + b * e ** (-1 * a * x)),2)
        s = score * (y1 - y10) / (1 - y10)
    elif value <= max_b:
        s = score
    else:
        if max_l - max_b == 0:
            print(f"Wrong: The value exceeds the upper limit {max_l}, {max_b}")
            exit(1)
        x = (max_b - value) / ((max_l - max_b)/10) + 5
        y2 = K / round((1 + b * e ** (-1 * a * x)),2)
        x = (max_b - max_v) / ((max_l - max_b)/10) + 5
        y20 = K / round((1 + b * e ** (-1 * a * x)),2)
        s = score * (y2 - y20) / (1 - y20)
    #print(s)
    return round(s, 2)

def adapter_deltaG(seq):
    """检测序列引物与接头的结合能"""
    adapter = "TTGTCTTCCTAAGACCGCTTGGCCTCCGACTT"
    deltaG = TmDeltaG.calDeltaG(adapter, seq, mono_conc=50, diva_conc=3, dntp_conc=0.8, deltaH=None, deltaS=None)
    return deltaG


def detect_hairpin(seq, strand):
    """检查引物是否存在反向互补结构，形成发夹结构"""
    x = lambda s: Seq.self_hairpin(seq) if s == '+' else Seq.self_hairpin(Seq.rev_com(seq))
    if x(strand):
        # print(seq, strand, Seq.rev_com(seq))
        return True
    return False


def cal_amp_gc_avg_median(res_df):
    # 计算扩增子的GC 和长度的均值及中位数
    gc_avg = round(res_df['Product GC percent'].mean(), 2)
    gc_median = round(res_df['Product GC percent'].median(), 2)
    amp_avg = round(res_df['product size'].mean(), 2)
    amp_median = round(res_df['product size'].median(), 2)

    return gc_avg, gc_median, amp_avg, amp_median



def design_summary_log(result_df, res_ct, res_ga, failed_cnt, opt):
    # 引物设计结果总结日志输出
    primer_cnt, ct_primer_cnt, ga_primer_cnt = len(result_df), len(res_ct), len(res_ga)
    fail_cnt_ct, fail_cnt_ga, fail_cnt, fail_len_ct, fail_len_ga, fail_len, cvg_ratio_ct, cvg_ratio_ga, cvg_ratio = coverage_ratio(opt.path, opt.path+'/Input.bed', opt.output + '.xlsx')
    # failed_ratio = round((float(int(fail_len) / gl.report_dic['mergeTargetLen']) * 100), 2)
    # cvg_ratio = round((1 - (failed_len + gl.report_dic['abandonHomoLen']) / gl.report_dic['mergeTargetLen']) * 100, 2)
    amp_gc_avg, amp_gc_median, amp_len_avg, amp_len_median = cal_amp_gc_avg_median(result_df)
    ct_amp_gc_avg, ct_amp_gc_median, ct_amp_len_avg, ct_amp_len_median = cal_amp_gc_avg_median(res_ct)
    ga_amp_gc_avg, ga_amp_gc_median, ga_amp_len_avg, ga_amp_len_median = cal_amp_gc_avg_median(res_ga)
    logging.info(f"\n{'=='*30}\n"
                 f"Total primer groups count: {primer_cnt}\n"
                 f"Total failed target regions count: {failed_cnt}\n"
                 f"Total uncover target regions length: {fail_len}\n"
                 # f"Failed targets region ratio: {failed_ratio}%\n"
                 f"Total target regions coverage ratio: {cvg_ratio}%\n"
                 f"Total amplicons gc content average: {amp_gc_avg}%\n"
                 f"Total amplicons gc content median:  {amp_gc_median}%\n"
                 f"Total amplicons length average: {amp_len_avg} bp\n"
                 f"Total amplicons length median: {amp_len_median} bp\n"
                 f"{'=='*30}\n"
                 f"C>T strain primer groups count: {ct_primer_cnt}\n"
                 f"C>T strain uncover target regions length: {fail_len_ct}\n"
                 f"C>T strain target regions coverage ratio: {cvg_ratio_ct}\n"
                 f"C>T strain amplicons gc content average: {ct_amp_gc_avg}%\n"
                 f"C>T strain amplicons gc content median: {ct_amp_gc_median}%\n"
                 f"C>T strain amplicons length average: {ct_amp_len_avg}bp\n"
                 f"C>T strain amplicons length median: {ct_amp_len_median}bp\n"
                 f"{'=='*30}\n"
                 f"G>A strain primer group counts: {ga_primer_cnt}\n"
                 f"G>A strain uncover target region length: {fail_len_ga}\n"
                 f"G>A strain primers cover target ratio: {cvg_ratio_ga}\n"
                 f"G>A strain amplicons gc content average: {ga_amp_gc_avg}%\n"
                 f"G>A strain amplicons gc content median: {ga_amp_gc_median}%\n"
                 f"G>A strain amplicons length average: {ga_amp_len_avg}bp\n"
                 f"G>A strain amplicons length median: {ga_amp_len_median}bp\n"
                 f"{'==' * 30}\n"
                 )
    # total primer on strain _CT and strain _GA
    gl.report_dic['primerCnt'] = primer_cnt
    gl.report_dic['FailedCnt'] = fail_cnt
    gl.report_dic['FailedLen'] = fail_len
    # gl.report_dic['FailedRatio'] = failed_ratio
    gl.report_dic['CoverageRatio'] = cvg_ratio
    gl.report_dic['gcMean'] = amp_gc_avg
    gl.report_dic['gcMedian'] = amp_gc_median
    gl.report_dic['ampliconMean'] = amp_len_avg
    gl.report_dic['ampliconMedian'] = amp_len_median
    # total primer on strain _CT
    gl.report_dic['CT_primerCnt'] = ct_primer_cnt
    gl.report_dic['CT_FailedCnt'] = fail_cnt_ct
    gl.report_dic['CT_FailedLen'] = fail_len_ct
    # gl.report_dic['CT_FailedRatio'] = failed_ratio
    gl.report_dic['CT_CoverageRatio'] = cvg_ratio_ct
    gl.report_dic['CT_gcMean'] = ct_amp_gc_avg
    gl.report_dic['CT_gcMedian'] = ct_amp_gc_median
    gl.report_dic['CT_ampliconMean'] = ct_amp_len_avg
    gl.report_dic['CT_ampliconMedian'] = ct_amp_len_median

    # total primer on strain _GA
    gl.report_dic['GA_primerCnt'] = ga_primer_cnt
    gl.report_dic['GA_FailedCnt'] = fail_cnt_ga
    gl.report_dic['GA_FailedLen'] = fail_len_ga
    # gl.report_dic['CT_FailedRatio'] = failed_ratio
    gl.report_dic['GA_CoverageRatio'] = cvg_ratio_ga
    gl.report_dic['GA_gcMean'] = ga_amp_gc_avg
    gl.report_dic['GA_gcMedian'] = ga_amp_gc_median
    gl.report_dic['GA_ampliconMean'] = ga_amp_len_avg
    gl.report_dic['GA_ampliconMedian'] = ga_amp_len_median


def load_target(target_file):
    """加载原始目标区域"""
    try:
        if '.xlsx' in target_file:
            data = pd.read_excel(target_file)
        else:
            data = pd.read_csv(target_file, sep='\s+', header=None, names=['chr', 'start', 'end', 'name'])
            data = data[['name', 'chr', 'start', 'end']]
            data['start'] = data['start'].astype(np.int)
            data['end'] = data['end'].astype(np.int)
    except:
        sys.stderr.write("Can't load the target region file when evaluate.")
        sys.exit()
    pos_dic = {}
    tar_len = 0
    for row in data.values:
        pos_dic.setdefault(str(row[0])+':'+str(row[1]), []).extend([x for x in range(row[2], row[3] + 1)])
        tar_len += row[3] - row[2] + 1
    return pos_dic, tar_len


def load_amplication(path, primers_file, tag_name):
    """加载扩增子区域"""
    df = pd.read_excel(path + '/' + primers_file, sheet_name=tag_name)

    df['amp_start'] = np.where(df['snGSP strand'] == '+', df['snGSP end'] + 1, df['tagGSP end'] + 1)  # correct uncover region, 20230410
    df['amp_end'] = np.where(df['snGSP strand'] == '+', df['tagGSP start'] - 1, df['snGSP start'] - 1)

    amp_dic = {}
    data = df[['name_id', 'chr', 'amp_start', 'amp_end']].values
    for row in data:
        if row[3] > row[2]:
            chr_name = row[1].split('_')[0].lower()
            if chr_name[-1] in ['x', 'y', 'm']:
                chr_name = chr_name[:-1] + chr_name[-1].upper()
            amp_dic.setdefault(chr_name, []).extend([x for x in range(row[2], row[3] + 1)])
        else:
            logging.info('cyclic chr', row)
    # df.to_excel(path + '/add_amp_snp_primers.xlsx', index=False)
    return amp_dic, df['product size'].sum()


def split_lst(f, chr_idx, lst):
    """拆分非连续的数组"""
    name = chr_idx.split(':')
    try:
        start = end = lst[0]
        for i in lst[1:]:
            if i - 1 > end:
                f.write(f"{name[0]}\t{name[-1]}\t{start}\t{end}\t{end - start + 1}\n")
                start = end = i
            else:
                end = i
        f.write(f"{name[0]}\t{name[-1]}\t{start}\t{end}\t{end - start + 1}\n")
    except :
        logging.info(f"{name[0]} uncovered")


def generate_bed(path, file, output):
    """generate the amplicon, cover region and uncover region bed file"""
    try:
        data = pd.read_csv(path+'/'+file, sep='\t', header=0)
    except IOError as e:
        raise e

    if len(data) > 0:
        data[['chr','start', 'end', 'name']].to_csv(path+'/'+output, sep='\t', index=False, header=False)


def get_uncover_region(path, failed_lst):
    print(failed_lst)
    failed = pd.DataFrame(data=failed_lst, columns=['name', 'chr', 'start', 'end', 'length', 'design_flag','sim.',
                                                       'HomoNum','failed reason'])
    try:
        un_cover_region = pd.read_csv(path+'/un_cover_target_pos.txt', sep='\t')
        # un_cover_region.drop_duplicates(subset=['start', 'end'], keep='last', inplace=True)
        un_cover_region.sort_values(by=['chr', 'start'], ascending=(True, True), inplace=True)
    except:
        logging.info("get un_cover_target_pos.txt failed")
        os._exit(1)

    map_failed = []
    for row in un_cover_region.values.tolist():
        # print(row)
        un_cover_set = set([ _ for _ in range(int(row[2]), int(row[3])+1)])
        failed_reason = []
        for row2 in failed[failed['chr'].apply(lambda x: x[:-3]) == row[1].upper()].values.tolist():
            failed_set = set([_ for _ in range(int(row2[2]), int(row2[3])+1)])
            diff = un_cover_set.intersection(failed_set)
            if 0 < len(diff) <= len(un_cover_set):
                un_cover_set = diff
                if not failed_reason:
                    failed_reason.append(row2[-1])
                elif row2[-1] not in failed_reason[0]:
                    failed_reason[0] = ",".join([failed_reason[0], row2[-1]])
        map_failed.append(row+failed_reason)
    show_failed = pd.DataFrame(data=map_failed, columns=['name', 'chr', 'start', 'end', 'length', 'failed reason'])
    return show_failed


def coverage_ratio(path, target_pos_file, primer_file):
    """计算扩增子在目标区的覆盖率"""
    target_pos, total_tar_len = load_target(target_pos_file)
    amp_pos_ct, ct_amp_len = load_amplication(path, primer_file, 'chain_CT_primer')
    amp_pos_ga, ga_amp_len = load_amplication(path, primer_file, 'chain_GA_primer')
    cover, cover_ct, cover_ga = 0, 0, 0
    logging.info(f'total target region base: {total_tar_len}')
    logging.info(f'total chain_CT amplicon region base: {ct_amp_len}')
    logging.info(f'total chain_CT amplicon region base: {ga_amp_len}')
    fct = open(path+'/gene_coverage_CT_ratio.txt', 'w')
    fga = open(path+'/gene_coverage_GA_ratio.txt', 'w')
    f   = open(path+'/gene_coverage_ratio.txt', 'w')
    f1 = open(path+'/un_cover_target_CT_pos.txt', 'w')
    f2 = open(path+'/un_cover_target_GA_pos.txt', 'w')
    f0 = open(path+'/un_cover_target_pos.txt', 'w')
    fct.write('name\tchr\tstart\tend\tcover_length\n')
    fga.write('name\tchr\tstart\tend\tcover_length\n')
    f.write('name\tchr\tstart\tend\tcover_length\n')
    f2.write("name\tchr\tstart\tend\tuncover_length\n")
    f0.write("name\tchr\tstart\tend\tuncover_length\n")
    f1.write("name\tchr\tstart\tend\tuncover_length\n")
    fail_cnt, fail_cnt_ct, fail_cnt_ga = 0, 0, 0
    for idx in target_pos.keys():
        chr_name = idx.split(':')[1]
        cover_set_ct = amp_pos_ct[chr_name] if chr_name in amp_pos_ct.keys() else set()
        cover_set_ga = amp_pos_ga[chr_name] if chr_name in amp_pos_ga.keys() else set()

        intersection_ct = list(set(target_pos[idx]).intersection(set(cover_set_ct)))
        intersection_ga = list(set(target_pos[idx]).intersection(set(cover_set_ga)))
        intersection = list(set(intersection_ct).union(set(intersection_ga)))
        intersection_ct.sort()
        intersection_ga.sort()
        intersection.sort()

        cover_ct += len(intersection_ct)
        cover_ga += len(intersection_ga)
        cover += len(intersection)
        #cover_length += max(len(intersection_ct), len(intersection_ga))
        split_lst(fct, idx, intersection_ct)
        split_lst(fga, idx, intersection_ga)
        split_lst(f, idx, intersection)
        # f.write(f"{idx}\t{len(intersection)}\t{len(target_pos[idx])}\t{round(len(intersection)/len(target_pos[idx]),4)}\n")
        failed_lst_ct = list(set(target_pos[idx]).difference(set(cover_set_ct)))
        failed_lst_ga = list(set(target_pos[idx]).difference(set(cover_set_ga)))
        failed_lst =  list(set(target_pos[idx]).difference(set(intersection)))
        failed_lst_ct.sort()
        failed_lst_ga.sort()
        failed_lst.sort()

        if len(failed_lst_ct) > 0:
            fail_cnt_ct += 1
            split_lst(f1, idx, failed_lst_ct)
        else:
            logging.info(f'{idx} all cover on chain_CT')

        if len(failed_lst_ga) > 0:
            fail_cnt_ga += 1
            split_lst(f2, idx, failed_lst_ga)
        else:
            logging.info(f'{idx} all cover on chain_GA')

        if len(failed_lst) > 0:
            fail_cnt += 1
            split_lst(f0, idx, failed_lst)

    f.close()
    f0.close()
    f1.close()
    f2.close()
    fct.close()
    fga.close()
    failed_len_ct = gl.report_dic['mergeTargetLen'] - cover_ct
    failed_len_ga = gl.report_dic['mergeTargetLen'] - cover_ga
    failed_len = gl.report_dic['mergeTargetLen'] - cover
    logging.info(f"failed_lst_CT: {fail_cnt_ct}, failed_len_CT: {failed_len_ct}")
    logging.info(f"failed_lst_GA: {fail_cnt_ga}, failed_len_GA: {failed_len_ga}")
    logging.info(f"failed_lst: {fail_cnt}, failed_len: {failed_len}")
    cover_ratio_ct = round(cover_ct / gl.report_dic['mergeTargetLen'] * 100, 2)
    cover_ratio_ga = round(cover_ga / gl.report_dic['mergeTargetLen'] * 100, 2)
    cover_ratio = round(cover / gl.report_dic['mergeTargetLen'] * 100, 2)
    #f.close()
    generate_bed(path, 'gene_coverage_CT_ratio.txt', 'cover_target_CT_region.bed')
    generate_bed(path, 'gene_coverage_GA_ratio.txt', 'cover_target_GA_region.bed')
    generate_bed(path, 'un_cover_target_CT_pos.txt', 'missing_target_CT_region.bed')
    generate_bed(path, 'un_cover_target_GA_pos.txt', 'missing_target_GA_region.bed')

    return fail_cnt_ct, fail_cnt_ga, fail_cnt, failed_len_ct, failed_len_ga, failed_len, cover_ratio_ct, cover_ratio_ga, cover_ratio


def load_target_region(path):
    """加载进行引物设计的目标区域文件"""
    try:
        with open(f'{path}/new_pos_tag.txt', 'r') as f:
            data = pd.read_csv(f'{path}/new_pos_tag.txt', sep='\t')
            data['name'] = data['name'].astype(str)
            return data.values
    except IOError as e:
        logging.info("No find target file", e)
        return None


def load_primer(pre_primer):
    """加载已经设计好的引物"""
    if re.search('csv', pre_primer):
        std = pd.read_csv(pre_primer)
    elif re.search('xlsx', pre_primer):
        std = pd.read_excel(pre_primer)
    else:
        logging.info("can't load the pre_primer file")
        sys.stderr.write("Load pre_primer file failed, please check the pre_primer file is exists.\n")
        sys.exit()
    try:
        del std['amp_start'], std['amp_end']
    except:
        logging.info("pre_primer file didn't have amp_start and amp_end column.")
    used_tr = std[std['tagGSP strand'] == '-'][['chr', 'tagGSP start', 'tagGSP end', 'tagGSP seq']]
    used_d = std[['chr', 'snGSP start', 'snGSP end', 'snGSP seq']]
    used_uf = std[std['GSP strand'] == '+'][['chr', 'GSP start', 'GSP end', 'GSP seq']]
    used_tf = std[std['tagGSP strand'] == '+'][['chr', 'tagGSP start', 'tagGSP end', 'tagGSP seq']]
    used_ur = std[std['GSP strand'] == '-'][['chr', 'GSP start', 'GSP end', 'GSP seq']]

    return used_d, used_tr, used_tf, used_uf, used_ur, std

def primer_design(opt, homo_0, ref_dic):
    """引物设计"""
    print(f"window size: {opt.window_size}")
    data = load_target_region(opt.path)
    # pos_idx = 0
    res_lst = Manager().list()
    used_d_lst = Manager().list()
    used_tr_lst = Manager().list()
    used_tf_lst = Manager().list()
    used_uf_lst = Manager().list()
    used_ur_lst = Manager().list()
    failed_lst = Manager().list()

    failed_lst.extend(list(homo_0))
    if opt.pre_primer:
        # load the pre_primer information.
        used_d, used_tr, used_tf, used_uf, used_ur, std = load_primer(opt.pre_primer.name)
        res_lst.extend(std.values)
        used_d_lst.extend(used_d.values)
        used_tr_lst.extend(used_tr.values)
        used_tf_lst.extend(used_tf.values)
        used_ur_lst.extend(used_ur.values)
        used_uf_lst.extend(used_uf.values)

    lock = Lock()
    lock_failed = Lock()
    thread = opt.thread if len(data) > opt.thread else len(data)
    logging.info(f"max threads: {opt.thread}")
    logging.info(f"true threads_count={thread}")
    p = Pool(processes=thread, initializer=gl.init_lock, initargs=(lock, lock_failed,))
    snp_flag = 1 if opt.fasta else 0
    for pos in data:
        if pos[5]==0:
            #print(pos)
            continue
        pos_idx = 1 if opt.main_chain == "CT" in pos[1] else 2  # for methylation, 20220319
        """
        designer(pos, opt.main_chain, opt.len_short, opt.len_long, opt.extend_range, opt.tm_min, opt.tm_max,
                            opt.extend, opt.amp_min, opt.amp_max, opt.software["primersDB"], opt.path,
                            opt.snp_dir, opt.dis_ud_min, opt.dis_ud_max, opt.snp_flag, opt.snGSP2tar,
                            opt.window_size, pos_idx, res_lst, ref_dic, used_d_lst, used_uf_lst, used_tr_lst,
                            used_ur_lst, used_tf_lst, failed_lst)
        """
        p.apply_async(func=designer,
                      args=(pos, opt.main_chain, opt.len_short, opt.len_long, opt.extend_range, opt.tm_min, opt.tm_max,
                            opt.extend, opt.amp_min, opt.amp_max, opt.software["primersDB"], opt.path,
                            opt.snp_dir, opt.dis_ud_min, opt.dis_ud_max, opt.snp_flag, opt.snGSP2tar,
                            opt.window_size, pos_idx, res_lst, ref_dic, used_d_lst, used_uf_lst, used_tr_lst,
                            used_ur_lst, used_tf_lst, failed_lst))

    p.close()
    p.join()
    #"""
    failed = list(failed_lst)
    if len(res_lst) == 0 and len(failed) > 0:
        logging.info("All target region design failed, can't find a qualified primer nearby the targets' region.")
        sys.stderr.write("All target region design failed, can't find a qualified primer nearby the targets' region.\n")
        sys.exit()
    elif len(res_lst) == 0 and len(failed) == 0:
        logging.info("All target region design failed because of some program error.")
        sys.stderr.write("All target region design failed because of some program error.\n")
        sys.exit()
    result = pd.DataFrame(data=list(res_lst), columns=['name_id', 'chr', 'Product GC percent', 'product size', 'GSP ID',
                                                       'GSP seq', 'GSP start', 'GSP end', 'GSP strand', 'GSP score',
                                                       'snGSP ID', 'snGSP seq', 'snGSP start', 'snGSP end', 'snGSP strand',
                                                       'snGSP score', 'tagGSP ID', 'tagGSP seq', 'tagGSP start', 'tagGSP end',
                                                       'tagGSP strand', 'tagGSP score', 'GSP snp', 'snGSP snp', 'tagGSP snp'])
    result['amp_start'] = np.where(result['snGSP strand'] == '+', result['snGSP start'], result['tagGSP start'])
    result['amp_end'] = np.where(result['snGSP strand'] == '+', result['tagGSP end'], result['snGSP end'])
    result.sort_values(by=['chr', 'amp_start', 'name_id'], ascending=(True, True, True), inplace=True)
    result_CT = result[result['chr'].str.endswith('_CT')]
    result_GA = result[result['chr'].str.endswith('_GA')]
    col1 = ['name_id', 'Product GC percent', 'product size', 'chr', 'snGSP seq', 'snGSP start', 'snGSP end', 'snGSP strand', 'chr', 'tagGSP seq', 'tagGSP start', 'tagGSP end', 'tagGSP strand']
    col2 = ['name_id', 'amp. GC percent(%)', 'amp. size(bp)', 'chr_name1', 'primer1', 'start1', 'end1', 'strand1', 'chr_name2', 'primer2', 'start2', 'end2', 'strand2']

    output_CT, output_GA = result_CT[col1], result_GA[col1]
    output_CT['snGSP seq'] = output_CT['snGSP seq'].apply(lambda x: x[21:])
    output_GA['snGSP seq'] = output_GA['snGSP seq'].apply(lambda x: x[21:])
    output_CT['tagGSP seq'] = output_CT['tagGSP seq'].apply(lambda x: x[17:])
    output_GA['tagGSP seq'] = output_GA['tagGSP seq'].apply(lambda x:x[17:])
    output_CT.columns = col2
    output_GA.columns = col2
    with pd.ExcelWriter(opt.path + '/' + opt.output + '.xlsx', mode='w', engine='openpyxl') \
            as writer1, pd.ExcelWriter(opt.path + '/Output.xlsx', mode='w', engine='openpyxl') as writer2:
        result_CT.to_excel(writer1, sheet_name='chain_CT_primer', index=False)
        result_GA.to_excel(writer1, sheet_name='chain_GA_primer', index=False)
        output_CT.to_excel(writer2, sheet_name='chain_CT_primer', index=False)
        output_GA.to_excel(writer2, sheet_name='chain_GA_primer', index=False)

    # result.to_excel(opt.path + '/' + opt.output + '.xlsx', index=False)
    result_CT[['chr', 'amp_start', 'amp_end', 'name_id']].to_csv(opt.path + '/amplicon_region_CT.bed', sep='\t', index=False,
                                                              header=False)
    result_GA[['chr', 'amp_start', 'amp_end', 'name_id']].to_csv(opt.path+'/amplicon_region_GA.bed', sep='\t', index=False, header=False)
    os.system(f"cat {opt.path}/amplicon_region_CT.bed {opt.path}/amplicon_region_GA.bed >{opt.path}/amplicon_region.bed")

    # failed_cnt, failed_len, failed_df1 = LoadData.cross_fail_region(opt.path, failed, merge_dic)
    design_summary_log(result, result_CT, result_GA, len(failed), opt)
    failed_df = get_uncover_region(opt.path, failed)

    return result, failed_df


def designer(pos, main_chain, len_short, len_long, extend_range, tm_min, tm_max, extend, amp_min, amp_max, script_dir, path,
             snp_dir, dis_ud_min, dis_ud_max, snp_flag, snGSP2tar, window_size, pos_idx,
             res_lst, ref_dic, used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic, used_tf_dic, failed_lst):
    """"""
    if pos[-1] > 0:
        pos[0] = f'Homo_{pos[0]}'
    # snp_dic = {} if fasta_flag else LoadData.load_snp(snp_dir, extend_range, pos[1], pos[2], pos[3])
    snp_dic = PreD.load_snp(snp_dir, extend_range, pos[1], pos[2], pos[3])  # for methylation, 20220319
    if pos[3] - pos[2] < 1.25 * window_size:
        res, flag, pos = point_design(len_short, len_long, extend_range, tm_min, tm_max, extend, amp_min, amp_max,
                                      script_dir, path, dis_ud_min, dis_ud_max, snp_flag, snGSP2tar,
                                      pos, snp_dic, ref_dic, pos_idx, 1,
                                      res_lst, used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic, used_tf_dic,failed_lst)
    else:
        res = region_design(main_chain, len_short, len_long, extend_range, tm_min, tm_max, extend, amp_min, amp_max,
                            script_dir, path, dis_ud_min, dis_ud_max, snp_flag, snGSP2tar, window_size,
                            pos, snp_dic, ref_dic, res_lst, used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic, used_tf_dic,
                            failed_lst)


def point_design(len_short, len_long, extend_range, tm_min, tm_max, extend, amp_min, amp_max, script_dir, path,
                 dis_ud_min, dis_ud_max, snp_flag, snGSP2tar, position, snp_dic, ref_dic, flag, count, res_lst,
                 used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic, used_tf_dic, failed_lst, cut_cnt=3000):
    """短位点/区域引物设计"""
    f_primer, r_primer, ref = split_seq(len_short, len_long, extend_range, tm_min, tm_max,
                                                 script_dir, path, position, snp_dic, ref_dic)
    # print("primer count: ", len(f_primer), len(r_primer))
    n_f_primer = f_primer[:100] if len(f_primer) > 100 else f_primer
    n_r_primer = r_primer[:100] if len(r_primer) > 100 else r_primer
    snGSP_tagGSP = CombinePrimer.build_snGSP_tagGSP(amp_min, amp_max, n_f_primer, n_r_primer)

    dt_len, is_cut = len(snGSP_tagGSP), 0
    # print("primer count: ", len(f_primer), len(r_primer), dt_len)
    if dt_len < 2000 and (len(f_primer) > 100 or len(r_primer) > 100):
        snGSP_tagGSP = CombinePrimer.build_snGSP_tagGSP(amp_min, amp_max, f_primer, r_primer)
        dt_len = len(snGSP_tagGSP)
    # print("fsp-rsp-snGSP_tagGSP count: ", len(f_primer), len(r_primer), dt_len)
    log_content = 'design'
    if dt_len == 0 and extend == 0:
        # gl.lock_fail.acquire()
        failed_lst.append([f'{position[0]}_{count}'] + list(position[1:]) + ['high/low GC low complicated region'])
        logging.info(f"snGSP_tagGSP doesn't enough, F{len(f_primer)}-R{len(r_primer)}-P{dt_len} {[f'{position[0]}_{count}'] + list(position[1:])} failed")
        # gl.lock_fail.release()
        return 'not ok', flag, None
    elif dt_len > cut_cnt:
        res, flag, pos = CombinePrimer.GSP_snGSP_tagGSP(dis_ud_min, dis_ud_max, amp_max, snp_flag, snGSP2tar, path,
                                                   ref, position[0], count, f_primer, r_primer,
                                                   snGSP_tagGSP[:cut_cnt], flag, position[2], position[3], res_lst,
                                                   used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic, used_tf_dic)
        is_cut = 1
    else:
        res, flag, pos = CombinePrimer.GSP_snGSP_tagGSP(dis_ud_min, dis_ud_max, amp_max, snp_flag, snGSP2tar, path,
                                                   ref, position[0], count, f_primer, r_primer,
                                                   snGSP_tagGSP, flag, position[2], position[3], res_lst,
                                                   used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic, used_tf_dic)
    if res == 'not ok' and is_cut:
        # 查找剩余组合
        log_content += ', find more...'
        new_snGSP_tagGSP = snGSP_tagGSP[cut_cnt:]
        # random.shuffle(new_snGSP_tagGSP)
        res, flag, pos = CombinePrimer.GSP_snGSP_tagGSP(dis_ud_min, dis_ud_max, amp_max, snp_flag, snGSP2tar, path,
                                                   ref, position[0], count, f_primer, r_primer,
                                                   new_snGSP_tagGSP, flag, position[2], position[3],
                                                   res_lst, used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic,
                                                   used_tf_dic)

    # gl.lock_fail.acquire()
    if res == 'not ok':
        if "find more" in log_content:
            failed_lst.append([f'{position[0]}_{count}'] + list(position[1:]) + ["primers' conflict"])
        else:
            failed_lst.append([f'{position[0]}_{count}'] + list(position[1:]) + ["high/low GC low complicated region"])
        logging.info(f"F{len(f_primer)}-R{len(r_primer)}-P{dt_len} {[f'{position[0]}_{count}'] + list(position[1:])} {log_content} failed.")
    else:
        logging.info(f"{[f'{position[0]}_{count}'] + list(position[1:])} {log_content} succeed.")
    # gl.lock_fail.release()

    return res, flag, pos


def region_design(main_chain, len_short, len_long, extend_range, tm_min, tm_max, extend, amp_min, amp_max, script_dir,
                  path, dis_ud_min, dis_ud_max, snp_flag, snGSP2tar, window_size, position, snp_dic, ref_dic, res_lst,
                  used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic, used_tf_dic, failed_lst):
    """连续区域引物设计"""
    flag ,count = 0, 0
    start, end, chr_name = position[2], position[3], position[1].split('_')[0]
    sub_pos = position
    sub_pos[1] = chr_name + '_CT' if main_chain == 'CT' else chr_name + '_GA'
    while start < end:
        flag += 1
        count += 1
        sub_pos[2] = start
        if start + 1.25 * window_size < end:
            sub_pos[3] = start + window_size
            print("for design pos: ",sub_pos)
            res, flag, pos = point_design(len_short, len_long, extend_range, tm_min, tm_max, extend, amp_min, amp_max,
                                          script_dir, path, dis_ud_min, dis_ud_max, snp_flag, snGSP2tar,
                                          sub_pos, snp_dic, ref_dic, flag, count,
                                          res_lst, used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic, used_tf_dic, failed_lst)
            if res == 'not ok':
                start += window_size
                sub_pos[1] = chr_name + '_GA' if sub_pos[1].split('_')[-1] == 'CT' else chr_name + '_CT'   # 2026-01-04, the last failed, need to change the ref chain next.
                print("failed, change after: ",sub_pos)
                continue
            else:
                start = pos
            sub_pos[1] = chr_name + '_GA' if sub_pos[1].split('_')[-1] == 'CT' else chr_name + '_CT'
            print("change after: ",sub_pos)
        else:
            print("change before: ",sub_pos)
            sub_pos[3] = end
            #sub_pos[1] = chr_name + '_GA' if sub_pos[1].split('_')[-1] == 'CT' else chr_name + '_CT'
            print("change after: ",sub_pos)
            res, flag, pos = point_design(len_short, len_long, extend_range, tm_min, tm_max, extend, amp_min, amp_max,
                                          script_dir, path, dis_ud_min, dis_ud_max, snp_flag, snGSP2tar,
                                          sub_pos, snp_dic, ref_dic, flag, count,
                                          res_lst, used_d_dic, used_uf_dic, used_tr_dic, used_ur_dic, used_tf_dic, failed_lst)
            start = end + 1

    return res
