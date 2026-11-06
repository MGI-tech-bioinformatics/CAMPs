"""
@File   : PrepareData.py
@Time   : 2023/10/26 14:57
@Author : Tian Chen
@Email  : chentian@genomics.cn
"""
import os
import sys
import time
import gzip
import logging
from functools import wraps
from collections import defaultdict
from multiprocessing import Manager

import pandas as pd

import globalVar as gl
import Similarity as Sim


def time_this(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        r = func(*args, **kwargs)
        end = time.perf_counter()
        logging.info(f"{func.__module__}.{func.__name__}: {round(end - start, 2)}")
        return r
    return wrapper

def load_reference(ref_path):
    try:
        if ref_path[-3:] == '.gz':
            with gzip.open(ref_path, 'r') as f:
                ref = f.read()
                ref_seq = ref.decode('utf-8').split('>')
        else:
            with open(ref_path, 'r') as f:
                ref = f.read()
                ref_seq = ref.split('>')
    except:
        logging.info(f"Load the {ref_path} reference failed, please check your reference's name and path.")
        sys.stderr.write(f"Load the {ref_path} reference failed, please check your reference's name and path.\n")
        exit()
    ref_dic = Manager().dict()
    ref_len_dic = dict()
    for seq in ref_seq[1:]:
        id_seq = seq.split('\n', 1)
        name = id_seq[0].split()[0]
        ref_dic[name] = id_seq[1].replace('\n', '')
        ref_len_dic[name] = len(ref_dic[name])
    # print(ref_dic.keys())
    return ref_dic, ref_len_dic

def check_input_target(opt, ref_len_dic):
    ''' 对用户输入的数据进行校验'''
    pos_data = pd.DataFrame()
    if opt.pos_type == 'cg_number':
        pos_data = get_cg_number_pos(opt)
    elif opt.pos_type == 'position':
        pos_data = get_input_pos(opt, ref_len_dic)

    # print(pos_data.head())
    pos_data, merge_dic = merge_pos_gap(opt, 0, pos_data)
    pos_data = pd.DataFrame(data=pos_data, columns=['name', 'chr', 'start', 'end', 'length'])
    pos_data = pos_data[['chr', 'start', 'end', 'name']]
    pos_data.to_csv(opt.path+'/Input.bed', sep='\t', header=None, index=False)


def get_cg_number_pos(opt):
    """根据输入的cg编号名称，检索对应的cg位点坐标信息"""
    try:
        with open(opt.infile.name, 'r') as f:
            cg_number = f.readlines()
    except UnicodeDecodeError:
        logging.info("load target cg_name failed.")
        sys.stderr.write("Load target cg_name failed. Invalid file formats, please check your input.\n")
        sys.exit()
    except IOError:
        logging.info("load target cg_name list failed.")
        sys.stderr.write("Load target cg_name failed! Invalid input formats may exist, please check your input.\n")
        sys.exit()
    cg_lst = [x.strip().lower() for x in set(cg_number)]
    res = []
    with open(opt.cg_index, 'r') as f:
        while True:
            line = f.readline().strip().split()
            if line[-1] in cg_lst:
                res.append(line)
                cg_lst.remove(line[-1])
            if len(cg_lst) == 0:
                break
    if len(res) < 0 :
        logging.error("can't find all target cg number in database.")
        sys.stderr.write("None of the input cg numbers were found in database, please check your input and"
                         " reference option.\n")
        sys.exit()
    if len(cg_lst) > 0 :
        logging.warning(f"{cg_lst} can't find in cg_index")
        sys.stderr.write(f"Warnings: The {cg_lst} can't be found in database and have been ignored when designing.\n")
    result = pd.DataFrame(data=res, columns=['chr', 'start', 'end', 'name'])

    # result.drop_duplicates(inplace=True)
    # result.sort_values(by=['chr', 'start', 'end', 'name'], ascending=(True, True, True, True), inplace=True)
    return result


def get_input_pos(opt, ref_len_dic):
    try:
        with open(opt.infile.name, 'r') as f:
            data = f.readlines()
    except UnicodeDecodeError:
        logging.info("load target position failed.")
        sys.stderr.write("Load target region decode failed, Invalid file formats. please check your input.\n")
        sys.exit()
    except IOError:
        logging.info("load target position failed.")
        sys.stderr.write("Load target region failed. The input file not found, please check your input.\n")
        sys.exit()
    name_set, res = set(), []
    rename, rename_id = defaultdict(list), 0
    try:
        for idx, line in enumerate(data):
            if ',' in line:
                pos = line.strip().split(',')
            elif '\t' in line:
                pos = line.strip().split('\t')
            else:
                pos = line.strip().split()
            # pos = re.split('\s+|,', line.strip())
            if not pos[0]:
                continue
            # pos = line.strip().split()
            if len(pos) == 0 or pos[0][0] == '#':
                continue
            try:
                pos[1], pos[2] = int(pos[1]), int(pos[2])
            except:
                if idx != 0:
                    logging.info(f"{pos} position error")
                    sys.stderr.write(
                        f"The target region {pos} meet position error, start or end not a numeric type,"
                        f"please check the position coordinates information.\n")
                continue
            print(pos)
            if len(pos) == 3:
                rename_id += 1
                sys.stdout.write(f"The line {idx + 1}:{pos} missing the target name, rename as target.{rename_id}.\n")
                pos.append(f'target.{rename_id}')

            if len(pos) < 3:
                logging.info(f"The {idx + 1} input line unqualified")
                sys.stderr.write(f"The {idx + 1} input line unqualified, please check it.")
                sys.exit()
            if opt.ref_name in ['hg19', 'hg38']:
                pos[0] = 'chr' + pos[0] if 'chr' not in pos[0].lower() else pos[0][:-1].lower() + pos[0][-1].upper()
            if pos[0].upper() + "_CT" in ref_len_dic.keys() and pos[2] <= ref_len_dic[pos[0].upper() + "_CT"]: # 甲基化BS处理后的参考序列带_CT/_GA后缀
                if int(pos[2]) - int(pos[1]) < 0:
                    logging.info(f"{pos} position error")
                    sys.stderr.write(f"The target region {pos} meet position error, the start coordinate is greater"
                                     f" than the end coordinate, please check the position coordinates information.\n")
                    continue
                if pos[3] not in name_set:
                    name_set.add(pos[3])
                else:
                    rename[pos[3]].append(1)
                    new_name = pos[3] + '.' + str(len(rename[pos[3]]))
                    logging.info(f"The line {idx + 1} target name '{pos[3]}' not unique or the target region missing"
                                 f" name, rename as {new_name}.")
                    sys.stdout.write(f"The line {idx + 1} target name '{pos[3]}' not unique or the target region missing"
                                     f" name, rename as {new_name}.\n")
                    pos[3] = new_name
            else:
                logging.info(f"The line {idx + 1} target region {pos} can't find on {opt.ref_name} reference, please check"
                             f" the chromosome index.")
                sys.stderr.write(f"The line {idx + 1} target region {pos} can't find in {opt.ref_name} reference,"
                                 f" please check the chromosome name and the position coordinates.\n")
                continue

            res.append(pos[:4])
    except:
        logging.info("input target position error or get wrong type.")
        sys.stderr.write("Load the input target region meet error, please check your input type and reference "
                         "option.\n")
        sys.exit(1)
    if len(res) == 0:
        sys.stderr.write("All inputs are invalid, please check it.\n")
        sys.exit()
    result = pd.DataFrame(data=res, columns=['chr', 'start', 'end', 'name'])

    return result


def merge_pos_gap(opt, min_gap, *args):
    data = pd.DataFrame()
    if args:
        for col in args:
            data = col[['name', 'chr', 'start', 'end']]
    else:
        try:
            data = pd.read_csv(opt.path+'/Input.bed', sep='\s+', header=None, names=['chr', 'start', 'end', 'name'])
            data = data[['name', 'chr', 'start', 'end']]
        except IOError as e:
            logging.error(f"load checked input data failed, please check it.")
            sys.stderr.write(f"load checked input data failed, please check it.")
            sys.exit()

    data.name = data.name.astype(str)
    data.start = data.start.astype(int)
    data.end = data.end.astype(int)
    data = data.sort_values(by=['chr', 'start'], ascending=(True, True))
    # 考虑输入错误的情况，start > end, or use try ... except
    data['length'] = data['end'] - data['start'] + 1
    data_lst = data.values
    head = 'merge-name\tgene-name\tchr\tstart\tend\tlen\n'
    write_merge_file(opt.path, 'w', head)
    merge_flag, merge_index = 1, 0
    new_pos, merge_dic = [], {}
    last = data_lst[0]
    merge_dic.setdefault(last[1].upper()+'_CT', []).extend([x for x in range(last[2], last[3] + 1)])
    merge_dic.setdefault(last[1].upper()+'_GA', []).extend([x for x in range(last[2], last[3] + 1)])
    # format: name | chr | start | end | length
    merge_length = 0
    for line in data_lst[1:]:
        merge_dic.setdefault(line[1].upper()+'_CT', []).extend([x for x in range(line[2], line[3] + 1)])
        merge_dic.setdefault(line[1].upper()+'_GA', []).extend([x for x in range(line[2], line[3] + 1)])
        if line[1] == last[1] and line[2] - last[3] <= min_gap:
            last[0] = f"{last[0]}-{line[0]}"
            if line[3] > last[3]:  # 考虑有些目标区域存在包含关系
                last[3] = line[3]
                last[4] += line[4]  # 目标长度累加
            merge_flag += 1
            continue
        else:
            merge_index, last = merge(opt.path, merge_flag, merge_index, last, min_gap)
            new_pos.append(last)
            merge_flag = 1
            merge_length += last[-1]
        last = line
    merge_index, last = merge(opt.path, merge_flag, merge_index, last, min_gap)
    # merge_dic.setdefault(last[1], []).extend([x for x in range(last[2], last[3] + 1)])
    new_pos.append(last)
    merge_length += last[-1]
    merge_summary_log(len(data), data['length'].sum(), len(new_pos), merge_length)
    if min_gap > 0:
        with open(opt.path + '/predict_time.txt', 'w') as f:
            f.write(f"{predict_time(opt, merge_length, len(new_pos))}")

    return new_pos, merge_dic


def write_merge_file(path, mode, content):
    """将合并后的区域位点索引写入merge_name_index.txt 文件"""
    try:
        with open(path+'/merge_name_index.txt', mode) as f:
            f.write(content)
    except IOError as e:
        logging.info(f"can't write to merge_name_index.txt")
        raise e


def merge(path, merge_flag, merge_index, last, min_gap):
    """merge操作后，目标区域名称转换, min_gap：两个邻近区域之间容忍的最小距离"""
    if merge_flag != 1:
        merge_index += 1
        content = f"Merge{str(merge_index)}_gap{min_gap}_{str(merge_flag)}\t{last[0]}\t{last[1]}\t" \
                  f"{str(last[2])}\t{str(last[-2])}\t{str(last[-1])}\n"
        write_merge_file(path, 'a+', content)
        last[0] =  f"Merge{str(merge_index)}_gap{min_gap}_{str(merge_flag)}"

    return merge_index, last

def merge_summary_log(cnt1, len1, cnt2, len2):
    """输出merge前后，目标位点及长度信息"""
    logging.info(f"\n"
                 f"original targets region number: {cnt1}\n"
                 f"original targets region total length: {len1} bp\n"
                 f"after merge targets region number: {cnt2}\n"
                 f"after merge targets region total length: {len2} bp\n"
                 f"{'=='*30}")
    # with open(opt.path+'/predict_time.txt', 'w') as f:
    #     f.write(f"{predict_time(opt, len2, cnt2)}")
    gl.report_dic['originTargetCnt'] = cnt1
    gl.report_dic['originTargetLen'] = len1
    gl.report_dic['mergeTargetCnt'] = cnt2
    gl.report_dic['mergeTargetLen'] = len2


def predict_time(opt, total_len, new_pos_cnt):
    # According the target region after merge to predict the design time.
    amp_size = opt.amp_max if total_len / new_pos_cnt < opt.amp_max * 2 else opt.amp_min
    predict_plex = (total_len / new_pos_cnt / amp_size + 1.5) * new_pos_cnt
    # threads = new_pos_cnt if new_pos_cnt < opt.thread else opt.thread
    # aln_t = round((total_len + new_pos_cnt * 400) / 8, 2)
    logging.info(f"predict_plex_nums: {predict_plex}")

    if total_len > 400000 or new_pos_cnt > 2000:
        logging.info("Your target regions' length over 400kbp or target counts over 2000")
        sys.stderr.write(
            "Your target regions' length over 400kbp or target counts over 2000, CAMPs currently supports at most 400 kbp of target regions and 2000 targets.")
        sys.exit()
    if predict_plex < 80:
        return "1 hour"
    elif predict_plex < 300:
        return "0.5 days"
    elif predict_plex < 800:
        return "1 day"
    elif predict_plex < 1500:
        return "1.5 days"
    elif predict_plex < 2000:
        return "2 days"
    else:
        return ">2 days"


def save_merge_pos(opt, merge_lst):
    with open(opt.path + '/merge_pos.csv', 'w') as f, open(opt.path + f'/merge_pos_CT.csv', 'w') as f1, open(opt.path+f'/merge_pos_GA.csv', 'w') as f2:
        f.write(",".join(['name', 'chr', 'start', 'end', 'length']) + "\n")
        f1.write(",".join(['name', 'chr', 'start', 'end', 'length', 'design_flag']) + "\n")
        # f2.write(",".join(['name', 'chr', 'start', 'end', 'length']) + "\n")
        design_CT = 1 if opt.cov == 2 or opt.main_chain == 'CT' else 0
        design_GA = 1 if opt.cov == 2 or opt.main_chain != 'CT' else 0
        for row in merge_lst:
            row = [str(x) for x in row]
            f.write(",".join(row) + "\n")
            #if int(row[-1]) < 1.25 * opt.window_size:
            if int(row[3]) - int(row[2]) + 1 < 1.25 * opt.window_size: # 2026-01-04, need to uset merge_length.
                f1.write(f"{row[0]}_CT,{row[1].upper()}_CT,{','.join(row[2:])},{design_CT}\n")
                f2.write(f"{row[0]}_GA,{row[1].upper()}_GA,{','.join(row[2:])},{design_GA}\n")
            else:
                region_CT = 1 if opt.main_chain == 'CT' else 0
                region_GA = 0 if opt.main_chain == 'CT' else 1
                f1.write(f"{row[0]}_CT,{row[1].upper()}_CT,{','.join(row[2:])},{region_CT}\n")
                f2.write(f"{row[0]}_GA,{row[1].upper()}_GA,{','.join(row[2:])},{region_GA}\n")
    try:
        os.system(f'cat {opt.path}/merge_pos_CT.csv {opt.path}/merge_pos_GA.csv >{opt.path}/merge_dup_pos.csv')
    except IOError as e:
        logging.error(f'loading merge_pos file meet error,please check it.')
        raise  e


def load_pos_csv(path, file_name):
    pos_df = pd.DataFrame()
    try:
        if '.csv' in file_name:
            pos_df = pd.read_csv(path + f'/{file_name}')
        elif '.txt' in file_name:
            pos_df = pd.read_csv(path + f'/{file_name}', sep='\t')
    except IOError as e:
        logging.info("can't load the target position after merge.")
        # sys.stderr.write("")
        raise e

    return pos_df.values


def homo_detect(opt, ref_dic, merge_dic, t=8):
    """对合并后的区域进行同源检测"""
    pos_lst = load_pos_csv(opt.path, 'merge_dup_pos.csv')
    if not os.path.exists(f'{opt.path}/homo_info.xlsx'):
        fa = open(opt.path + '/Product_pos.fa', 'w')
        for pos in pos_lst:
            name, chr_id, start, end, design_flag = pos[0], pos[1], int(pos[2]), int(pos[3]), int(pos[5])
            length = end - start + 1
            if int(pos[-1]) < 1.25 * opt.window_size:  # 同源评估与设计时的window_size 保持一致
                write_fa_file(name, chr_id, start, end, length, design_flag, opt.extend_range, ref_dic, fa)
            else:
                end = start + opt.homo_window_size
                while end < pos[3]:
                    write_fa_file(name, chr_id, start, end, opt.window_size, design_flag, opt.extend_range, ref_dic, fa)
                    start = end
                    end = start + opt.homo_window_size
                end = pos[3]
                length = end - start + 1
                write_fa_file(name, chr_id, start, end, length, design_flag, opt.extend_range, ref_dic, fa)
        fa.close()
        cmd1 = f"{opt.software['bwa']} mem {opt.ref_path} {opt.path}/Product_pos.fa -t {t} -a 1>{opt.path}/Product_aln.sam " \
               f"2>>{opt.path}/logs/bwa.log"
        os.system(cmd1)
        os.system(f'grep "^@" {opt.path}/Product_aln.sam >{opt.path}/sam_head')
        all_pos = Sim.HomoSim(opt.path + '/Product_aln.sam', merge_dic)
        all_pos.to_excel(opt.path + '/homo_info.xlsx', index=False)
    else:
        all_pos = pd.read_excel(opt.path + '/homo_info.xlsx')
    region = merge_un_homo(all_pos)
    # region.to_excel(prj_path+'/homo_info_merge.xlsx', index=False)
    tar_region = region[region['sim_cnt'] < opt.homo_num]
    tar_region.drop(['sim_region'], axis=1, inplace=True)
    homo_0, homo_1 = homo_summary_log(opt.path, region, tar_region, max_homo_cnt=opt.homo_num)
    return tar_region, homo_0.values.tolist(), homo_1


def write_fa_file(name, chr_id, start, end, length, flag, extend_range, ref_dic, fa):
    """构建fasta文件"""
    name_id = f'>{name}:{chr_id}:{start}:{end}:{length}:{flag}'
    seq_start = start - extend_range - 1 if start - extend_range - 1 >= 0 else 0
    seq_end = end + extend_range - 1 if end + extend_range - 1 < len(ref_dic[chr_id]) else None
    seq = ref_dic[chr_id][seq_start: seq_end]
    fa.write(f"{name_id}\n{seq}\n")


def merge_un_homo(region_lst):
    """
        【合并非同源区域】函数
        :parameter region_lst: 同源检测后的位点信息
        :return: 新的合并位点信息
        2019-11-7, chen tian
        """
    data = region_lst.values
    last, cnt = data[0], 0
    new_pos = []
    for row in data[1:]:
        if row[0] == last[0]:
            if row[-3] == last[-3] and row[2] == last[3]:
                last[3] = row[3]
                last[-2] += (row[-2] - 1)
                # last[-1] = np.append(last[-1], row[-1])
            else:
                cnt += 1
                last[0] += f'_sub{cnt}'
                new_pos.append(last)
                last = row
        else:
            if cnt > 0:
                last[0] = last[0] + f'_sub{cnt + 1}'
            new_pos.append(last)
            cnt = 0
            last = row
    if cnt > 0:
        last[0] = last[0] + f'_sub{cnt + 1}'
    new_pos.append(last)
    design_region = pd.DataFrame(data=new_pos, columns=['name', 'chr', 'start', 'end', 'length', 'design_flag', 'sim',
                                                        'sim_cnt', 'tar_len', 'sim_region'])
    design_region['length'] = design_region['tar_len']
    del design_region['tar_len']

    return design_region


def homo_summary_log(prj_path, all_region, design_region, max_homo_cnt=10):
    """ 同源检测结果信息输出及记录"""
    homo_1 = all_region[all_region['sim_cnt'] > 0]
    homo_cnt_1 = homo_1['length'].sum()
    homo_0 = all_region[all_region['sim_cnt'] >= max_homo_cnt]
    homo_cnt_0 = homo_0['length'].sum()
    homo_ratio = round(homo_cnt_1 / all_region['length'].sum() * 100, 2)
    abandon_ratio = round(homo_cnt_0 / all_region['length'].sum() * 100, 2)

    logging.info(f"\n"
                 f"Homology targets region count: {homo_1.shape[0]}\n"
                 f"Homology targets region total length: {homo_cnt_1} bp\n"
                 f"Homology targets region ratio: {homo_ratio}%\n"
                 f"Abandoned Homology region count: {homo_0.shape[0]}\n"
                 f"Abandoned Homology region length: {homo_cnt_0} bp\n"
                 f"Abandoned Homology region ratio: {abandon_ratio} %\n"
                 f"{'=='*30}\n"
                 f"Design region number: {design_region.shape[0]}\n"
                 f"Design region total length: {design_region['length'].sum()} bp\n"
                 f"{'=='*30}")

    gl.report_dic['homoRegionCnt'] = homo_1.shape[0]
    gl.report_dic['homoRegionLen'] = homo_cnt_1
    gl.report_dic['homoRegionRatio'] = homo_ratio
    gl.report_dic['abandonHomoCnt'] = homo_0.shape[0]
    gl.report_dic['abandonHomoLen'] = homo_cnt_0
    gl.report_dic['abandonHomoRation'] = abandon_ratio
    gl.report_dic['designRegionCnt'] = design_region.shape[0]
    gl.report_dic['designRegionLen'] = design_region['length'].sum()

    homo_0['failed reason'] = f'≥{max_homo_cnt} homology regions'
    homo_0 = homo_0[['name', 'chr', 'start', 'end', 'length', 'sim', 'sim_cnt', 'failed reason']]
    homo_0.to_csv(prj_path + '/Abandon_Homo_pos.txt', sep='\t',index=False)
    # design_resion.to_csv(prj_path + '/new_pos_tag.txt', sep='\t', index=False)
    return homo_0, homo_1


def aln_sequence(opt):
    """获取候选短序列比对的特异性信息"""
    [run_exe, t_v, para_K] = ['primersDB_CL', 8, ' -K '+opt.software['primersDB']+'/kernel.cl'] if opt.kernel == 'GPU' else ['primersDB', opt.thread, '']
    pos_lst = load_pos_csv(opt.path, 'new_pos_tag.txt')
    generate_aln_pos(opt, pos_lst)
    is_write = 'B' if os.path.exists(f'{os.path.dirname(opt.ref_path)}/{opt.ref_name}_trans_{opt.k_mer}mer.primerdb') else 'w'
    command = f'{opt.software["primersDB"]}{run_exe} {opt.ref_path}{para_K}' \
              f' -r {opt.path}/aln_pos.txt' \
              f' -{is_write} {os.path.dirname(opt.ref_path)}/{opt.ref_name}_trans_{opt.k_mer}mer.primerdb' \
              f' -k {opt.k_mer}' \
              f' -t {t_v} -pp -3 -c' \
              f' -z {opt.path}/Hits_result 2>{opt.path}/logs/primersDB.log'
    logging.info(f"run primersDB as follows command:\n{command}")
    os.system(command)
    with open(opt.path+'/Hits_result.ln.gz.OK', 'w') as f:
        f.write(f'{time.asctime(time.localtime(time.time()))}')
    return True


def generate_aln_pos(opt, tar_lst):
    """生成短序列候选区间"""
    try:
        with open(opt.path+'/aln_pos.txt', 'w') as fw:
            for line in tar_lst:
                start = line[2] - opt.extend_range if line[2] - opt.extend_range > 0 else 0
                end = line[3] + opt.extend_range
                fw.write(f"{line[1]}\t{start}-{end}\n")
    except IOError as e:
        logging.info("create aln_pos failed")
        sys.exit()


def load_snp(snp_dir, extend_range, chr_name, start, end):
    """
    load the chr snp info
    :param extend_range:
    :param snp_dir:
    :param chr_name: the target chr
    :param start: target start
    :param end: target end
    :return: chr_snp _info_list
    """
    chr_name = chr_name.split('_')[0].lower()  # for methylation, 20220319
    if chr_name[-1] in ['x', 'y', 'm']:
        chr_name = chr_name[:-1]+chr_name[-1].upper()
    file = snp_dir + chr_name + '.snp'
    # file = 'chr1'+'.snp'
    snp_dic = {}
    try:
        # cmd = f'awk' + " '$2>" + str(start - opt.extend_range) + '&& $2<' + str(end + opt.extend_range) + "' " + file
        cmd = f"awk '$2>{start - extend_range} && $2<{end + extend_range}' {file}"
        snp_info = os.popen(cmd).readlines()
    except IOError:
        logging.info(f"can't find {file}")
        snp_info = []
    for i in snp_info:
        k = i.split()
        # v = k[5][1:-1].split(',')[1:]
        v = k[-1].split(',')
        if k[-1] != '-':
            if float(max(v)) >= 0.01:
                # if float(k[-1]) >= 0.01:
                snp_dic[k[0] + ':' + k[1]] = k
    return snp_dic  # return a snp dict

"""
def homo_res_to_target(path, filename):
    try:
        with open(path+'/'+filename, 'r') as f:
            data = f.readlines()
    except IOError as e:
        logging.info(f"No found {filename} file, please check it.")
        raise e
    last = data[1].strip().split()
    tar_name = set()
    for line in data[2:]:
        row = line.strip().split()
"""
