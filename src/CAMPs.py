#!/usr/bin/env python
# -*- coding: utf-8 -*- #
"""CAMPs: multiplex PCR primer design software for DNA methylation."""
import configparser

Program = 'CAMPs'
Author = 'Tian Chen, MGI, China'
AuthorEmail = 'Tian Chen <chentian@mgi-tech.com>'
Date = 'Oct-25-2023'
License = 'Please contact Tian Chen <chentian@mgi-tech.com>'
Version = 'V1.0'

import os
import time
import logging
import argparse
import pandas as pd
from datetime import datetime

import globalVar as gl
import MethPrimerDesign as MPD
import PrimerEvaluate as PE
import htmlReporter as hR
import PrimerAln as PA
import PrepareData


def get_opt():
    '''check and parsing the opts.'''
    parser = argparse.ArgumentParser(prog='CAMPs',
                                     description="CAMPs: a multiplex-PCR primer design program for methylation.",
                                     usage="%(prog)s.py [options] -i target_position.txt -o output_name", )
    parser.add_argument('--version', action='version', version=f'{Program} {Version}')
    parser.add_argument('--config', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'conf.ini'),
                        help='INI configuration file (default: src/conf.ini).')
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'), required=True,
                        help="[Required] Target region for primer design, format: chr\tstart\tend\tname. [file]")
    parser.add_argument('-o', '--output', type=str, required=False, default='test_project',
                        help="[optional] The project name for output. [string]")
    parser.add_argument('--ref_name', type=str, required=False, default='hg19', choices=['hg19', 'hg38', 'contol_ref'],
                        help="[optional] The reference name, default = 'hg19'. [str]")
    parser.add_argument('--pos_type', type=str, required=False, default='position', choices=['position', 'cg_number'],
                        help="[optional] The input target position type, default is 'position'. ['position', 'cg_number']")
    parser.add_argument('--flanking', type=int, required=False, default=0, choices=[0, 5, 10, 15],
                        help="[optional] The flanking value [0, 5, 10, 15] with cg position, default = 0.[int]")
    parser.add_argument('--order_num', type=str, required=False, default=f"CAMPs{datetime.now().date().strftime('%Y%m%d')}{datetime.now().time().strftime('%H%M%S')}",
                        help="[optional] The order number for the primer design, default='CAMPs+Date'.[string]")
    parser.add_argument('--pre_primer', nargs='?', type=argparse.FileType('r'), required=False, default=None,
                        help="[optional] The Primer .xlsx file have been designed, default is None. [file]")
    parser.add_argument('--ref_path', type=str, required=False, default=None,
                        help="[optional] Reference library path. Required for --ref_name contol_ref.")
    parser.add_argument('--snp_dir', type=str, required=False, default=None,
                        help="[optional] The snp lib path, default = hg19 snp information directory .[str]")
    parser.add_argument('--main_chain', type=str, required=False, default='CT', choices=['CT', 'GA'],
                        help="[optional] choose the main chain to design primer for target region at first, choose 'CT' or 'GA chain. default='CT'.[str]")
    parser.add_argument('--cov', type=int, required=False, default=2, choices=[1, 2],
                        help="[optional] The coverage for target region, choose 1: only design primer based on main_chain, choose 2: design primer based on CT and GA chain. default=2.[int]")


    # parser.add_argument('--cg_path', type=str, required=False, default=None,
    #                     help="[optional] The snp lib path, default = hg19 rs number information directory .[str]")

    parser.add_argument('-t', '--thread', type=int, required=False, default=36,
                        help="[optional] The max thread number, default = 36. [int]")
    parser.add_argument('--window_size', type=int, required=False, default=80,
                        help="[optional] The sliding window size in continuous target region. default = 80 bp. [int]")
    parser.add_argument('--homo_window_size', type=int, required=False, default=200,
                        help="[optional] The homology detect sliding window size. default = 200 bp. [int]")
    parser.add_argument('-k', '--k_mer', type=int, required=False, default=10,
                        help="[optional] The k-mer size during short sequence alignment, default = 10. [int]")
    parser.add_argument('--kernel', type=str, required=False, default='GPU', choices=['GPU', 'CPU'],
                        help="[optional] The primersDB running parameter ['GPU', 'CPU'], default = 'GPU'. [str]")
    parser.add_argument('--software', type=str, required=False, default={},
                        help="The software configuration used by CAMPs.")

    filter_g = parser.add_argument_group('Results filter settings', 'set these arguments for filtering the results.')


    filter_g.add_argument('--sim_rate', type=float, required=False, default=0.9,
                          help="[optional] The threshold of homology for target position, default = 0.9. [float]")
    filter_g.add_argument('--homo_num', type=int, required=False, default=10,
                          help="[optional] The threshold of Homology count for target position, default = 10. [int]")
    filter_g.add_argument('--tm_max', type=float, required=False, default=62,
                          help="[optional] The allowed max tm value, default = 62.0. [float]")
    filter_g.add_argument('--tm_min', type=float, required=False, default=58,
                          help="[optional] The allowed min tm value, default = 58.0.  [float]")
    filter_g.add_argument('--poly', type=int, required=False, default=5,
                          help="[optional] The acceptable length of poly, default = 5. [int]")
    filter_g.add_argument('--max_maf', type=float, required=False, default=0.1,
                          help="[optional] The max maf value for snp. default = 0.1. [float]")
    filter_g.add_argument('--gc_max', type=float, required=False, default=60,
                          help="[optional] The primers allowed max gc percent, default = 60. [int]")
    filter_g.add_argument('--gc_min', type=float, required=False, default=35,
                          help="[optional] The primers allowed min gc percent, default = 35. [int]")
    filter_g.add_argument('--extend_range', type=int, required=False, default=200,
                          help="[optional] The candidate primers interval range, default = 200. [int]")
    filter_g.add_argument('--amp_max', type=int, required=False, default=200,
                          help="[optional] The amplicon's max size, default = 200. [int]")
    filter_g.add_argument('--amp_min', type=int, required=False, default=100,
                          help="[optional] The amplicon's min size, default = 150. [int]")
    filter_g.add_argument('-s', '--snp_flag', action='store_true', required=False,
                          help="[optional] The target is SNP position flag, the default is False.")
    filter_g.add_argument('--snGSP2tar', type=int, required=False, default=50,
                          help="[optional] The max distance between snGSP start and target end, default = 50 bp. [int]")
    filter_g.add_argument('--se_strict', type=int, choices=[0, 1], default=1,
                          help="[optional] Whether strict for snGSP2tar, default = 1.")
    filter_g.add_argument('-e', '--extend', action="store_true", required=False,
                          help="[optional] Add alternative amp2 parameters.")
    filter_g.add_argument('--amp_max2', type=int, required=False, default=200,
                          help="[optional] The another amplicon's max size, default = 200. [int]")
    filter_g.add_argument('--amp_min2', type=int, required=False, default=100,
                          help="[optional] The another amplicon's min size, default = 100. [int]")
    filter_g.add_argument('--pos_gap', type=int, required=False, default=None,
                          help="[optional] The merged distance between targets, default is amp_min. [int]")
    filter_g.add_argument('--len_long', type=int, required=False, default=35,
                          help="[optional] The primer sequence's longest length, default = 35 bp. [int]")
    filter_g.add_argument('--len_short', type=int, required=False, default=18,
                          help="[optional] The primer sequence's longest length, default = 18 bp. [int]")
    filter_g.add_argument('--dis_ud_max', type=int, required=False, default=15,
                          help="[optional] The max distance between GSP and snGSP, default = 15 bp. [int]")
    filter_g.add_argument('--dis_ud_min', type=int, required=False, default=-10,
                          help="[optional] The min distance between GSP and snGSP, default = -10 bp. [int]")
    filter_g.add_argument('--amp_type', type=str, required=False, default=None,
                          choices=['60-120', '100-200', '145-165', '200-300', '300-400', "SE50"],
                          help="[optional] The amplicon size, default='60-120'bp. [str]")

    therm_g = parser.add_argument_group('Experimental settings:', 'set these arguments for real experimental values.')
    therm_g.add_argument('--mono', nargs='?', type=float, required=False, default=50,
                         help='[Optional] Concentration of monovalent cations [mM], default = 50.0 [Float]')
    therm_g.add_argument('--diva', nargs='?', type=float, required=False, default=3,
                         help='[Optional] Concentration of divalent cations [mM], default = 3 [Float]')
    therm_g.add_argument('--oligo', nargs='?', type=float, required=False, default=50,
                         help='[Optional] Concentration of annealing oligos [nM], default = 50.0. [Float]')
    therm_g.add_argument('--dntp', nargs='?', type=float, required=False, default=0.25,
                         help='[Optional] Concentration of dNTPs [nM], default = 0.25. [Float]')

    options = parser.parse_args()
    options.scr_dir = os.path.dirname(os.path.abspath(__file__))
    config = configparser.ConfigParser()
    if not config.read(options.config, encoding='utf-8'):
        parser.error(f'Configuration file not found: {options.config}')
    for section in ('ref', 'snp', 'cg', 'software'):
        if not config.has_section(section):
            parser.error(f'Configuration is missing [{section}]')
    for section, keys in {'ref': ('hg19', 'hg38'), 'snp': ('hg19', 'hg38'),
                          'cg': ('hg19', 'hg38'), 'software': ('primersDB', 'bwa', 'samtools')}.items():
        for key in keys:
            if not config.has_option(section, key):
                parser.error(f'Configuration is missing [{section}] {key}')
    options.software = config['software']

    if not options.pos_gap:
        options.pos_gap = options.amp_min if not options.snp_flag else 30 # SE50 default merge gap is 30bp

    if options.ref_name == 'hg19':
        options.ref_path = options.ref_path or config['ref']['hg19']
        options.snp_dir = options.snp_dir or config['snp']['hg19']
        options.cg_index = config['cg']['hg19']
    elif options.ref_name == 'hg38':
        options.ref_path = options.ref_path or config['ref']['hg38']
        options.snp_dir = options.snp_dir or config['snp']['hg38']
        options.cg_index = config['cg']['hg38']
    else:
        if not options.ref_path:
            parser.error('--ref_path is required when --ref_name is contol_ref')
        options.ref_path =  options.ref_path+'/'+ options.ref_name + '.tran.fa'
        options.snp_dir = options.snp_dir or os.path.dirname(options.ref_path)

    if not os.path.isfile(options.ref_path):
        parser.error(f'Reference FASTA not found: {options.ref_path}. Update --config or --ref_path.')
    if options.pos_type == 'cg_number' and not os.path.isfile(getattr(options, 'cg_index', '')):
        parser.error('CpG index not found; configure [cg] for the selected reference.')
    options.snp_dir = os.path.join(options.snp_dir, '')
    options.software = dict(options.software)
    options.software['primersdb'] = os.path.join(options.software['primersdb'], '')
    # Downstream code historically accesses this mixed-case key.
    options.software['primersDB'] = options.software['primersdb']
    options.path = mkdir(options.output)

    if options.amp_type == '60-120':
        options.amp_max, options.amp_min = 120, 60
        options.window_size = 40
        options.extend_range = 110
    elif options.amp_type == '100-200':
        options.amp_max, options.amp_min = 195, 100
        options.window_size = 80
        options.extend_range = 190
    elif options.amp_type == '145-168':
        options.amp_max, options.amp_min = 168, 145
        options.window_size = 70
        options.extend_range = 158
    elif options.amp_type == '200-300':
        options.amp_max, options.amp_min = 295, 200
        options.window_size = 160
        options.extend_range = 280
    elif options.amp_type == '300-400':
        options.amp_max, options.amp_min = 395, 300
        options.window_size = 260
        options.extend_range = 350
    elif options.amp_type == 'SE50':
        options.amp_max, options.amp_min = 90, 60
        options.window_size = 20
        options.extend_range = 80

    return options


def mkdir(project_name):
    """创建项目目录"""
    pwd = os.getcwd()
    path = f'{pwd}/' + project_name
    is_exists = os.path.exists(path)
    if not is_exists:
        os.mkdir(path)
    else:
        print(f"{path} already exists.")
    if not os.path.exists(path+'/logs'):
        os.mkdir(path+'/logs')
        # print(f'{path} created.')
    else:
        print(f"{path}/logs already exists.")

    return path

def bulid_ref_index(path, ref_name):
    ref_path = path + '/data'
    if not os.path.exists(ref_path):
        os.mkdir(ref_path)
        # considered the '\r\n' in gzip file can't be identified, try to decompress the gzip file. 20220-07-29
        if '.gz' in ref_name:
            os.system(f'gzip -dc {ref_name} >{ref_path}/{ref_name[:-3]}')
            os.system(f'bwa index {ref_path}/{ref_name[:-3]} 2>>{path}/logs/bwa.log')
        else:
            os.system(f'ln -s ../../{ref_name} {ref_path}')
            os.system(f'bwa index {ref_path}/{ref_name} 2>>{path}/logs/bwa.log')

    return ref_path + '/', ref_name[:-3] if '.gz' in ref_name else ref_name

def fasta_format_valid(ref_name):
    with open(ref_name, 'r') as f:
        data = f.read().split('>')
    for line in data:
        chr_id, seq = line.split('\n', 1)
        if len(seq) <= 0 :
            logging.error(f"load the reference meet error, please input the correct fasta file as the reference.")
            raise EOFError
    return True


def zip_result(script_path, result_path, order_num):
    pwd = os.getcwd()
    os.system(f'cp {script_path}/README.md {result_path}')
    os.chdir(result_path)
    os.system('cp Input.bed input_pos.bed')
    os.system('cp Output.xlsx output.xlsx')
    os.system(f"zip -q {order_num}.zip  README.md Input.bed amplicon_region*.bed missing_target_*region.bed cover_target_*region.bed Output.xlsx {order_num}.html")
    os.chdir(pwd)



def main():
    gl.init()
    opt = get_opt()
    print("window_size: ",opt.window_size)
    logging.basicConfig(format='%(asctime)s: %(message)s', level=logging.DEBUG,
                        filename=f'{opt.path}/CAMPs.log',
                        filemode='a')
    ref_dic, ref_len_dic = PrepareData.load_reference(opt.ref_path)
    logging.info(f'step0: check input data ...')
    if not os.path.exists(opt.path + '/Input.bed'):
        PrepareData.check_input_target(opt, ref_len_dic)
    #sys.exit(1)

    logging.info(f'step1: Merge adjacent base within  {opt.pos_gap} bp ...')
    if not os.path.exists(opt.path + '/merge_pos.csv'):
        merge_lst, merge_dic = PrepareData.merge_pos_gap(opt, opt.pos_gap)
        PrepareData.save_merge_pos(opt, merge_lst)
        # data = pd.DataFrame(data=merge_lst, columns=['name', 'chr', 'start', 'end', 'length'])
        # data.to_csv(opt.path + '/merge_pos.csv', index=False)
    else:
        merge_dic = {}
        data = pd.read_csv(opt.path + '/merge_pos.csv')
        gl.report_dic['mergeTargetLen'] = data.length.sum()

    logging.info('step2: HomoDetection...')
    # if opt.cov == 2:
    design_pos, homo_0, homo_1 = PrepareData.homo_detect(opt, ref_dic, merge_dic)
    design_pos.to_csv(opt.path + '/new_pos_tag.txt', sep='\t', index=False)


    logging.info('step3: Start short sequence Alignment...')
    if not os.path.exists(opt.path + '/Hits_result.ln.gz.idx'):
        Aln = PrepareData.aln_sequence(opt)

    logging.info('step4: load the hits index from alignment...')
    MPD.load_hit_index(opt.path)

    logging.info('step5: Start to design primers...')
    result_df, failed_df = MPD.primer_design(opt, homo_0, ref_dic)

    logging.info(f"\n{'=' * 30}\nstep6: Start to evaluate ...")
    hR.write_html_report(opt, result_df, failed_df, homo_1)
    dimer = PE.dimers_detect(opt.path, opt.output + '.xlsx')
    PA.create_sam_file(opt.ref_path, opt.path, result_df)
    head = PA.load_head(f'{opt.path}/sam_head')
    zip_result(opt.scr_dir, opt.path, opt.order_num)
    logging.info(f"All CAMPs designs have completed.")
    with open(opt.path + '/DONE', 'w') as f:
        f.write(f'{time.asctime(time.localtime(time.time()))}')




if __name__ == '__main__':
    main()


