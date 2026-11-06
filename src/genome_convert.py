#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from threading import Thread
import subprocess
from os.path import getsize

name = sys.argv[1]
#name = 'hg19'

with open(name+'.fa') as r:
    genome = r.read().upper()
glen = len(genome)

def CT_tran(genome,glen):
    genome_list = list(genome)
    rename = False
    for i in range(glen-2):
        if genome[i] == '>':
            rename = True
        elif rename:
            if genome[i] != '\n':
                continue
            else:
                genome_list[i] = '_CT\n'
                rename = False
        elif genome[i] == 'C':
            if genome[i+1] == 'G' or genome[i+1:i+3] == '\nG':
                continue
            genome_list[i] = 'T'
        # elif genome[i:i+2] == 'TG' or genome[i+1:i+3] == 'T\nG':
            # genome_list[i] = 'C'
        else:
            pass
    with open(name+'_CT.fa','w')as w:
        w.write(''.join(genome_list))

def GA_tran(genome,glen):
    genome_list = list(genome)
    rename = False
    for i in range(glen):
        if genome[i] == '>':
            rename = True
        elif rename:
            if genome[i] != '\n':
                continue
            else:
                genome_list[i] = '_GA\n'
                rename = False
        elif genome[i] == 'G':
            if genome[i-1] == 'C' or genome[i-2:i] == 'C\n':
                continue
            genome_list[i] = 'A'
        # elif genome[i-1:i+1] == 'CA' or genome[i-2:i+1] == 'C\nA':
            # genome_list[i] = 'G'
        else:
            pass
    with open(name+'_GA.fa','w')as w:
        w.write(''.join(genome_list))

t1= Thread(target=CT_tran,args=(genome,glen))
t1.start()
t2= Thread(target=GA_tran,args=(genome,glen))
t2.start()

t1.join()
t2.join()

subprocess.run('cat {n}_CT.fa {n}_GA.fa > {n}.tran.fa'.format(n=name), shell=True)
