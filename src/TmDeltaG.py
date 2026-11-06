#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Time  : 2019/4/29 17:37
# @Author: Chentian
# @File  : TmDeltaG.py

import math
import Seq
import sys
import ThermodynamicsParameters as TP


class Cal():
    def __init__(self, qseq, sseq, mono_conc= 50, diva_conc=1.5, oligo_conc=50, dntp_conc=0.25):
        deltaH, deltaS = calDeltaHS(qseq, sseq)

        self.Tm = calTm(qseq, sseq, mono_conc=mono_conc, diva_conc=diva_conc, oligo_conc=oligo_conc,
                        dntp_conc=dntp_conc, deltaH=deltaH, deltaS=deltaS)
        self.DeltaG = calDeltaG(qseq, sseq, mono_conc=mono_conc, diva_conc=diva_conc, dntp_conc=dntp_conc,
                                deltaS=deltaS, deltaH=deltaH)


def calDeltaHS(qseq, sseq):
    init_begin = 'init%s%s' % (qseq[0], sseq[0])
    init_end = 'init%s%s' % (qseq[-1], sseq[-1])

    if init_begin in TP.dH_full and init_end in TP.dH_full:
        deltaH = TP.dH_full[init_begin] + TP.dH_full[init_end]
        deltaS = TP.dS_full[init_begin] + TP.dS_full[init_end]
    else:
        deltaH = 0
        deltaS = 0

    for i in range(len(qseq) - 1):
        dinuc = '%s%s' % (qseq[i:(i+2)], sseq[i:(i+2)])

        if dinuc in TP.dS_full and dinuc in TP.dS_full:
            deltaH = deltaH + TP.dH_full[dinuc]
            deltaS = deltaS + TP.dS_full[dinuc]

    return round(deltaH,2), round(deltaS,2)


def calDeltaG(qseq, sseq, mono_conc=50, diva_conc=1.5, dntp_conc=0.25, deltaH=None, deltaS=None):
    """ Calculate the free Gibbs energy"""
    mono_conc = float(mono_conc)
    diva_conc = float(diva_conc)
    dntp_conc = float(dntp_conc)

    if not (deltaH and deltaS):
        deltaH, deltaS = calDeltaHS(qseq, sseq)
    tao = 273.15 + 37

    mono_conc = mono_conc + divalent2monovalen(diva_conc, dntp_conc)
    mono_conc = mono_conc / 1000

    deltaS_adjust = deltaS + 0.368 * (len(sseq) - 1) * math.log(mono_conc, math.e)

    deltaG = (deltaH * 1000 - tao * deltaS_adjust) / 1000

    return round(deltaG,2)


def divalent2monovalen(divalent, dntp):
    """ Divalent to monovalent"""
    if divalent == 0:
        dntp = 0
    if divalent < 0 or dntp <0:
        print >> sys.stderr, 'Error: divalent2monovalen'
        exit()
    if divalent < dntp:
        # According to the theory, melting temperature does not depends on divalent cations
        divalent = dntp

    return 120 * (math.sqrt(divalent - dntp))


def calTm(qseq, sseq, mono_conc=50, diva_conc=1.5, oligo_conc=50, dntp_conc=0.25, deltaH=None, deltaS=None):
    """ Calculate Tm value of amplicon """

    mono_conc = float(mono_conc)
    diva_conc = float(diva_conc)
    oligo_conc = float(oligo_conc)
    dntp_conc = float(dntp_conc)

    if not (deltaH and deltaS):
        deltaH, deltaS = calDeltaHS(qseq, sseq)

    deltaH = deltaH * 1000
    oligo_conc = oligo_conc / 1e9

    mono_conc = mono_conc + divalent2monovalen(diva_conc, dntp_conc)
    mono_conc = mono_conc / 1000

    deltaS = deltaS + 0.368 * (len(sseq) - 1) * math.log(mono_conc, math.e)

    Tm = deltaH / (deltaS + 1.987 * math.log(oligo_conc / 4, math.e)) - 273.15

    return round(Tm,2)

def main():
    qseq = 'ATAACTCACGCCTATAATTCCAAC'
    # qseq = 'TGCGGAGCGCGCACGCGGGAACCCGCGCTGGAGGCGGGCGAGGGCCGAGGGGCAGCTAGGGAGCGCGGCTTGAGGAGGGCGGGGCCGCCC'
    sseq = Seq.complement(qseq)
    # sseq = Seq.reverse(sseq)
    print(sseq.upper())
    print(qseq,sseq)
    print(Seq.rev_com(qseq))

    mono = 50
    diva = 3  # our mg2+ is 3mM, default=1.5
    oligo = 50
    dntp = 0.8  #our dntp is 0.2*4=0.8, default=0.25
    seq = Cal(qseq.upper(), sseq.upper(), mono_conc=mono, oligo_conc=oligo, diva_conc=diva, dntp_conc=dntp)

    print('Tm: ', seq.Tm)
    print('DeltaG: ', seq.DeltaG)


if __name__ == '__main__':
    main()