#!/usr/bin/python3
__author__ = 'sunxin'

from lib_sam_line import *
from lib_sam_file import *

import os, sys
import argparse


def sam2con(in_fh, ancient = False, filt_db = False, niter = 3, out_ambi = False, filter_freq = 0.3) :
    # generate consensus
    s_fh = open(in_fh, 'r')

    # read in sam file
    if filt_db :
        sfh = sam_file(s_fh, filter_file = filt_db)
    else :
        sfh = sam_file(s_fh)

    sfh.iter_con(iter_time = niter, filter_ancient= ancient, low_freq = filter_freq)
    seq_con = sfh.conseq()

    print(seq_con[0])
    print(seq_con[1])

    if out_ambi :
        sfh.out_ambi()

    s_fh.close()


def sam_snp_stat(in_fh , read_length = 150) :
    # generate samfile snp status

    s_fh = open(in_fh, 'r')

    sfh = sam_file(s_fh)
    snpstat = sfh.snp_stat(read_length = read_length)

    pr_order = ['AT', 'AC', 'AG','AN', 'TA', 'TC', 'TG', 'TN', 'CA', 'CG', 'CT', 'CN', 'GA', 'GT', 'GC', 'GN' ]

    print(snpstat[2])
    # output result
    print('\t'.join(['snp_L'] + pr_order))
    for i in range(0,150) :
        outl = [str(i)]
        for j in range(0, 16) :
            outl.append(str(snpstat[0][str(i)][pr_order[j]]))
        print('\t'.join(outl))

    print('\t'.join(['snp_R'] + pr_order))
    for i in range(0,150) :
        outl = [str(i)]
        for j in range(0, 16) :
            outl.append(str(snpstat[1][str(i)][pr_order[j]]))
        print('\t'.join(outl))

    s_fh.close()






def arg():
    # argument managemnt

    parse = argparse.ArgumentParser(prog='tools_sunxin_sam')

    sub = parse.add_subparsers(help='see subcommand')

    sub1 = sub.add_parser('consensus', help='build consensus sequence')
    sub1.add_argument('-i', '--input', help='input sam file, without S', required=True,
                      type=str)
    sub1.add_argument('-an','--ancient', help='filter for ancient marker : default=False',
                      action='store_true')
    sub1.add_argument('-f', '--fdb', help='filter database file, default=False',
                      default='False', type=str)
    sub1.add_argument('-n', '--iter', help='iteration times, default=1', type=int, default=1)
    sub1.add_argument('-x', '--out_ambi', help='output ambigous base frequency, default=False',
                      action='store_true')
    sub1.add_argument('-c', '--cutoff', help='filter low frequency cutoff, default=0.3',
                      default=0.3, type=float)

    sub2 = sub.add_parser('snp_stat', help = 'generate snp status for sam file')
    sub2.add_argument('-i', '--input', help='input sam file, without S', required=True,
                      type=str)
    sub2.add_argument('-l', '--length', help='read length, default=150', type=int,
                      default=150)

    args = parse.parse_args()


    if 'cutoff' in vars(args).keys() :
        sam2con(
            str(args.input),
            ancient=args.ancient,
            filt_db= args.fdb if args.fdb != 'False' else False,
            niter=args.iter,
            out_ambi=args.out_ambi,
            filter_freq=args.cutoff
        )
    if 'length' in vars(args).keys() :
        sam_snp_stat(
            str(args.input),
            read_length=args.length
        )


if __name__ == '__main__' :

    arg()




