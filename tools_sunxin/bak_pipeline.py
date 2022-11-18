#!/usr/bin/python3
__author__ = 'sunxin'

import multiprocessing
import os, argparse
from lib_pipeline import *



def run_indi(run_arg) :
    sample_name = str(run_arg[0])
    s_name = str(run_arg[1])
    wdir = str(run_arg[2])
    threads = run_arg[3]
    do_all = run_arg[4]
    cuta = run_arg[5]
    cutqual = run_arg[6]
    cutmin = run_arg[7]
    fastqc = run_arg[8]
    mt_map = run_arg[9]
    wg_map = run_arg[10]
    cuta_ss = run_arg[11]

    file_obj = file(name=sample_name, short_name=s_name, work_dir=wdir)

    if do_all :
        file_obj.cutada(qual=cutqual, min=cutmin)
        file_obj.fqc()
        file_obj.mt_map(t=threads)
        file_obj.wg_map(t=threads)

    if cuta :
        file_obj.cutada(qual=cutqual, min=cutmin)
    if cuta_ss :
        file_obj.cutada_ss(qual=cutqual, min=cutmin)
    if fastqc :
        file_obj.fqc()
    if mt_map :
        file_obj.mt_map(t=threads)
    if wg_map :
        file_obj.wg_map(t=threads)



def run(indi=False, s_name=False, indi_file=False,
        np = 10, wdir=False, threads=False,
        do_all=False,
        cuta=False, cuta_ss=False, cutqual=False, cutmin=False,
        fastqc=False, mt_map=False, wg_map=False) :
    sample_list = []
    s_name_list = []
    if indi_file :
        sample_file = open(indi_file, 'r')
        while 1 :
            sample_file_line = sample_file.readline().strip().split(' ')
            if len(sample_file_line) == 1 :
                break
            sample_list.append(sample_file_line[0])
            s_name_list.append(sample_file_line[1])

    if indi :
        sample_list = indi.split(',')
        s_name_list = s_name.split(',')

    if len(sample_list) == 0 :
        return 0

    list_len = len(s_name_list)

    if len(s_name_list) <= 10 :
        mp_pool = multiprocessing.Pool(len(s_name_list))
    else :
        mp_pool = multiprocessing.Pool(int(np))

    run_arg = zip(sample_list,
                  s_name_list,
                  [wdir] * list_len,
                  [threads] * list_len,
                  [do_all] * list_len,
                  [cuta] * list_len,
                  [cutqual] * list_len,
                  [cutmin] * list_len,
                  [fastqc] * list_len,
                  [mt_map] * list_len,
                  [wg_map] * list_len,
                  [cuta_ss] * list_len)

    mp_pool.map(run_indi, run_arg)



def arg() :
    '''
    argument management
    '''

    parse = argparse.ArgumentParser(prog='pipeline')

    parse.add_argument('-i', '--input', help='input indivdual, xxx,xxx', type=str)
    parse.add_argument('-s', '--s_name', help='short name, xxx,xxx', type=str)
    parse.add_argument('-f', '--file', help='input file list, default=False', type=str, default=False)
    parse.add_argument('--wdir', help='work directory (base), default=. ', type=str, default='.')
    parse.add_argument('-t', '--threads', help='threads, default= 10', type=int, default=10)

    parse.add_argument('--all', help='conduct all analysis, default=False',
                       action='store_true')
    parse.add_argument('--cutadapt', help='cut adaptor and trim sequence, default=False',
                       action='store_true')
    parse.add_argument('--cutadapt_ss', help='cut adaptor for ss library, default=False',
                       action='store_true')
    parse.add_argument('--cutqual', help='cut adaptor min qual, default=30', type=int, default=30)
    parse.add_argument('--cutmin', help='cut adaptor min length, default=40', type=int, default=40)
    parse.add_argument('--fastqc', help='run fastqc on trimmed files, default=False',
                       action='store_true')
    parse.add_argument('--mt_map', help='map to mitochondria sequence, default=False',
                       action='store_true')
    parse.add_argument('--wg_map', help='map to whole genome, default=False',
                       action='store_true')

    args =parse.parse_args()
    run(indi=args.input if args.file == False else False,
        s_name=args.s_name if args.file == False else False,
        indi_file= args.file if args.file != False else False,
        do_all=args.all if args.all != False else False,
        cuta=args.cutadapt if args.cutadapt != False else False,
        cuta_ss=args.cutadapt_ss if args.cutadapt_ss != False else False,
        cutqual=args.cutqual,
        cutmin=args.cutmin,
        fastqc=args.fastqc if args.fastqc != False else False,
        mt_map=args.mt_map if args.mt_map != False else False,
        wg_map=args.wg_map if args.wg_map != False else False,
        threads=args.threads,
        wdir=args.wdir

        )



if __name__ == '__main__' :
    arg()


