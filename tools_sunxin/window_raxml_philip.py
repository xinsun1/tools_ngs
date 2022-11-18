#!/usr/bin/python3

__author__ = 'sunxin'

import sys, os
import multiprocessing
from lib_split_file import *
from lib_cmd import *


def run_indi(arg_indi) :
    '''
    read BED file as a list
    '''

    bed_file = arg_indi['child_bed']
    wdir = arg_indi['wdir']

    os.chdir(wdir)

    b_fh = open(bed_file, 'r')

    while 1 :
        l_bfh = b_fh.readline().stri\
            ().split('\t')

        if len(l_bfh) == 1 :
            break

        philip_name = '_'.join(l_bfh)
        run_raxml = cmd(cmd_line='raxmlHPC-HYBRID-AVX -f a -x 1234 -p 1234 ' \
         '-T 2 -# 100 -m GTRGAMMA -s CMDARG1.philip -n CMDARG1',
                        in_p=[philip_name],
                        work_dir=wdir,
                        shell=False,
                        wait=True,
                        )
        run_raxml.run()
    b_fh.close()


def run(np, bed_file, wdir):

    # read current dir
    start_dir = os.getcwd()

    # split bed file
    child_bed =  split_file(file=bed_file,
                            np=np,
                            wdir=start_dir)
    child_bed_list = child_bed.split()

    # child file mp run
    arg_list = []
    arg_indi = {
        'child_bed' : '',
        'wdir' : wdir,
    }

    for i in child_bed_list :
        arg_indi['child_bed'] = start_dir + '/tmp_' + bed_file + '/' + i
        arg_app = arg_indi.copy()
        arg_list.append(arg_app)

    mp_pool = multiprocessing.Pool(int(np))
    mp_pool.map(run_indi, arg_list)

    # reduce result
    child_bed.destory()


if __name__ == '__main__' :
    print('Run Raxml for windowed philip file\n'
          'Each individual run will use 2 threads\n'
          'Usage: \n'
          'window_raxml_philip.py NP window.Bed philip_file_dir\n')

    if len(sys.argv) == 1 :
        print('PLEASE INPUT')
    else :
        run(sys.argv[1],       # NP
            sys.argv[2],       # window.bed
            sys.argv[3],       # philip directory
            )
