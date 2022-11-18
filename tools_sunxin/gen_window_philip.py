#!/usr/bin/python3
__author__ = 'sunxin'


from lib_split_file import *
import vcf
import sys, os
import multiprocessing
from lib_cmd import *
from pyfasta import Fasta


def read_bed(bed_file, wdir) :
    '''
    read BED file as a list
    '''

    os.chdir(wdir)

    b_fh = open(bed_file, 'r')

    bed_list = []
    while 1 :
        l_bfh = b_fh.readline().strip().split('\t')

        if len(l_bfh) == 1 :
            break

        # BED format [0,100)
        # position follow BED format
        chr = str(l_bfh[0])
        p_s = int(l_bfh[1])
        p_e = int(l_bfh[2])

        bed_list.append([chr, p_s, p_e]) # as python position

    b_fh.close()

    return bed_list


def gen_philip(sample_ID_list, bed_list, fa_list, wdir) :
    '''
    generate philip file for each bed position
    '''

    os.chdir(wdir)
    ofh = open('_'.join([str(i) for i in bed_list]) + ".philip", 'w')

    # print philip head
    print('\t'.join([str(len(sample_ID_list)),
                     str(int(bed_list[2]) - int(bed_list[1]))]),
          file=ofh)
    # print philip body
    sample_ID_list.sort()
    for i in sample_ID_list :
        print('\t'.join([i,
                         fa_list[i][bed_list[0]][bed_list[1]:bed_list[2]]]),
              file=ofh)

    ofh.close()


def run_indi(arg_indi) :

    ''' individual run '''

    wdir = arg_indi['wdir']
    # read bed file
    bed_file = arg_indi['child_bed']
    bed_list = read_bed(bed_file, wdir)

    # read fasta file
    sample_list = arg_indi['sample_list']
    sample_ID_list = arg_indi['sample_ID_list']
    sample_ID_list.sort()

    fa_list = {}
    for i in sample_ID_list :
        fa_list[i] = Fasta(sample_list[i])

    # generate philip for each pos
    for i in bed_list :
        gen_philip(sample_ID_list, i, fa_list, wdir)


def main(np, sample_list_file, bed_file, out_dir) :

    # read current dir
    start_dir = os.getcwd()

    # read sample list
    sample_list = {}
    sample_ID_list = []
    s_fh = open(sample_list_file, 'r')

    while 1 :
        l_sfh = s_fh.readline().strip().split(' ')

        if len(l_sfh) != 2 :
            break

        sample_list[str(l_sfh[0])] = str(l_sfh[1])
        sample_ID_list.append(l_sfh[0])
    s_fh.close()

    # split bed
    child_bed = split_file(file=bed_file,
                           np=np,
                           wdir=start_dir)
    child_bed_list = child_bed.split()

    # child file mp run
    arg_list = []
    arg_indi  = {
        'child_bed' : '',
        'wdir' : '',
        'sample_list' : '',
        'sample_ID_list' : '',

    }

    arg_indi['wdir'] = start_dir + out_dir
    arg_indi['sample_list'] = sample_list
    arg_indi['sample_ID_list'] = sample_ID_list

    for i in child_bed_list :
        arg_indi['child_bed'] = start_dir + '/tmp_' + bed_file + '/' + i
        arg_app = arg_indi.copy()
        arg_list.append(arg_app)

    mp_pool = multiprocessing.Pool(int(np))
    mp_pool.map(run_indi, arg_list)

    # reduce result
    child_bed.destory()


if __name__ == '__main__' :
    print('Window generate Philip file from WG FASTA. \n'
          'Parallel running supported. \n'
          'sample_list contain : sample_ID \\s absolute_fa_file_position. \n'
          '\n'
          'Usage: \n'
          'gen_window_philip.py NP sample_list window.Bed out_dir\n')

    if len(sys.argv) == 1 :
        print('PLEASE INPUT')
    else :
        main(sys.argv[1],       # NP
             sys.argv[2],       # sample_list
             sys.argv[3],       # window.bed
             sys.argv[4],       # out_dir
             )