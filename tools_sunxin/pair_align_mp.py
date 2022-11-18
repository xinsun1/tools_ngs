#!/usr/bin/python3
__author__ = 'sunxin'

from pair_align import *
from lib_split_file import *
import sys
import subprocess
import multiprocessing
import time

def split(input, np, chunk) :
    '''
    split file

    input
        input file, INPUT_1.fq, INPUT_2.fq
    np
        number of processor
    chunk
        file chunksize, for
    return
        [tmp_xxx_INPUT]
    '''

    f1 = split_file(input + '_1.fq', np, chunksize=chunk)
    f2 = split_file(input + '_2.fq', np, chunksize=chunk)
    o1 = f1.split()
    f2.split()
    out = [ i.split('_1.fq')[0] for i in o1 ]
    return out

def merge_file(out_list) :
    out_name = '_'.join(out_list[0].split('_')[2:])

    out_aligned = [ i + '_aligned.fa' for i in out_list]
    out_un1 = [ i + '_unaligned_1.fq' for i in out_list]
    out_un2 = [ i + '_unaligned_2.fq' for i in out_list]
    tmp_fq1 = [ i + '_1.fq' for i in out_list]
    tmp_fq2 = [ i + '_2.fq' for i in out_list]


    with open(out_name + '_aligned.fa', 'w') as ofh1 :
        subprocess.call(['cat'] + out_aligned, stdout=ofh1)
    ofh1.close()
    with open(out_name + '_unaligned_1.fq', 'w') as ofh2 :
        subprocess.call(['cat'] + out_un1, stdout=ofh2)
    ofh2.close()
    with open(out_name + '_unaligned_2.fq', 'w') as ofh3 :
        subprocess.call(['cat'] + out_un2, stdout=ofh3)
    ofh3.close()

    for i in out_aligned :
        subprocess.Popen(['rm', i])
    for i in out_un1 :
        subprocess.Popen(['rm', i])
    for i in out_un2 :
        subprocess.Popen(['rm', i])
    for i in tmp_fq1 :
        subprocess.Popen(['rm', i])
    for i in tmp_fq2 :
        subprocess.Popen(['rm', i])


if __name__ == '__main__' :
    inputfile = sys.argv[1]
    np = sys.argv[2]
    chunksize = sys.argv[3]

    start_time = time.time()
    file_list = split(input=inputfile, np=np, chunk=chunksize)
    split_time = time.time()
    print('Split file time :' + str(split_time - start_time))

    mp_pool = multiprocessing.Pool(int(np))
    align_list = mp_pool.map(align, file_list)
    merge_file(align_list)
    finish_time = time.time()
    print('Finished.')
    print('Total time :' + str(finish_time - start_time))





