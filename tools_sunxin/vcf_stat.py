#!/usr/bin/python3

__author__ = 'sunxin'

'''
Generate summary statistics for vcf files
Pre-process for using R.

Multi-threads supported.
'''


from lib_cmd import *
import argparse, os
import multiprocessing
import numpy
from lib_split_file import *
import vcf


def v_stat(arg_indi) :
    #### get args ####
    vcf_file = arg_indi['child_vcf']
    wdir     = arg_indi['wdir']

    #### mut type ####
    mut = {
        'AT' : 0,
        'AC' : 1,
        'AG' : 2,
        'TA' : 3,
        'TC' : 4,
        'TG' : 5,
        'CA' : 6,
        'CT' : 7,
        'CG' : 8,
        'GA' : 9,
        'GT' : 10,
        'GC' : 11
     }

    os.chdir(wdir)
    fh = vcf.Reader(open(vcf_file), 'r')

    sample_list = fh.samples
    sample_list.sort()
    nsample = len(sample_list)

    ar_indi = numpy.zeros((4, 12, nsample), dtype=numpy.dtype(int))      # [gt_index, mut_index, ID]

    for lh in fh:
        lh_mut = str(lh.REF) + str(lh.ALT[0])

        for i in range(nsample) :
            gt_i = lh.genotype(sample_list[i]).gt_type
            if gt_i == None :
                gt_index = 3
            else :
                gt_index = int(gt_i)
            ar_indi[gt_index, mut[lh_mut], i] += 1

    return [ar_indi, sample_list]


def run(args) :
    #### get args ####
    np     = args.np
    in_vcf = args.vcf
    o_name = args.out

    #### in vcf split ####
    start_dir = os.getcwd()
    in_vcf_fh = split_file(file = in_vcf,
                           np   = np,
                           wdir = start_dir)
    child_vcf_list = in_vcf_fh.split_head()

    #### child file mp run ####
    arg_list = []
    arg_indi = {
        'child_vcf' : '',
        'wdir'      : '',
    }
    arg_indi['wdir'] = start_dir + '/tmp_' + in_vcf
    for i in child_vcf_list :
        arg_indi['child_vcf'] = i
        arg_app = arg_indi.copy()
        arg_list.append(arg_app)
    mp_pool = multiprocessing.Pool(int(np))
    re_v2p = mp_pool.map(v_stat, arg_list)

    #### reduce result ####
    os.chdir(start_dir)
    o_fh = open(o_name, 'w')

    sample_list = re_v2p[0][1]
    ar_total    = re_v2p[0][0]

    for i in range(1, len(re_v2p)) :
        ar_total += re_v2p[i][0]

    mut = ['AT', 'AC', 'AG', 'TA', 'TC', 'TG',
           'CA', 'CT', 'CG', 'GA', 'GT', 'GC']

    dim_ar = ar_total.shape
    print('\t'.join(['Count', 'GT', 'Mut_type', 'Sample_ID']), file=o_fh)
    for index_gt in range(dim_ar[0]) :
        for index_mut in range(dim_ar[1]) :
            for index_sample in range(dim_ar[2]) :
                ol = [str(ar_total[index_gt, index_mut, index_sample]),
                      str(index_gt),
                      mut[index_mut],
                      sample_list[index_sample]
                      ]
                print('\t'.join(ol), file=o_fh)

    o_fh.close()
    in_vcf_fh.destory()


def arg() :
    '''
    argument
    '''

    parse = argparse.ArgumentParser(prog = 'vcf_stat')

    parse.add_argument('-p', '--np',  help=' maximum number of pool workers, default= 10',
                       type=int, default=10)
    parse.add_argument('-i', '--vcf', help='input vcf file', type=str, default=False)
    parse.add_argument('-o', '--out', help='out file name',  type=str, required=True)

    args = parse.parse_args()

    run(args)


if __name__ == '__main__' :
    arg()
