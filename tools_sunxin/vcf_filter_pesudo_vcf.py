#!/usr/bin/python3

__author__ = 'sunxin'

'''

Take filtered pesudo vcf.
Generate vcf with all homozygous sites.

Strategy:
    Random take 1 read from each loci.

'''

from lib_split_file import *
import vcf
import os
import multiprocessing
from lib_cmd import *
import argparse
import random

def line_filter(vl, sample_list) :


    for i in sample_list :
        if vl.genotype(i)['GT'] == '0/1' :
            ad = vl.genotype(i)['AD']
            if ad == [0, 0] :
                return 0
            gt_list = int(ad[0]) * ['0'] + int(ad[1]) * ['1']
            random.seed()
            random.shuffle(gt_list)
            gt_homo_random = random.sample(gt_list, 1)[0]
            gt_homo = gt_homo_random + '/' + gt_homo_random
            vl.genotype(i).data = vl.genotype(i).data._replace(GT=gt_homo)

    return 1


def v_f_indi(arg_indi) :
    os.chdir(arg_indi['wdir'])
    fh = vcf.Reader(open(arg_indi['child_vcf'], 'r'))

    ovcf = open('tmp_' + arg_indi['child_vcf'] + '.vcf', 'w')
    tmp_ovcf = vcf.Writer(ovcf, fh)
    sample_list = fh.samples

    for lh in fh :
        olf = line_filter(lh, sample_list)
        if olf == 0 :
            continue
        else :
            tmp_ovcf.write_record(lh)
    ovcf.close()

    get_body = cmd_out(cmd_line='grep -v ^# CMDARG1',
                           in_p=['tmp_' + arg_indi['child_vcf'] + '.vcf'],
                           work_dir=arg_indi['wdir'],
                           shell=True,
                           wait=True,
                           outfile='./tmp_body_' + arg_indi['child_vcf'] + '.vcf')
    get_body.run()


def v_f(args) :
    '''
    vcf filter
    '''

    # decode args
    in_vcf = args.vcf
    np = args.np

    # in vcf split
    if args.continue_run :
        start_dir = os.getcwd()
        child_vcf_list = []
        for i in range(0,np) :
            if i <= 9 :
                child_vcf_list.append("tmp_0" + str(i) + "_" + in_vcf)
            else :
                child_vcf_list.append("tmp_" + str(i) + "_" + in_vcf)

    else :
        start_dir = os.getcwd()
        in_vcf_fh = split_file(file=in_vcf,
                               np=np,
                               wdir=start_dir)
        child_vcf_list = in_vcf_fh.split_head()

    # child file mp run
    arg_list = []
    arg_indi  = {
        'child_vcf'     : '',
        'wdir'          : start_dir + '/tmp_' + in_vcf
    }

    arg_indi['wdir'] = start_dir + '/tmp_' + in_vcf
    for i in child_vcf_list :
        arg_indi['child_vcf'] = i
        arg_app = arg_indi.copy()
        arg_list.append(arg_app)
    mp_pool = multiprocessing.Pool(int(np))
    mp_pool.map(v_f_indi, arg_list)

    # reduce result
    child_out_vcf_file_list = []

    for i in child_vcf_list :
        child_out_vcf_file_list.append('tmp_body_' + i + '.vcf')



    reduce_child_vcf = cmd_out(cmd_line='cat CMDARG1',
                           in_p=[' '.join(['tmp_head_' + in_vcf] + child_out_vcf_file_list)],
                           work_dir=start_dir + '/tmp_' + in_vcf,
                           outfile=start_dir + '/' + args.name + in_vcf,
                           shell=True,
                           wait=True)
    reduce_child_vcf.run()

    if not args.no_rmdir :
        in_vcf_fh.destory()


def arg() :
    '''
    argument
    '''

    parse = argparse.ArgumentParser(prog = 'vcf_filter_pesudo_vcf')

    parse.add_argument('-p', '--np',
                       help=' maximum number of pool workers, default= 10',
                       type=int, default=10)
    parse.add_argument('-i', '--vcf',
                       help='input vcf file', required=True, type=str)
    parse.add_argument('--name', help="name affix to input vcf filename, default=info_", type=str,
                       default='info_')
    parse.add_argument('--continue_run', help="continue running without new split file", action='store_true',
                       default=False)
    parse.add_argument('--no_rmdir', help="running without rmdir", action='store_true',
                       default=False)

    args = parse.parse_args()

    v_f(args)


if __name__ == '__main__' :
    arg()
