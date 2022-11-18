#!/usr/bin/python3

__author__ = 'sunxin'

'''
Remove filtered snp for no outgroup snp dataset.

Strategy:

1) read file, Out group already removed
2) calculate allele freq
3) if freq 0 or 1 == 0, remove snp
'''

from lib_split_file import *
import vcf
import os
import multiprocessing
from lib_cmd import *
import argparse


def line_filter(vl, sample_list) :

    # not consider OUT_GROUP snps
    # count allele freq, only support for 2 allele
    vl_c = {"." : 0,
            "0" : 0,
            "1" : 0}

    for i in sample_list :
        gt_list = vl.genotype(i)['GT'].split('/')
        vl_c[str(gt_list[0])] += 1
        vl_c[str(gt_list[1])] += 1

    # apply filter for snps not exist in these samples
    # check if the SNP exist in OUT removed samples
    if vl_c["0"] == 0 or vl_c["1"] == 0 :
        return 0
    else :
        return 1

   # 0 : is a out group snp
   # 1 : is a tiger SNP

def v_f_indi(arg_indi) :
    os.chdir(arg_indi['wdir'])
    fh = vcf.Reader(open(arg_indi['child_vcf']), 'r')

    sample_list_all = fh.samples

    ovcf = open('tmp_' + arg_indi['child_vcf'] + '.vcf', 'w')
    tmp_ovcf = vcf.Writer(ovcf, fh)

    for lh in fh :
        f_out = line_filter(lh, sample_list_all)
        if f_out == 1 :
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
    out_vcf = args.out

    # in vcf split
    start_dir = os.getcwd()
    in_vcf_fh = split_file(file=in_vcf,
                           np=np,
                           wdir=start_dir)
    child_vcf_list = in_vcf_fh.split_head()

    # child file mp run
    arg_list = []
    arg_indi  = {
        'child_vcf' : '',
        'wdir' : '',
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
                           outfile=start_dir + '/' + out_vcf,
                           shell=True,
                           wait=True)
    reduce_child_vcf.run()

    in_vcf_fh.destory()


def arg() :
    '''
    argument
    '''

    parse = argparse.ArgumentParser(prog = 'vcf_filter')

    parse.add_argument('-p', '--np',
                       help=' maximum number of pool workers, default= 10',
                       type=int, default=10)
    parse.add_argument('-i', '--vcf',
                       help='input vcf file', required=True, type=str)
    parse.add_argument('-o', '--out',
                       help="output file name", required=True, type=str)

    args = parse.parse_args()

    v_f(args)


if __name__ == '__main__' :
    arg()
