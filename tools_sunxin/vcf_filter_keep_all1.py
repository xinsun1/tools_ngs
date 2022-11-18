#!/usr/bin/python3

__author__ = 'sunxin'

'''
Lasted update : keep all 1 site for OUT snp

Raw vcf file filter generated from GATK.

Take input from bcftools filtered vcf with only biallelic SNP allowed.
bcftooll remove indel and allow biallelic site.

Filter strategy :
1) Read sample population file for MODERN and AMO pop sample.
2) check if --out to allow for SNPs only from OUT pop.
   If allowed, check MIN_DP_OUT of OUT snp.
3) For Modern, check 



To do :
modify out vcf file name with OUT information added.
'''

from lib_split_file import *
import vcf
import os
import multiprocessing
from lib_cmd import *
import argparse


def line_filter(vl, sample_list, sample_dict, out_list=False) :

    MIN_DP_MODERN = 7
    MIN_DP_AMO = 3
    MIN_DP_OUT = 20
    MIN_DP_PTV = 4
    

    out_list = out_list
    # if allow OUT_GROUP snps
    out_snp = 0
    if out_list :

        for i in out_list :
            gt = vl.genotype(i)['GT']

            # if SNP in outgroup, check DP
            if gt == "0/1" or gt == "1/1" :
                dp = vl.genotype(i)['DP']
                if dp == None :
                    continue

                elif int(dp) >= MIN_DP_OUT :
                    out_snp = 1
                    break

    # get mutation type
    mut_type = str(vl.REF) + str(vl.ALT[0])

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
    if vl_c["1"] == 0 :
        return [2, out_snp]
    if vl_c["0"] == 0 :
        return [3, out_snp]


    # apply filter for snps do not sequenced in 80% of sample_list
    # check for NA rate in these samples, frequency allowed 50%
    if float(vl_c["0"] + vl_c["1"]) / 2 / len(sample_list) <= 0.5 :
        return [0, out_snp]

    # apply filter for DP
    amo_count = 0
    ptv_count = 0
    dp_pass_ptv = 0
    dp_pass_amo = 0
    ptv06_pass = 0
    for i in sample_list :
        gt = vl.genotype(i)['GT']
        dp = vl.genotype(i)['DP']

        if gt != "0/0" and gt != "./." :
            if dp == None :
                continue
            if sample_dict[i] == "MODERN" :
                if int(dp) >= MIN_DP_MODERN :
                    return [1, out_snp]
            # special filter for amo samples :
            # if SNP was only found in amo samples,
            # this snp should exist in at least 3 samples
            if sample_dict[i] == "AMO" :
                amo_count += 1
                if int(dp) >= MIN_DP_AMO :
                    dp_pass_amo = 1
                    if amo_count >= 3 and dp_pass_amo == 1 :
                        return [1, out_snp]


            # specific for PTV samples
            # if heterozygous SNP was only found in PTV samples,
            # this SNP should exist in at least 2 samples or
            # GQ >60


            if sample_dict[i] == "PTV" :
                ptv_count += 1
                if int(dp) >= MIN_DP_PTV :
                    dp_pass_ptv = 1
                    # exist in PTV06 and another one
                    if i == 'PTV06' :
                        ptv06_pass = 1
                    if ptv_count >= 2 and dp_pass_ptv == 1 and ptv06_pass == 1:
                        return [1, out_snp]


    return [0, out_snp]    # 0 : is a tiger SNP but failed
                           # 1 : is a tiger SNP pass
                           # 2 : not a tiger SNP
                           # 3 : all 1 site


def v_f_indi(arg_indi) :
    os.chdir(arg_indi['wdir'])
    fh = vcf.Reader(open(arg_indi['child_vcf']), 'r')

    sample_dict = arg_indi['filter_list']
    sample_list = arg_indi['filter_list'].keys()

    sample_list_all = fh.samples
    out_list = arg_indi['out_list']

    ovcf = open('tmp_' + arg_indi['child_vcf'] + '.vcf', 'w')
    tmp_ovcf = vcf.Writer(ovcf, fh)

    for lh in fh :
        f_out = line_filter(lh, sample_list, sample_dict, out_list)
        if f_out[0] == 1 :
            tmp_ovcf.write_record(lh)
        elif f_out[0] == 3 :
            # write a all 1 snp
            tmp_ovcf.write_record(lh)
        elif f_out[0] == 0 and f_out[1] == 1 and arg_indi['mod_gt'] == True:
            for i in sample_list_all :
                if i not in out_list :
                    # for snps only in outgroup, modify GT to ./.
                    # modify GT to ./. if not 0/0
                    if i in sample_list :
                        if lh.genotype(i)['GT'] != '0/0' and sample_dict[i] != "MODERN" :
                            lh.genotype(i).data = lh.genotype(i).data._replace(GT='./.')
                    else :
                        if lh.genotype(i)['GT'] != '0/0' :
                            lh.genotype(i).data = lh.genotype(i).data._replace(GT='./.')
            tmp_ovcf.write_record(lh)
        elif f_out[0] == 2 and f_out[1] == 1 :
            # write a out group snp
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
    sample_file = args.sample


    # read in sample list
    sample_list = {}
    s_fh = open(sample_file, 'r')
    while 1 :
        s_l = s_fh.readline().strip().split('\t')

        if len(s_l) == 1 :
            break
        sample_list[s_l[0]] = s_l[1]

    # read in outgroup list , if True
    if args.out :
        out_list = args.out.split(',')
    else :
        out_list = False

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
        'child_vcf' : '',
        'filter_list' : '',
        'wdir' : '',
        'out_list' : False,
        'mod_gt' : False
    }

    arg_indi['wdir'] = start_dir + '/tmp_' + in_vcf
    arg_indi['filter_list'] = sample_list
    if out_list :
        arg_indi['out_list'] = out_list
    arg_indi['mod_gt'] = args.mod_gt
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
                           outfile=start_dir + '/' + 'dbSNP_' + in_vcf,
                           shell=True,
                           wait=True)
    reduce_child_vcf.run()

    if not args.no_rmdir :
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
    parse.add_argument('-s', '--sample',
                       help='sample list for filter, SAMPLE\\tGROUP',
                       required=True, type=str)
    parse.add_argument('--mod_gt',
                       help='turn on will modify GT to ./. for SNPs only in outgroup '
                            'if GT is not 0/0, make as UNKNOWN',
                       action='store_true', default=False)
    parse.add_argument('--out', help='outgroup list', type=str, default=False)
    parse.add_argument('--name', help="name affix to input vcf filename", type=str,
                       default=False)
    parse.add_argument('--continue_run', help="continue running without new split file", action='store_true',
                       default=False)
    parse.add_argument('--no_rmdir', help="running without rmdir", action='store_true',
                       default=False)


    args = parse.parse_args()

    v_f(args)


if __name__ == '__main__' :
    arg()
