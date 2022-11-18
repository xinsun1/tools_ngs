#!/usr/bin/python3

__author__ = 'sunxin'

'''

Raw vcf file filter generated from GATK.

Take input from vcf_filter_new_step2.py
Hard filter passed, MaxDP marked, Informative site allowed only.
This filter will included :
    mark snp with low DP as missing;
    remove snp with too much missing data.


Filter strategy :
1) Read sample population file for MODERN and AMO pop sample.
    SampleID    Group_Info(OUT,lowDP, Degraded, Modern, etc.) Sample_quality(fresh, degraded), MeanDP
3) Filter DP minimum :
                                    fresh           degraded
    Pesudohaploid                   >= 4               >= 1
    High quality (conserved)        >= 4               >= 4 
    High quality                    >= 4               >= 3
    
    max(1/3 of meanDP, minDP)

4) Missing allowed :
                                    fresh           degraded
    Pesudohaploid                     1                 18
    High quality (conserved)          1                 6
    High quality                      1                 15

    Missing allowed = 1(fresh) + degraded - expected / 2
        
5) Output SNP :
    DP_FILTER = 1
    MISSING_FILTER = 1 
    output

'''

from lib_split_file import *
import vcf
import os
import multiprocessing
from lib_cmd import *
import argparse


def line_filter(vl, sample_list, sample_dict, keep_out, minDP, maxMiss) :

    # generate filter DP for each sample
    filterDP_dict = {}
    for i in sample_list :
        if sample_dict[i]['quality'] == 'fresh' :
            filterDP_dict[i] = max(round(sample_dict[i]['meanDP'] / 3, 2), 4.0)
        else :
            filterDP_dict[i] = max(round(sample_dict[i]['meanDP'] / 3, 2), minDP)

    # Filter DP minimum, fresh >= 4 , degrade >= minDP, mark as missing
    for i in sample_list :
        if vl.genotype(i)['DP'] == None :
            vl.genotype(i).data = vl.genotype(i).data._replace(GT='./.')
            continue
        if float(vl.genotype(i)['DP']) < filterDP_dict[i] :
            vl.genotype(i).data = vl.genotype(i).data._replace(GT='./.')


    # calculate missing rate
    vl_freq = {'.' : 0,   # allele freq
               '0' : 0,
               '1' : 0}
    for i in sample_list :
        gt_list = vl.genotype(i)['GT'].split('/')
        vl_freq[str(gt_list[0])] += 1
        vl_freq[str(gt_list[1])] += 1
    miss_count = vl_freq['.'] / 2

    filter_result = 0
    if vl_freq['1'] != 0 :
        # is a snp
        if float(miss_count) <= maxMiss :
            filter_result = 1

    return filter_result


def v_f_indi(arg_indi) :
    os.chdir(arg_indi['wdir'])
    fh = vcf.Reader(open(arg_indi['child_vcf']), 'r')

    sample_dict = arg_indi['sample_info']
    sample_list = arg_indi['sample_info'].keys()

    keep_out = arg_indi['keep_out']

    ovcf = open('tmp_' + arg_indi['child_vcf'] + '.vcf', 'w')
    tmp_ovcf = vcf.Writer(ovcf, fh)

    for lh in fh :
        f_out = line_filter(lh, sample_list, sample_dict,
                            keep_out, arg_indi['minDP'], arg_indi['maxMiss'])

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
    sample_file = args.sample_info

    # read in sample info
    sample_info = {}
    s_fh = open(sample_file, 'r')
    while 1 :
        s_l = s_fh.readline().strip().split('\t')

        if len(s_l) == 1 :
            break
        sample_info[s_l[0]] = {
            'group' : str(s_l[1]),
            'quality' : str(s_l[2]),
            'meanDP' : float(s_l[3])
        }

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
        'sample_info'   : sample_info,
        'wdir'          : start_dir + '/tmp_' + in_vcf,
        'keep_out'      : False,
        'minDP'         : float(args.dp),
        'maxMiss'       : float(args.miss)
    }

    arg_indi['wdir'] = start_dir + '/tmp_' + in_vcf
    arg_indi['filter_list'] = sample_info
    if args.out :
        arg_indi['keep_out'] = True
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

    parse = argparse.ArgumentParser(prog = 'vcf_filter')

    parse.add_argument('-p', '--np',
                       help=' maximum number of pool workers, default= 10',
                       type=int, default=10)
    parse.add_argument('-i', '--vcf',
                       help='input vcf file', required=True, type=str)
    parse.add_argument('-s', '--sample_info',
                       help='sample information provided for filter.'
                            'Sample_ID\\tGROUP\\tSample_quality\\tmean_DP',
                       required=True, type=str)
    parse.add_argument('--out', help='Keep Outgroup SNP, default : TRUE', action='store_true', default=True)
    parse.add_argument('--name', help="name affix to input vcf filename, default=info_", type=str,
                       default='info_')
    parse.add_argument('--continue_run', help="continue running without new split file", action='store_true',
                       default=False)
    parse.add_argument('--no_rmdir', help="running without rmdir", action='store_true',
                       default=False)
    parse.add_argument('--dp', help='minDP for degraded samples', type=int)
    parse.add_argument('--miss', help='missing allowed for a loci', type=int)


    args = parse.parse_args()

    v_f(args)


if __name__ == '__main__' :
    arg()
