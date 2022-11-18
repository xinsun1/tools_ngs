#!/usr/bin/python3

__author__ = 'sunxin'

'''

Raw vcf file filter generated from GATK.

Take input from bcftools filtered vcf with only biallelic SNP allowed.
bcftooll remove indel and allow biallelic site.
Hard filter passed


Filter strategy :
1) Read sample population file for MODERN and AMO pop sample.
    SampleID    Group_Info(OUT,lowDP, Degraded, Modern, etc.) Sample_quality(fresh, degraded), MeanDP
2) Read sample information of fresh and degraded sample.
3) Filter DP maximum : fresh 3x meanDP, degraded 4x meanDP. 
   Mark failed SNPs as unknown.
4) Informative filter :
        Check if alt allele exist in Ingroup samples; return IS_INGROUP_SNP_SIGNAL
        
        fresh sample : at least one sample with MinDP >= 4;
        Degraded with USER treatment and no significant C2T : at least one sample with MinDP >= 4 (RUSA21, pti220);
        Degraded with C2T : Min 3 samples and with 1 with DP >= 4;
        return INGROUP_SNP_FILTER_SIGNAL
        
        Outgroup : if allowed, MinDP >= 8, 1/3 meanDP;
        return OUTSNP_FILTER_SIGNAL
        
5) Output SNP :
        INGROUP_SNO_FILTER == pass; OUTPUT
        
        OUTGROUP_SNP_FILTER == pass && IS_INGROUP_SNP == failed; OUTPUT
        modify lowDP sample GT to ./. if not 0/0  

        
Special attention :
MinDP was not filtered for any sample.
Further filter will included :
    mark snp with low DP as missing;
    remove snp with too much missing data.

'''

from lib_split_file import *
import vcf
import os
import multiprocessing
from lib_cmd import *
import argparse


def line_filter(vl, sample_list, sample_dict, keep_out) :

    out_list = []
    in_list = []
    in_list_nolowDP = []
    for i in sample_list :
        if sample_dict[i]['group'] == 'OUT' :
            out_list.append(i)
        else :
            in_list.append(i)
            if sample_dict[i]['group'] !=  'lowDP' :
                in_list_nolowDP.append(i)

    # Filter DP maximum, fresh 3x , degrade 4x, mark as missing
    for i in sample_list :
        if vl.genotype(i)['DP'] == None :
            vl.genotype(i).data = vl.genotype(i).data._replace(GT='./.')
            continue
        if sample_dict[i]['quality'] == 'fresh' :
            if float(vl.genotype(i)['DP']) > (sample_dict[i]['meanDP'] * 3) :
                vl.genotype(i).data = vl.genotype(i).data._replace(GT='./.')
        elif sample_dict[i]['quality'] == 'degraded' :
            if float(vl.genotype(i)['DP']) > (sample_dict[i]['meanDP'] * 4) :
                vl.genotype(i).data = vl.genotype(i).data._replace(GT='./.')

    # check if snp exist in ingroup samples
    vl_freq = {'.' : 0,   # allele freq
               '0' : 0,
               '1' : 0}
    for i in in_list :
        gt_list = vl.genotype(i)['GT'].split('/')
        vl_freq[str(gt_list[0])] += 1
        vl_freq[str(gt_list[1])] += 1

    if vl_freq['1'] == 0 :
        is_ingroup_snp = 0
    else :
        is_ingroup_snp = 1

    #  Ingroup SNP filter
    ingroup_snp_filter = 0
    degrade_sample_count = 0
    degrade_sample_DP_pass = 0

    for i in in_list_nolowDP :
        gt = vl.genotype(i)['GT']
        dp = vl.genotype(i)['DP']

        if gt == '0/1' or gt == '1/1' :
            if dp == None :
                continue
            if float(dp) >= 4.0 :
                if sample_dict[i]['group'] == 'Modern' :
                    ingroup_snp_filter = 1
                    break
                elif sample_dict[i]['group'] == 'Degraded' :
                    degrade_sample_DP_pass = 1
                    degrade_sample_count += 1
            elif sample_dict[i]['group'] == 'Degraded' :
                degrade_sample_count += 1

        if degrade_sample_count >= 3 and degrade_sample_DP_pass == 1 :
            ingroup_snp_filter = 1
            break

    # Outgroup snp filter
    outgroup_snp_filter = 0

    if keep_out :
        for i in out_list :
            gt = vl.genotype(i)['GT']
            dp = vl.genotype(i)['DP']

            if gt == '0/1' or gt == '1/1':
                if dp == None :
                    continue
                if float(dp) >= 8.0 :
                    outgroup_snp_filter = 1
                    break

    return [is_ingroup_snp, ingroup_snp_filter, outgroup_snp_filter]


def v_f_indi(arg_indi) :
    os.chdir(arg_indi['wdir'])
    fh = vcf.Reader(open(arg_indi['child_vcf']), 'r')

    sample_dict = arg_indi['sample_info']
    sample_list = arg_indi['sample_info'].keys()

    keep_out = arg_indi['keep_out']


    ovcf = open('tmp_' + arg_indi['child_vcf'] + '.vcf', 'w')
    tmp_ovcf = vcf.Writer(ovcf, fh)

    for lh in fh :
        f_out = line_filter(lh, sample_list, sample_dict, keep_out)
        # f_out = [IS_INGROUP_SNP, INGROUP_SNP_FILTER, OUTGROUP_SNP_filter]
        if f_out[1] == 1 :
            # is ingroup snp
            tmp_ovcf.write_record(lh)
        elif f_out[0] == 0 and f_out[2] == 1 :
            # not ingroup snp and is outgroup snp, modify lowDP sample GT to ./. if not 0/0
            for i in sample_list :
                if sample_dict[i]['group'] == 'lowDP':
                    lh.genotype(i).data = lh.genotype(i).data._replace(GT='./.')
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
        'child_vcf' : '',
        'sample_info' : sample_info,
        'wdir' : '',
        'keep_out' : False
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


    args = parse.parse_args()

    v_f(args)


if __name__ == '__main__' :
    arg()
