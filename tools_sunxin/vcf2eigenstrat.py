#!/usr/bin/python3

__author__ = 'sunxin'

from lib_split_file import *
import vcf
import sys, os
import multiprocessing
from lib_cmd import *


def v2e_indi(arg_indi) :
    os.chdir(arg_indi['wdir'])
    fh = vcf.Reader(open(arg_indi['child_vcf']), 'r')

    samples = fh.samples
    samples.sort()
    nsample = len(samples)

    tmp_snp = open(arg_indi['child_vcf'] + '.snp', 'w')
    tmp_geno = open(arg_indi['child_vcf'] + '.geno', 'w')

    for lh in fh :
        o_geno_str = ''
        o_snp_list = []

        # snp[ID, chr, g_dis, physic dis, REF, ALT]
        o_snp_list.append(str(lh.CHROM) + ':' + str(lh.POS))
        o_snp_list.append(str('1'))
        o_snp_list.append(str('0'))
        o_snp_list.append(str('0'))
        o_snp_list.append(str(lh.REF[0]))
        o_snp_list.append(str(lh.ALT[0]))

        # o_geno_str num of ref

        for i in range(0, nsample):
            gt_list = lh.genotype(samples[i])['GT'].split('/')

            if gt_list[0] == "." :    # ALLOW MISSING DATA
                n = 9
            else :
                n = 0
                for j in gt_list :
                    if str(j) == '0' :
                        n += 1

            o_geno_str = o_geno_str + str(n)

        print(o_geno_str, file=tmp_geno)
        print('\t'.join(o_snp_list), file=tmp_snp)

    tmp_geno.close()
    tmp_snp.close()

    return samples

def v2e(in_vcf, out_name, np) :
    '''

    :param in_vcf:
    in vcf file name
    :param out_name:
     out eigenstrat file name
    np :
     numbef of threads
    :return:
    '''


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
    re_samples = mp_pool.map(v2e_indi, arg_list)

    # reduce result
    child_snp_file_list = []
    child_geno_file_list = []

    for i in child_vcf_list :
        child_snp_file_list.append(i + '.snp')
        child_geno_file_list.append(i + '.geno')


    reduce_child_snp = cmd_out(cmd_line='cat CMDARG1',
                           in_p=[' '.join(child_snp_file_list)],
                           work_dir=start_dir + '/tmp_' + in_vcf,
                           outfile=start_dir + '/' + out_name + '.snp',
                           shell=True,
                           wait=True)
    reduce_child_snp.run()

    reduce_child_geno = cmd_out(cmd_line='cat CMDARG1',
                           in_p=[' '.join(child_geno_file_list)],
                           work_dir=start_dir + '/tmp_' + in_vcf,
                           outfile=start_dir + '/' + out_name + '.geno',
                           shell=True,
                           wait=True)
    reduce_child_geno.run()

    samples = re_samples[0]

    os.chdir(start_dir)
    o_indiv = open(out_name + '.ind', 'w')
    for i in range(len(samples)) :
        print(samples[i], file=o_indiv)

    o_indiv.close()

    in_vcf_fh.destory()

if __name__ == '__main__' :
    print('vcf format to eigenstrat format conversion for admixtools ')
    print('in_vcf, out_name, np')
    v2e(sys.argv[1], sys.argv[2], sys.argv[3])

