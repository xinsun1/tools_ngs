#!/usr/bin/python3

__author__ = 'sunxin'

'''
vcf to allele frequency spectrum claculation

Multi-threads supported.
'''

from lib_split_file import *
import vcf
import sys, os
import multiprocessing
from lib_cmd import *
import numpy


def v2a_pop(arg_indi) :
    vcf_file = arg_indi['child_vcf']        # child vcf file
    pop_dict = arg_indi['pop_dict']         # population information
    wdir = arg_indi['wdir']                 # work directory
    use_out = int(arg_indi['use_out'])           # whether use OUT 0/1

    os.chdir(wdir)
    fh = vcf.Reader(open(vcf_file), 'r')

    pop_list = list(pop_dict.keys())
    pop_list.sort()
    len_pop = []
    dim_pop = []
    pop_list_mod = []
    if use_out == 1 :
        num_pop = len(pop_list)         # include OUT,  number of populations for calculate AFS
    else :
        num_pop  = len(pop_list) - 1    # not include OUT


    for i in range(len(pop_list)) :
        if pop_list[i] == 'OUT' :
            if use_out == 0 :
                continue
        len_pop.append(len(pop_dict[pop_list[i]]))
        dim_pop.append(2 * len(pop_dict[pop_list[i]]) + 1)
        pop_list_mod.append(pop_list[i])

    pop_list_mod.sort()

    ar_indi = numpy.zeros(dim_pop, dtype=numpy.dtype(int))

    sum_count = 0       # number of snps used
    sum_drop = 0        # number of snps dropped
    for lh in fh :
        # is ancestral allele
        ancestral_allele = ''
        for i in pop_dict['OUT'] :
            gt_i = lh.genotype(i).gt_type
            if gt_i == 0 or gt_i == 2 :
                if ancestral_allele == '' :
                    ancestral_allele = gt_i
                else :
                    if ancestral_allele != gt_i :
                        ancestral_allele = False
                        break
            else :
                ancestral_allele = False
                break

        if not ancestral_allele :
            sum_drop += 1               # ancestral allele undetermined
            continue

        derived_freq_count = []                 # list to store allele count
        is_na = 0
        for i in range(num_pop) :
            # if ./. exist
            if is_na == 1 :
                break

            pop_derived_count = 0
            for j in pop_dict[pop_list_mod[i]] :
                gt_j = lh.genotype(j).gt_type   # individual gt_type 0, 1, 2, None
                if gt_j == None :
                    is_na = 1
                    break
                if ancestral_allele == 0 :
                    pop_derived_count += int(gt_j)
                else :
                    pop_derived_count += (2 - int(gt_j))
            derived_freq_count.append(pop_derived_count)

        if is_na == 1 :
            sum_drop += 1               # ./. exist for individual
            continue
        else :
            # if everything looks fine
            ar_indi[tuple(derived_freq_count)] += 1
            sum_count += 1

    return [ar_indi, sum_count, sum_drop]


def v2a(in_vcf, pop_file, use_out, out_name, np) :
    '''

    :param in_vcf:
    in vcf file name
    :param out_name:
     out phylip file name
    np :
     numbef of threads
    :return:
    '''

    #### read pop file
    pop_fh = open(pop_file, 'r')
    pop_dict = {}
    while 1 :
        pop_lh = pop_fh.readline().strip().split(' ')

        if len(pop_lh) != 2 :
            break

        if pop_lh[1] in pop_dict.keys() :
            pop_dict[pop_lh[1]].append(pop_lh[0])
        else :
            pop_dict[pop_lh[1]] = [pop_lh[0]]

    #### in vcf split
    start_dir = os.getcwd()
    in_vcf_fh = split_file(file=in_vcf,
                           np=np,
                           wdir=start_dir)
    child_vcf_list = in_vcf_fh.split_head()

    #### child file mp run
    arg_list = []
    arg_indi  = {
        'child_vcf' : '',
        'pop_dict'  : '',
        'wdir'      : '',
        'use_out'   : 0,
    }

    arg_indi['wdir'] = start_dir + '/tmp_' + in_vcf
    arg_indi['pop_dict'] = pop_dict
    arg_indi['use_out'] = int(use_out)
    for i in child_vcf_list :
        arg_indi['child_vcf'] = i
        arg_app = arg_indi.copy()
        arg_list.append(arg_app)
    mp_pool = multiprocessing.Pool(int(np))
    re_v2p = mp_pool.map(v2a_pop, arg_list)

    #### reduce result

    os.chdir(start_dir)

    # get population order
    pop_list = list(pop_dict.keys())
    pop_list.sort()
    pop_list_mod = []                       # consider the input from USE_OUT
    for i in range(len(pop_list)) :
        if pop_list[i] == 'OUT':
            if int(use_out) == 0 :
                continue
        pop_list_mod.append(pop_list[i])
    pop_list_mod.sort()

    ar_total = re_v2p[0][0].copy()
    sum_count = re_v2p[0][1]
    sum_drop = re_v2p[0][2]

    for i in range(1, len(re_v2p)) :
        ar_total += re_v2p[i][0]
        sum_count += re_v2p[i][1]
        sum_drop += re_v2p[i][2]

    ofh = open(out_name, 'w')
    print('Total SNPs used for calculation : ' + str(sum_count))
    print('Total SNPs dropped : ' + str(sum_drop))
    print('Population order :')
    print('\t'.join(pop_list_mod))

    dim_ar = ar_total.shape

    if len(dim_ar) == 2 :
        out_ar_len = dim_ar[1]              # set output ar len for 2D
    else :
        out_ar_len = 20                     # default output ar_len for multiD : 20

    out_len_index = 0
    out_line = []                           # int as str to store in the list
    for i in numpy.nditer(ar_total) :
        if out_len_index == out_ar_len :
            print('\t'.join(out_line), file=ofh)
            out_len_index=0
            out_line = []
        out_len_index += 1
        out_line.append(str(i))
    print('\t'.join(out_line), file=ofh)

    ofh.close()
    in_vcf_fh.destory()

if __name__ == '__main__' :
    print('vcf file to allele frequency spectrum claculation. \n'
          'Usage:'
          'vcf2afs.py in_vcf_file population_list use_out(1/0) out_file_name Number_threads'
    )

    if len(sys.argv) == 1 :
        print('Please input')
    else :
        v2a(sys.argv[1],        # vcf file
            sys.argv[2],        # population list of individuals
            sys.argv[3],        # use outgroup or not
            sys.argv[4],        # out file name
            sys.argv[5])        # Number of threads to use
