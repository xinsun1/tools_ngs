#!/usr/bin/python3
__author__ = 'sunxin'

from lib_split_file import *
import vcf
import sys, os
import multiprocessing
from lib_cmd import *


iupac_dict = {'AA' : 'A',
              'GG' : 'G',
              'CC' : 'C',
              'TT' : 'T',
              'AG' : 'R',
              'AC' : 'M',
              'AT' : 'W',
              'GA' : 'R',
              'GC' : 'S',
              'GT' : 'K',
              'CA' : 'M',
              'CG' : 'S',
              'CT' : 'Y',
              'TA' : 'W',
              'TG' : 'K',
              'TC' : 'Y',
              '..' : 'N',  # ALLOW MISSING DATA
              }

def get_gt(l_vcf, sample) :
    '''
    '''
    global iupac_dict

    tmp_seq = ''.join(l_vcf.genotype(sample)['GT'].split('/'))
    tmp_replace = tmp_seq.replace('0', str(l_vcf.REF)).replace('1', str(l_vcf.ALT[0]))
    gt = iupac_dict[tmp_replace]

    return gt

def v2p_indi(arg_indi) :
    os.chdir(arg_indi['wdir'])
    fh = vcf.Reader(open(arg_indi['child_vcf']), 'r')

    samples = fh.samples
    samples.sort()
    nsample = len(samples)

    dict_tmp = {}
    for i in range(nsample) :
        dict_tmp[samples[i]] = []

    for lh in fh :
        for i in range(nsample) :
            dict_tmp[samples[i]].append(get_gt(lh, samples[i]))

    list_return = []
    for i in range(nsample) :
        list_return.append(''.join(dict_tmp[samples[i]]))

    return [samples, list_return]


def v2p(in_vcf, out_name, np) :
    '''

    :param in_vcf:
    in vcf file name
    :param out_name:
     out phylip file name
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
    re_v2p = mp_pool.map(v2p_indi, arg_list)

    # reduce result

    os.chdir(start_dir)
    out_list = []

    samples = re_v2p[0][0]
    for i in range(len(samples)) :
        out_list.append([samples[i], ''])

    for i in range(len(re_v2p)) :
        for j in range(len(samples)) :
            out_list[j][1] = out_list[j][1] + re_v2p[i][1][j]

    ofh = open(out_name + '.phylip', 'w')

    print('\t'.join([str(len(samples)), str(len(out_list[0][1]))]), file=ofh)

    for i in range(len(out_list)) :
        print('\t'.join(out_list[i]), file=ofh)

    ofh.close()
    in_vcf_fh.destory()

if __name__ == '__main__' :
    print('vcf format to phylip format conversion ')
    print('in_vcf, out_name, np')
    v2p(sys.argv[1], sys.argv[2], sys.argv[3])
