#!/usr/bin/python3
__author__ = 'sunxin'

import numpy as np
import vcf
import os



def dis(str1, str2) :
    # dis calculation between two GT
    # allow missing data, processed in pdis_indi

    # with no missing data
    l1 = str1.split('/')
    l2 = str2.split('/')

    n = 0
    for i in l1 :
        if i in l2 :
            n += 1
    for j in l2 :
        if j in l1 :
            n += 1

    if n == 0 :
        return 1
    elif n == 4 :
        return 0
    else :
        return 0.5

def pdis_indi(arg_indi) :

    os.chdir(arg_indi['wdir'])
    file = arg_indi['child_vcf']


    fh = vcf.Reader(open(file, 'r'))

    samples = fh.samples
    nsample = len(samples)

    dis_array = np.zeros((nsample, nsample))
    pair_array = np.zeros((nsample, nsample))


    for lh in fh :
        for i in range(nsample):
            gt1 = lh.genotype(samples[i])['GT']

            # if missing data, pair discard
            if gt1 == './.' :
                continue
            for j in range(i + 1, nsample) :
                gt2 = lh.genotype(samples[j])['GT']
                # if miss data, pair discard
                if gt2 == './.' :
                    continue

                dis_pair = dis(gt1, gt2)
                dis_array[i, j] += dis_pair
                pair_array[i, j] += 1

    return [samples, dis_array, pair_array]

