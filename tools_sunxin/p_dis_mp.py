#!/usr/bin/python3
__author__ = 'sunxin'

import multiprocessing as mp
import numpy as np
import vcf
import sys


# dis calculation between two GT
def dis(str1, str2) :
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

def getp(lh) :
    global samples, dis_array, nsample

    lh_array = np.zeros((nsample, nsample))

    for i in range(nsample):
        for j in range(i + 1, nsample) :
            lh_array[i, j] += dis(lh.genotype(samples[i])['GT'], lh.genotype(samples[j])['GT'])

    dis_array += lh_array


if __name__ == '__main__' :
    # arg1 = filename
    # arg2 = np

    file = sys.argv[1]
    nprocess = sys.argv[2]

    fh = vcf.Reader(open(file, 'r'))

    global samples, dis_array, nsample
    samples = fh.samples
    nsample = len(samples)

    dis_array = np.zeros((nsample, nsample))

    with mp.Pool(processes=int(nprocess)) as pool :
        pool.map(getp, fh)

    outfile = open('dis_array' ,'w')

    print('\t'.join(samples), file=outfile)

    for i in range(dis_array.shape[0]) :
        print('\t'.join(map(str, list(dis_array[i,]))), file=outfile)





