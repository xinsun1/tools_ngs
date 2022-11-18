#!/usr/bin/python3

__author__ = 'sunxin'

import subprocess
import sys

# sort vcf file by order from chr_order
# too slow

def sort_ref(in_vcf, chr_order):

    fh = open(chr_order, 'r')

    chr_list = []

    while 1 :
        lh = fh.readline().strip()

        if len(lh) == 0 :
            break

        chr_list.append(str(lh))

    fh.close()

    vfh = open(in_vcf, 'r')
    while 1 :
        lh = vfh.readline().strip()
        if len(lh) == 0 :
            break

        llh = lh.split('\t')

        if llh[0][0] == '#' :
            with open('head_' + in_vcf, 'a') as ofh :
                print(lh, file=ofh)
            ofh.close()
        else :
            with open('tmp_' + llh[0], 'a') as ofh :
                print(lh, file=ofh)
            ofh.close()

    chr_olist = ['tmp_' + i for i in chr_list]

    with open('sort_' + in_vcf, 'w') as ofh :
        subprocess.call('cat head_' + in_vcf + ' ' + ' '.join(chr_olist),
                        stdout=ofh,
                        shell=True)

    vfh.close()




if __name__ == '__main__' :
    print('in_vcf, chr_order')
    sort_ref(sys.argv[1], sys.argv[2])
