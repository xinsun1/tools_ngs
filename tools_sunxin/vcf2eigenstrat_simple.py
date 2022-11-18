#!/usr/bin/python3

__author__ = 'sunxin'

import vcf
import sys, os

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


    fh = vcf.Reader(open(in_vcf), 'r')

    samples = fh.samples
    samples.sort()
    nsample = len(samples)

    out_snp = open(out_name + '.snp', 'w')
    out_geno = open(out_name + '.geno', 'w')

    for lh in fh:
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

            if gt_list[0] == ".":  # ALLOW MISSING DATA
                n = 9
            else:
                n = 0
                for j in gt_list:
                    if str(j) == '0':
                        n += 1

            o_geno_str = o_geno_str + str(n)

        print(o_geno_str, file=out_geno)
        print('\t'.join(o_snp_list), file=out_snp)

    out_geno.close()
    out_snp.close()

    o_indiv = open(out_name + '.ind', 'w')
    for i in range(len(samples)) :
        print(samples[i], file=o_indiv)

    o_indiv.close()


if __name__ == '__main__' :
    print('vcf format to eigenstrat format conversion for admixtools ')
    print('in_vcf, out_name, np')
    v2e(sys.argv[1], sys.argv[2], sys.argv[3])

