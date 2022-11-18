#!/usr/bin/python3
__author__ = 'sunxin'

'''
This script generate G-PhoCS seq file based on the loci given.

Last modified : 20181017
'''


import vcf, sys
from pyfasta import Fasta

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

def in_loci(loci_list, pos) :
    '''

    :param chr:
     query chr
    :param loci:
     loci list in chr 0based
    :param pos:
     query pos 1based
    :return:
    '''

    q_pos = pos - 1
    for i in range(0, len(loci_list)) :
        if q_pos < int(loci_list[i][0]) :
            return 0
        if q_pos < int(loci_list[i][1]) :
            return loci_list[i]
    return 0

def update_fa(sub_fa, l_vcf, seq_loci) :
    '''

    :param sub_fa:
     sub_fa_dict
    l_vcf
        vcf line
    seq_loci
        loci info for this seq

    '''

    global iupac_dict

    up_dict = {}
    pos = int(l_vcf.POS) - 1
    x_pos = pos - int(seq_loci[0]) + 1

    for i in sub_fa.keys() :
        seq = sub_fa[i]
        tmp_seq = ''.join(l_vcf.genotype(i)['GT'].split('/'))
        tmp_replace = tmp_seq.replace('0', str(l_vcf.REF)).replace('1', str(l_vcf.ALT[0]))
        gt = iupac_dict[tmp_replace]
        new_seq = seq[0 : x_pos - 1] + gt + seq[x_pos : 1000]
        up_dict[i] = new_seq

    return up_dict

def gen_fa(vcf_file, fa, bed, out_name) :
    '''

    :param vcf:
     vcf file
    :param fa:
     in loci fa
    :param bed:
     in loci bed file
    :param out_name:
     out file prefix name
    :return:
    '''

    vcf_fh = vcf.VCFReader(open(vcf_file, 'r'))

    ofh = open(out_name + '_gphocs_seq', 'w')

    fa_fh = Fasta(fa)

    bed_fh = open(bed, 'r')

    bed_dict = {}

    while 1 :
        bed_l = bed_fh.readline().strip().split('\t')

        if len(bed_l) == 1 :
            break

        if not bed_l[0] in bed_dict.keys() :
            bed_dict[bed_l[0]] = []

        bed_dict[bed_l[0]].append([bed_l[1], bed_l[2]])

    samples = vcf_fh.samples
    nsample = len(samples)

    fa_dict = {}
    for i in fa_fh.keys() :
        fa_dict[i] = {}
        for j in range(0, nsample) :
            fa_dict[i][samples[j]] = str(fa_fh[i])

    for lh in vcf_fh :
        lh_chr = lh.CHROM
        lh_pos = int(lh.POS)

        if not lh_chr in bed_dict.keys() :
            continue

        is_loci = in_loci(bed_dict[lh_chr], lh_pos)

        if is_loci == 0 :
            continue

        up = update_fa(sub_fa=fa_dict[str(lh_chr) + ':' + str(is_loci[0]) +
                                      '-' + str(is_loci[1])],
                       l_vcf=lh,
                       seq_loci=is_loci)

        fa_dict[str(lh_chr) + ':' + str(is_loci[0]) + '-' + str(is_loci[1])] = up

    print(len(fa_dict.keys()), file=ofh)
    print('', file=ofh)

    for i in fa_dict.keys() :
        print('\t'.join([i, str(nsample), str(1000)]), file=ofh)

        for j in range(0, nsample) :
            print('\t'.join([samples[j], fa_dict[i][samples[j]]]), file=ofh)
        print('', file=ofh)

    ofh.close()
    bed_fh.close()

if __name__ == '__main__':
    print('in_vcf, loci_fa, loci_bed, out_name')
    gen_fa(vcf_file=sys.argv[1],
           fa=sys.argv[2],
           bed=sys.argv[3],
           out_name=sys.argv[4])




























