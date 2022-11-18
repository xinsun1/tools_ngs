#!/usr/bin/python3
__author__ = 'sunxin'


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

def get_gt(l_vcf, sample) :
    '''
    '''
    global iupac_dict

    tmp_seq = ''.join(l_vcf.genotype(sample)['GT'].split('/'))
    tmp_replace = tmp_seq.replace('0', str(l_vcf.REF)).replace('1', str(l_vcf.ALT[0]))
    gt = iupac_dict[tmp_replace]

    return gt

def gen_fa(vcf_file, fa, out_name, sample) :
    '''
    :param vcf:
     vcf file
    :param fa:
     in loci fa
    :param out_name:
     out file prefix name
    :return:
    '''

    vcf_fh = vcf.VCFReader(open(vcf_file, 'r'))

    ofh = open(out_name + '.fa', 'w')

    fa_fh = Fasta(fa)
    fa_list = fa_fh.keys()
    fa_dict = {}
    for i in fa_list :
        fa_dict[i] = False

    sample = sample


    c_chr = ''
    c_fa = ''
    c_pos = 0
    is_start = True

    for lh in vcf_fh :
        lh_chr = lh.CHROM
        lh_pos = int(lh.POS)        # 1-based

        if is_start :
            is_start = False
            c_chr = lh_chr


        if lh_chr == c_chr :
            # sample chr
            # update fa, pos
            c_fa += fa_fh[c_chr][c_pos:lh_pos - 1] # 0-based
            c_pos = lh_pos
        else :
            # new chr
            # update pre_chr, print pre_chr
            c_fa += fa_fh[c_chr][c_pos:]
            print(">" + c_chr, file=ofh)    # output fa
            print(c_fa, file=ofh)
            fa_dict[c_chr] = True           # remove fa from list
            # update fa, pos as start
            c_pos = 0
            c_chr = lh_chr
            c_fa = ''
            c_fa += fa_fh[c_chr][c_pos:lh_pos - 1]
            c_pos = lh_pos

        # update info for this SNP
        c_fa += get_gt(lh, sample)

    c_fa += fa_fh[c_chr][c_pos:]
    print(">" + c_chr, file=ofh)            # output fa
    print(c_fa, file=ofh)
    fa_dict[c_chr] = True                   # remove fa from list

    # if no SNP in chr
    for i in fa_list :
        if fa_dict[i] == False :
            print(">" + i, file=ofh)
            print(fa_fh[i], file=ofh)

    ofh.close()


if __name__ == '__main__':
    print('VCF format to FASTA file transformation. \n'
          'Per sample transformation for parallelism. \n'
          '\n'
          'Usage: \n'
          'vcf2fa.py IN_VCF IN_FA SAMPLE_ID OUT_NAME \n')

    gen_fa(vcf_file=sys.argv[1],
           fa=sys.argv[2],
           sample=sys.argv[3],
           out_name=sys.argv[4])

