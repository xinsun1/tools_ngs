#!/usr/bin/python3

__author__ = 'sunxin'

from pyfasta import Fasta
import sys

def is_CpG(fa, GC=0.5, len=200, ) :
    '''

    :param chr:
     scaffold num
    :param st:
     start pos, 0 based
    :param fa:
     fasta fh by pyfasta
    :param per_N:
     percentage of N allowed
    :return:
    '''

    seq = fa[str(chr)][int(st) : int(st) + 1000]
    seq_N = float(seq.upper().count('N')) / 1000.0
    if seq_N > float(per_N) :
        return False
    else :
        return True


def get_loci(fa, out_name, chr_len, per_N=0.05) :
    '''

    :param fa:
     in fasta file
    :param mark_file:
     in mark file bed format, 0 based
    :param chr_len:
     in chr length file, chr'\t'len
    :return:
    '''

    fa_fh = Fasta(fa)

    chr_fh  = open(chr_len, 'r')

    fa_loci_fh = open('loci_' + out_name + '.fa', 'w')
    bed_loci_fh = open('loci_' + out_name + '.bed', 'w')

    chr_dict = {}

    while 1 :
        chr_fh_l = chr_fh.readline().strip().split('\t')

        if len(chr_fh_l) == 1 :
            break

        chr_dict[chr_fh_l[0]] = chr_fh_l[1]

    for i in chr_dict.keys() :
        st_pos = random.randint(0, 50000)

        pos = st_pos
        while pos + 1000 <= int(chr_dict[i]) :
            if is_loci(chr=i, st=pos, fa=fa_fh, per_N=per_N) :
                print('>' + str(i) + ':' + str(pos) + '-' + str(int(pos) + 1000),
                      file=fa_loci_fh)
                print(fa_fh[str(i)][int(pos) : int(pos) + 1000],
                      file=fa_loci_fh)
                print('\t'.join([str(i), str(pos), str(pos + 1000)]),
                      file=bed_loci_fh)

                pos = pos + 1000 + 50000
            else :
                pos = pos + 1

    chr_fh.close()
    fa_loci_fh.close()
    bed_loci_fh.close()


if __name__ == '__main__':
    print('Find CpG islands in the genome' \
          'Reference https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=660100749_awTzDUQrFX6aB7XnZsqNrH9hF6Gu&c=chrA2&g=cpgIslandSuper')
    print('In Fasta, out_name, chr_len file')
    get_loci(fa=sys.argv[1],
             out_name=sys.argv[2],
             chr_len=sys.argv[3],
             per_N=sys.argv[4])


