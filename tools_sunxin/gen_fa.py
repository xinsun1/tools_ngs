__author__ = 'sunxin'

from pyfasta import Fasta
import sys

def gen_fa(in_fa, out_fa, seq_len, gap = 5) :
    fh = Fasta(in_fa)
    ofh = open(out_fa, 'w')
    seq_len = int(seq_len)
    gap = int(gap)

    for i in fh.keys() :
        for j in range(0, len(fh[i]) - int(seq_len), int(gap)) :
            seq = str(fh[i][j : j + int(seq_len)]).upper()
            if not 'N' in seq :
                print('>' + str(i) + ':' + str(j + 1) + '-' + str(j + seq_len + 1), file = ofh)
                print(seq, file=ofh)


if __name__ == '__main__':
    print('input file, output file, seq_len, gap')
    gen_fa(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


