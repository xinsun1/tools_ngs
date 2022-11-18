#!/bin/python3
__author__ = 'sunxin'

# lib for pairwise sequence alignment

import parasail
import argparse
from lib_fastaq import *
import time

def trace_back(ar, seq1, seq2, gopen, gextend, match, mismatch) :

    # function to trace back a DP table from parasail alignment result
    # return two aligned sequence

    dim = ar.shape      # DP table shape
    nseq1 = int(dim[0]) - 1      # seq1 index
    nseq2 = int(dim[1]) - 1      # seq2 index

    outseq1 = ''
    outseq2 = ''

    score = lambda a : match if a[0]==a[1] else mismatch

    while nseq1 > 0 and nseq2 > 0 :
        if ar[nseq1, nseq2] - ar[nseq1 - 1, nseq2 - 1] == score([seq1[nseq1], seq2[nseq2]]) :
            # diagonal situation

            outseq1 += seq1[nseq1]
            outseq2 += seq2[nseq2]

            nseq1 -= 1
            nseq2 -= 1

        else :

            found = False
            for i in range(1,min(nseq1, nseq2) + 1) :
                # min(nseq1, nseq2)

                if ar[nseq1, nseq2] - ar[nseq1 - i, nseq2] == (i - 1) * gextend + gopen :
                    # vertical i gap, gap open in i

                    outseq1 += seq1[(nseq1 - i + 1): (nseq1 + 1)][::-1]
                    outseq2 += '-' * i

                    nseq1 -= i
                    found = True

                elif ar[nseq1, nseq2] - ar[nseq1, nseq2 - i] == (i - 1) * gextend + gopen :
                    # horizontal i gap, gap open in i

                    outseq1 += '-' * i
                    outseq2 += seq2[(nseq2 - i + 1) : (nseq2 + 1)][::-1]

                    nseq2 -= i
                    found = True

                elif ar[nseq1, nseq2] - ar[nseq1 - i, nseq2] == i * gextend and \
                    ar[nseq1 - i, nseq2] - ar[nseq1 - i - 1, nseq2] == (gopen or gextend) :
                    # vertical i gap, gap open not in i

                    outseq1 += seq1[(nseq1 - i + 1): (nseq1 + 1)][::-1]
                    outseq2 += '-' * i

                    nseq1 -= i
                    found = True

                elif ar[nseq1, nseq2] - ar[nseq1, nseq2 - i] == i * gextend and \
                    ar[nseq1, nseq2 - i] - ar[nseq1, nseq2 - i - 1] == (gopen or gextend) :
                    # horizontal i gap, gap open not in i

                    outseq1 += '-' * i
                    outseq2 += seq2[(nseq2 - i + 1) : (nseq2 + 1)][::-1]

                    nseq2 -= i
                    found = True

                if found == True :
                    break

            if found == False :
                if max(nseq1, nseq2) == nseq1 :
                    for i in range(nseq2 + 1, nseq1) :
                        if ar[nseq1, nseq2] - ar[nseq1 - i, nseq2] == (i - 1) * gextend + gopen :
                            # vertical i gap, gap open in i

                            outseq1 += seq1[(nseq1 - i + 1): (nseq1 + 1)][::-1]
                            outseq2 += '-' * i

                            nseq1 -= i
                            found = True
                        elif ar[nseq1, nseq2] - ar[nseq1 - i, nseq2] == i * gextend and \
                            ar[nseq1 - i, nseq2] - ar[nseq1 - i - 1, nseq2] == (gopen or gextend) :
                            # vertical i gap, gap open not in i

                            outseq1 += seq1[(nseq1 - i + 1): (nseq1 + 1)][::-1]
                            outseq2 += '-' * i

                            nseq1 -= i
                            found = True

                        if found == True :
                            break
                else :
                    for i in range(nseq1 + 1, nseq2) :
                        if ar[nseq1, nseq2] - ar[nseq1, nseq2 - i] == (i - 1) * gextend + gopen :
                            # horizontal i gap, gap open in i

                            outseq1 += '-' * i
                            outseq2 += seq2[(nseq2 - i + 1) : (nseq2 + 1)][::-1]

                            nseq2 -= i
                            found = True
                        elif ar[nseq1, nseq2] - ar[nseq1, nseq2 - i] == i * gextend and \
                            ar[nseq1, nseq2 - i] - ar[nseq1, nseq2 - i - 1] == (gopen or gextend) :
                            # horizontal i gap, gap open not in i

                            outseq1 += '-' * i
                            outseq2 += seq2[(nseq2 - i + 1) : (nseq2 + 1)][::-1]

                            nseq2 -= i
                            found = True

                        if found == True :
                            break

            if found == False :
                print(seq1)
                print(seq2)
                print(nseq1)
                print(nseq2)
                print(ar)
                print (ar[nseq1 - 1 : nseq1+1, nseq2 - 1: nseq2+1])
                print(ar[75:85, 75:85])

                raise ValueError


    if nseq1 == 0 and nseq2 == 0 :
        # [0,0]

        outseq1 += seq1[0]
        outseq2 += seq2[0]

    elif nseq1 == 0 :
        # seq1 used out, horizontal border
        # let gap in the border

        outseq1 += (seq1[0] + '-' * nseq2)
        outseq2 += seq2[0 : nseq2 + 1][::-1]

    elif nseq2 == 0 :
        # seq2 used out, vertical border
        # let gap in the border

        outseq1 += seq1[0 : nseq1 + 1][::-1]
        outseq2 += (seq2[0] + '-' * nseq1)

    else :
        raise ValueError

    return [outseq1[::-1], outseq2[::-1]]


def merge_seq(align_seq1, align_seq2) :

    out_seq = ''
    is_N = False
    for i in range(0,len(align_seq1)) :
        if align_seq1[i] == align_seq2[i] :
            out_seq += align_seq1[i]
        elif align_seq1[i] == '-' :
            out_seq += align_seq2[i]
        elif align_seq2[i] == '-' :
            out_seq += align_seq1[i]
        else :
            out_seq += 'N'
            is_N = True

    return [out_seq , is_N]


def reverse(seq) :
    'reverse seq'

    a = {
        'A' : 'T',
        'G' : 'C',
        'C' : 'G',
        'T' : 'A',
        'N' : 'N'
         }

    out = ''
    for i in range(0,len(seq)) :
        out += a[seq[::-1][i]]

    return out


def align(seq_file, match = 2, mismatch = -1, gapopen = -5, gapextend = -1) :

    match = match
    mismatch =  mismatch
    gapopen = gapopen
    gapextend = gapextend
    out_name = seq_file

    matrix = parasail.matrix_create("AGCTN", match, mismatch)

    seq1_file = fasta(file = seq_file + '_1.fq', is_fastq = True, is_pair = True)
    seq2_file = fasta(file = seq_file + '_2.fq', is_fastq = True, is_pair = True)

    out_fh_aligned = open(out_name + '_aligned.fa', 'w')
    out_fh_un1 = open(out_name + '_unaligned_1.fq', 'w')
    out_fh_un2 = open(out_name + '_unaligned_2.fq', 'w')

    while 1 :

        seq1 = seq1_file.next()
        seq2 = seq2_file.next()

        if seq1 == 0 :
            break

        if str(seq1.name) != str(seq2.name) :
            raise ValueError(0)

        result = parasail.nw_stats_table_striped_64(seq1.line, reverse(seq2.line),
                                                            -gapopen, -gapextend, matrix)

        alignseq = trace_back(ar = result.score_table,
                              seq1 = seq1.line,
                              seq2 = reverse(seq2.line),
                              gopen = gapopen,
                              gextend = gapextend,
                              match = match,
                              mismatch = mismatch)

        merge_result = merge_seq(alignseq[0], alignseq[1])

        if merge_result[1] :
            seq1.print(out_fh_un1)
            seq2.print(out_fh_un2)
        else :
            print('>' + seq1.name[1:len(seq1.name)], file=out_fh_aligned)
            print(merge_result[0], file=out_fh_aligned)

    return out_name

def arg():

    # arguements management

    parser = argparse.ArgumentParser(prog='Pair to Single')

    parser.add_argument('-s', '--pair', help='input pair-end name, PAIR_1.fq, PAIR_2.fq', required=True, type=str)
    parser.add_argument('-m', '--match', help='match score, default=2', default=2, type=int)
    parser.add_argument('-mis', '--mismatch', help='mismatch score, default=-1', default=-1, type=int)
    parser.add_argument('-gop', '--gapopen', help='gap open score, default=-5', default=-5, type=int)
    parser.add_argument('-gex', '--gapextend', help='gap extend score, default=-1', default=-1, type=int)

    args = parser.parse_args()

    align(seq_file= args.pair,
          match=args.match,
          mismatch=args.mismatch,
          gapopen=args.gapopen,
          gapextend=args.gapextend)



if __name__ == '__main__' :
    start_time = time.time()
    arg()
    finish_time = time.time()
    print('Finished.')
    print('Total time :' + str(finish_time - start_time))





