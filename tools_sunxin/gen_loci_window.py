#!/usr/bin/python3

__author__ = 'sunxin'

'''
This script generate window based loci for phylogeny and divergence time estimate.
Output philip file at the same time.

IN :
    fasta files , already marked
    bed file for chr used as index
    window size, step size, missing rate for each loci
    
Out :
    loci bed
    philip file for each loci
    
1) Read in bed file and fasta file
2) Start move, estimate N rate.

    Pass : out put to bed, output philip file.
           move a step size
           continue
    Fail : move a small step
           continue

File structure :
----|----Start_Dir
        |----output_name_dir
            |----loci_output_name.bed
            |----scaffold_start_end
                |----scaffold1_start_end.philip
            |----...
'''

from pyfaidx import Fasta
import argparse, os, subprocess, time



def gen_loci(args) :
    '''
    '''

    # parse arguments
    fasta_file_name = args.fa
    bed_file        = args.bed
    chr_str         = args.chr
    window_size     = args.window
    step_size       = args.step
    small_size      = args.small        # small step size
    n_rate          = args.n_rate
    out_name        = args.out
    start_dir       = os.getcwd()

    # read in fasta file
    sample_list = []
    fa_dict = {}

    with open(fasta_file_name, 'r') as ffn_fh :
        for line in ffn_fh :
            line_list = line.strip().split(' ')
            sample_list.append(str(line_list[0]))
            fa_dict[str(line_list[0])] = Fasta(str(line_list[1]))
    sample_list.sort()

    # read in bed file
    chr_list = []
    chr_dict = {}
    if chr_str :
        chr_str_list = chr_str.split('_')
        chr_list.append(chr_str_list[0])
        chr_dict[chr_str_list[0]] = int(chr_str_list[2])

    else :
        with open(bed_file, 'r') as bed_fh :
            for line in bed_fh :
                line_list = line.strip().split('\t')
                chr_list.append(str(line_list[0]))
                chr_dict[str(line_list[0])] = int(line_list[2])
    chr_list.sort()

    # mkdir for running
    subprocess.Popen(['mkdir', out_name])
    time.sleep(1)           # wait for the new folder create
    os.chdir(start_dir + '/' + out_name)


    # start move

    for chr_idx in chr_list :
        bed_ofh = open('loci_' + chr_idx + '.bed', 'w')  # loci list bed file
        chr_len = chr_dict[chr_idx]
        # if length(chr) < window, drop
        if chr_len < window_size :
            continue
        # find window in this chr/scaffold
        w_start = 1             # start 1 based
        w_end = window_size     # end 1 based, [] , convert for faidx
                                # [w_start - 1, w_end]

        while w_end <= chr_len :
            tmp_seq_list = []
            all_pass = True
            for n in range(0, len(sample_list)) :
                tmp_seq = fa_dict[sample_list[n]][chr_idx][w_start - 1 : w_end].seq
                if tmp_seq.count('N') / window_size <= n_rate :
                    # if sequence pass n_rate
                    tmp_seq_list.append(tmp_seq)
                else :
                    # if sequence fail n_rate
                    w_start += small_size
                    w_end += small_size
                    all_pass =False
                    break       # break or loop

            if all_pass :
                # pass loci filter
                # output bed file
                print('\t'.join([chr_idx, str(w_start), str(w_end)]),
                      file=bed_ofh)

                # output philip file
                subprocess.Popen(['mkdir',
                                  '_'.join([chr_idx, str(w_start), str(w_end)])])
                time.sleep(1)       # wait for new folder create
                os.chdir('./' + '_'.join([chr_idx, str(w_start), str(w_end)]))

                with open('_'.join([chr_idx, str(w_start), str(w_end)]) + '.philip', 'w') as tmp_ofh :
                    print('\t'.join([str(len(sample_list)), str(window_size)]), file=tmp_ofh)

                    for i in range(0, len(sample_list)):
                        print('\t'.join([sample_list[i], tmp_seq_list[i]]), file=tmp_ofh)
                os.chdir('../')
                w_start += step_size
                w_end += step_size
        bed_ofh.close()




def arg() :
    '''
    argument parser
    '''

    parse = argparse.ArgumentParser(prog='gen_loci_window')

    parse.add_argument('-f', '--fa',
                       help='fasta file name list, ID\sFA_FILE', required=True, type=str)
    parse.add_argument('-w', '--window', help='window size', required=True, type=int)
    parse.add_argument('-s', '--step', help='step size', required=True, type=int)
    parse.add_argument('-l', '--small', help='small step size', required=True, type=int)
    parse.add_argument('-b', '--bed', help='chr size bed file', required=False, type=str)
    parse.add_argument('--chr', help='chr_start_end, default :Fasle', default=False, type=str)
    parse.add_argument('-n', '--n_rate', help='missing rate allow for each sequence', required=True, type=float)
    parse.add_argument('-o', '--out', help='output directory name', required=True, type=str)



    args = parse.parse_args()

    gen_loci(args)


if __name__ == '__main__':
    arg()
