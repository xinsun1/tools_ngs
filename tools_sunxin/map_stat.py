#!/usr/bin/python3

__author__ = 'sunxin'

from lib_cmd import *
import argparse, os
import multiprocessing

def get_stat(bam) :
    WDIR = os.getcwd()

    b_stat = cmd_out(cmd_line=('samtools stats CMDARG1.bam'), #[short_name, type]
                     in_p=[bam],
                     work_dir=WDIR,
                     wait=True,
                     outfile='stat_' + bam,
                     shell=False)
    b_stat.run()

    run_base_c = cmd_sh(cmd_line='grep \'bases mapped (cigar):\' stat_CMDARG1',
                        in_p=[bam],
                        work_dir=WDIR,
                        shell=True)
    base_c = int(run_base_c.run().split('\t')[2])

    # bam to bed
    b2bed = cmd_out(cmd_line=('bedtools bamtobed -i CMDARG1.bam '), #[short_name, type]
                     in_p=[bam],
                     work_dir=WDIR,
                     wait=True,
                     outfile=bam + '.bed',
                     shell=False)
    b2bed.run()

    # bed merge
    b_merge = cmd_out(cmd_line=('bedtools merge -i CMDARG1.bed '), #[short_name, type]
                     in_p=[bam],
                     work_dir=WDIR,
                     wait=True,
                     outfile='merge_' + bam + '.bed',
                     shell=False)
    b_merge.run()

    # count coverage
    TIGER_WG = 2332967107

    cov_count = cmd_sh(cmd_line=('awk \'BEGIN{a=0}{a+=($3-$2)}END{print a}\' merge_CMDARG1.bed'),
                       in_p=[bam],
                       work_dir=WDIR,
                       shell=True)
    cov_count_out= cov_count.run()
    return [bam, round(float(base_c)/float(TIGER_WG), 3),
            round(float(base_c)/float(cov_count_out)),
            round(float(cov_count_out) / float(TIGER_WG), 3)]


def run(args):
    sample_list = []

    if args.file :
        fh = open(args.file, 'r')

        while 1 :
            lh = fh.readline().strip()

            if len(lh) == 0 :
                break

            sample_list.append(lh)

        fh.close()
    elif args.bam :
        sample_list.append(args.bam)

    if len(sample_list) == 0 :
        return 0

    mp_pool = multiprocessing.Pool(args.np)
    stat_list = mp_pool.map(get_stat, sample_list)

    ofh = open(args.out, 'w')
    for i in stat_list :
        print('\t'.join(str(j) for j in i) , file=ofh)

    ofh.close()


def arg() :
    '''
    argument
    '''

    parse = argparse.ArgumentParser(prog = 'map_stat')

    parse.add_argument('-p', '--np', help=' maximum number of pool workers, default= 10',
                       type=int, default=10)
    parse.add_argument('-i', '--bam', help='input bam file, without .bam', type=str, default=False)
    parse.add_argument('-o', '--out', help='out file name', required=True, type=str)
    parse.add_argument('-f', '--file', help='file contain sample list, without .bam,ONE per LINE',
                       type=str, default=False)

    args = parse.parse_args()

    run(args)


if __name__ == '__main__' :
    arg()
