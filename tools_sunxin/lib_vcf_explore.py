#!/usr/bin/python3
__author__ = 'sunxin'

import vcf
import os
import argparse
import subprocess
from lib_cmd import *


def get_info_vcf_line(info, line, pr_line) :
    if info in line.INFO.keys() :
        if type(line.INFO[info]) is list :
            pr_line.append(','.join([str(i) for i in line.INFO[info]]))
        else :
            pr_line.append(line.INFO[info])
    else :
        pr_line.append("NA")


def run_explore_indi(arg_in) :
    '''
    vcf explore for each child vcf
    '''

    # change to work directory
    os.chdir(arg_in['wdir'])

    FH = vcf.VCFReader(open(arg_in['child_vcf'], 'r'))

    out_file_name = 'explore_' + arg_in['child_vcf']
    OFH = open(out_file_name, 'w')

    out_head = []

    is_begin = True

    for l in FH :

        pr_line = []

        if arg_in['qual'] :
            if is_begin :
                out_head.append('QUAL')

            pr_line.append(l.QUAL)

        if arg_in['info_field'] :
            info_list = arg_in['info_field'].strip().split(',')

            if is_begin :
                out_head = out_head + info_list

            for i in info_list :
                get_info_vcf_line(i, l, pr_line)

        if arg_in['format_sample'] :
            sample_list = arg_in['format_sample'].strip().split(',')
            format_list = arg_in['format_field'].strip().split(',')

            for i in  sample_list :
                for j in format_list :

                    if is_begin :
                        out_head.append(str(i) + '_' + str(j))

                    pr_line.append(l.genotype(i)[j])

        if is_begin :
            if arg_in['head'] :
                print(*out_head, file=OFH, sep='\t')
            is_begin = False

        print(*pr_line, file=OFH, sep='\t')

    if OFH :
        OFH.close()

    return str(out_file_name)

    # if plot :
    #     plot_arg = ['Rscript', '/share/users/sunx/9-program/tools_sunxin/plot_explore.R',
    #                 str(input_vcf)]
    #     subprocess.Popen(plot_arg)

def explore_reduce(child_list, wdir, out_name, out_dir, dir_R) :
    '''
    reduce child list to final result
    '''

    os.chdir(wdir)


    reduce_child = cmd_out(cmd_line='cat CMDARG1',
                           in_p=[' '.join(child_list)],
                           work_dir=wdir,
                           outfile=out_dir + '/' + out_name,
                           shell=True,
                           wait=True)
    reduce_child.run()

    R_plot = cmd(cmd_line='Rscript CMDARG1 CMDARG2',
                 in_p=[dir_R + '/plot_explore.R', out_name],
                 work_dir=out_dir,
                 shell=False,
                 wait=True)
    R_plot.run()

    return 0

