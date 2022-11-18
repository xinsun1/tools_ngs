#!/usr/bin/python3
__author__ = 'sunxin'

'''
main control for pipeline

Last modified : 190321
'''

import multiprocessing
import os, argparse
from lib_file import *
from lib_cmd import *
from lib_split_file import *
from lib_vcf_explore import *
from lib_p_dis import *
from lib_dbSNP import *

# define environment
dir_R = '/home/sunx/0-share/9-program/tools_sunxin/R_script_sunxin'
REF_WG = '$TIGER_WG'
REF_MT = '$TIGER_MT'


def run_map(arg_in) :
    '''
    map function
    '''

    file_obj = file(name=str(arg_in['name']),
                    short_name=str(arg_in['short_name']),
                    work_dir=str(arg_in['wdir']))

    if arg_in['map']['se'] :
        if arg_in['map']['cuta']:
            file_obj.cuta_se(trim_start = arg_in['map']['trim_ad'])

        file_obj.map_se(t=int(arg_in['threads']), wg_ref=arg_in['map']['ref_wg'],
                        mt_ref=arg_in['map']['ref_mt'], mem_arg=arg_in['map']['wg_arg'])
        if arg_in['map']['stat'] :
            file_obj.map_stat(se=True)
        return 0

    if arg_in['map']['all'] :
        if arg_in['map']['single_strand'] :
            file_obj.cuta_ss(trim_start = arg_in['map']['trim_ad'])
        else :
            file_obj.cuta(trim_start = arg_in['map']['trim_ad'])

        if arg_in['map']['clp'] :
            file_obj.adar(threads=int(arg_in['threads']))
            file_obj.collapse_mt(t=int(arg_in['threads']), ref=arg_in['map']['ref_mt'])
            file_obj.collapse_wg(t=int(arg_in['threads']), ref=arg_in['map']['ref_wg'])
            file_obj.map_stat(clp=True)
        else :
            file_obj.mt_map(t=int(arg_in['threads']), ref=arg_in['map']['ref_mt'])
            file_obj.wg_map(t=int(arg_in['threads']), ref=arg_in['map']['ref_wg'],
                            mem_arg=arg_in['map']['wg_arg'])
            file_obj.map_stat()
    else :
        if arg_in['map']['cuta']:
            if arg_in['map']['single_strand'] :
                file_obj.cuta_ss(trim_start = arg_in['map']['trim_ad'])
            else :
                file_obj.cuta(trim_start = arg_in['map']['trim_ad'])

        if arg_in['map']['clp'] :
            if arg_in['map']['clp_only'] :
                file_obj.adar(threads=int(arg_in['threads']))
                return 0
            file_obj.adar(threads=int(arg_in['threads']))
            if arg_in['map']['mt'] :
                file_obj.collapse_mt(t=int(arg_in['threads']), ref=arg_in['map']['ref_mt'])
            if arg_in['map']['wg'] :
                file_obj.collapse_wg(t=int(arg_in['threads']), ref=arg_in['map']['ref_wg'])
            if arg_in['map']['stat'] :
                file_obj.map_stat(clp=True)
        else :
            if arg_in['map']['mt'] :
                file_obj.mt_map(t=int(arg_in['threads']), ref=arg_in['map']['ref_mt'])
            if arg_in['map']['wg'] :
                file_obj.wg_map(t=int(arg_in['threads']), ref=arg_in['map']['ref_wg'],
                                mem_arg=arg_in['map']['wg_arg'])
            if arg_in['map']['stat'] :
                file_obj.map_stat()


def run_snp(arg_in) :
    '''
    snp function
    '''
    pass


def run_diy(arg_in) :
    '''
    diy function
    '''
    file_obj = file(name=str(arg_in['name']),
                    short_name=str(arg_in['short_name']),
                    work_dir=str(arg_in['wdir']),
                    cwdir=False)
    file_obj.diy(work_dir=str(arg_in['diy']['work_dir']),
                 cmd_str=str(arg_in['diy']['cmd_str']),
                 in_p=str(arg_in['diy']['in_p']),
                 wait=arg_in['diy']['wait'],
                 out_file=arg_in['diy']['out_file'],
                 shell=arg_in['diy']['shell'])


def run_vcf_explore(args) :
    '''
    vcf explore function
    args:
    args from argparser
    '''

    # in file split
    in_vcf = args.vcf
    in_vcf_fh = split_file(file=in_vcf,
                           np=args.np,
                           wdir=args.wdir
                           )
    child_vcf_list = in_vcf_fh.split_head()

    # child_file mp run
    arg_list = []
    arg_indi  = {
        'child_vcf' : '',
        'info_field' : False,
        'qual' : False,
        'format_sample' : False,
        'format_field' : False,
        'wdir' : '',
        'head' : False
    }

    arg_indi['info_field'] = args.info_field
    arg_indi['qual'] = args.qual
    arg_indi['format_sample'] = args.format_sample
    arg_indi['format_field'] = args.format_field
    arg_indi['wdir'] = args.wdir + '/tmp_' + in_vcf

    write_head = 1
    for i in child_vcf_list :
        arg_indi['child_vcf'] = i
        if write_head == 1 :
            arg_indi['head'] = True
            write_head = 0
        arg_app = arg_indi.copy()
        arg_list.append(arg_app)
    mp_pool = multiprocessing.Pool(args.np)
    child_explore_list = mp_pool.map(run_explore_indi, arg_list)

    global dir_R
    # reduce result
    explore_reduce(child_list=child_explore_list,
                   wdir=arg_indi['wdir'],
                   out_dir=args.wdir,
                   out_name='explore_' + in_vcf,
                   dir_R=dir_R)
    in_vcf_fh.destory()

def run_pdis(args) :
    '''
    p-distance calculation function
    args:
    args from argparser
    '''

    # in file split
    in_vcf = args.vcf
    in_vcf_fh = split_file(file=in_vcf,
                           np=args.np,
                           wdir=args.wdir
                           )
    child_vcf_list = in_vcf_fh.split_head()


    # child_file mp run
    arg_list = []
    arg_indi  = {
        'child_vcf' : '',
        'wdir' : '',
    }

    arg_indi['wdir'] = args.wdir + '/tmp_' + in_vcf

    for i in child_vcf_list :
        arg_indi['child_vcf'] = i
        arg_app = arg_indi.copy()
        arg_list.append(arg_app)

    mp_pool = multiprocessing.Pool(args.np)
    child_pdis = mp_pool.map(pdis_indi, arg_list)

    # reduce result
    os.chdir(args.wdir)
    if len(child_pdis) == 0 :
        return 0

    od_fh = open('p_dis_array_' + in_vcf, 'w')
    op_fh = open('p_dis_pair_' + in_vcf, 'w')
    print('\t'.join(child_pdis[0][0]), file=od_fh)
    print('\t'.join(child_pdis[0][0]), file=op_fh)

    odis = child_pdis[0][1]
    opair = child_pdis[0][2]

    for i in range(1, len(child_pdis)) :
        odis = odis + child_pdis[i][1]
        opair = opair + child_pdis[i][2]

    for i in range(odis.shape[0]) :
        print('\t'.join(map(str, list(odis[i,]))), file=od_fh)
        print('\t'.join(map(str, list(opair[i,]))), file=op_fh)
    od_fh.close()
    op_fh.close()


    in_vcf_fh.destory()



def run(args) :
    '''
    process args, mp_run
    '''

    if args.subparser_name == 'pdis' :
        run_pdis(args)
    elif args.subparser_name == 'vcf_explore' :
        run_vcf_explore(args)

    arg_list = de_arg(args)

    if len(arg_list) == 0 :
        return 0

    if args.subparser_name == 'map' :
        mp_pool = multiprocessing.Pool(int(arg_list[0]['np']))
        mp_pool.map(run_map, arg_list)

    elif args.subparser_name == 'diy' :
        mp_pool = multiprocessing.Pool(int(arg_list[0]['np']))
        mp_pool.map(run_diy, arg_list)

    elif args.subparser_name == 'snp' :
        pass




def de_arg(args) :
    '''
    decode args from argparser
    :return
    arg list
    '''

    arg_list = []
    arg_indi = {
        'name' : '',
        'short_name' : '',
        'wdir' : args.wdir,
        'threads' : args.threads,
        'np' : 0,
        'run_type' : args.subparser_name,
        'map' : {'all' : True,
                 'cuta' : False,
                 'trim_ad' : False,
                 'mt' : False,
                 'wg' : False,
                 'wg_arg' : '',
                 'single_strand' : False,
                 'clp' : False,
                 'stat' : False,
                 'ref_mt' : '$TIGER_MT',
                 'ref_wg' : '$TIGER_WG',
                 'clp_only' : False,
                 'se' : False,
                 },
        'diy' : {'cmd_str' : '',
                 'work_dir' : '',
                 'in_p' : '',
                 'wait' : True,
                 'out_file' : False,
                 'shell' : False,
                 },
        'snp' : {},
        }

    if args.subparser_name == 'map':

        # read sample names and short names
        if args.file :
            sample_list = []
            s_name_list = []
            sample_file = open(args.file, 'r')
            while 1 :
                sample_file_line = sample_file.readline().strip().split(' ')
                if len(sample_file_line) == 1 :
                    break
                sample_list.append(sample_file_line[0])
                s_name_list.append(sample_file_line[1])
            sample_file.close()
        else :
            sample_list = args.indi.split(',')
            s_name_list = args.s_name.split(',')
        # read map args
        arg_indi['map']['all'] = args.all
        arg_indi['map']['cuta'] = args.cuta
        arg_indi['map']['trim_ad'] = args.trim_ad
        arg_indi['map']['mt'] = args.mt
        arg_indi['map']['wg'] = args.wg
        arg_indi['map']['wg_arg'] = args.wg_arg
        arg_indi['map']['single_strand'] = args.single_strand
        arg_indi['map']['clp'] = args.clp
        arg_indi['map']['stat'] = args.stat
        arg_indi['map']['ref_mt'] = args.ref_mt
        arg_indi['map']['ref_wg'] = args.ref_wg
        arg_indi['map']['clp_only'] = args.clp_only
        arg_indi['map']['se'] = args.se


        if len(sample_list) < args.np :
            arg_indi['np'] = len(sample_list)
        else :
            arg_indi['np'] = args.np

        for i in range(0, len(sample_list)) :
            arg_indi['name'] = sample_list[i]
            arg_indi['short_name'] = s_name_list[i]
            arg_app = arg_indi.copy()           # USE dict.copy() a=b do not work for dict
            arg_list.append(arg_app)

    elif args.subparser_name == 'diy' :
        # read sample names and short names
        if args.file :
            sample_list = []
            s_name_list = []
            sample_file = open(args.file, 'r')
            while 1 :
                sample_file_line = sample_file.readline().strip().split(' ')
                if len(sample_file_line) == 1 :
                    break
                sample_list.append(sample_file_line[0])
                s_name_list.append(sample_file_line[1])
            sample_file.close()
        else :
            sample_list = args.indi.split(',')
            s_name_list = args.s_name.split(',')
        # read snp args
        arg_indi['diy']['cmd_str'] = args.cmd_str
        arg_indi['diy']['in_p'] = args.in_p
        arg_indi['diy']['wait'] = args.wait
        arg_indi['diy']['out_file'] = args.out_file
        arg_indi['diy']['shell'] = args.shell
        arg_indi['diy']['work_dir'] = args.work_dir

        if len(sample_list) < args.np :
            arg_indi['np'] = len(sample_list)
        else :
            arg_indi['np'] = args.np

        for i in range(0, len(sample_list)) :
            arg_indi['name'] = sample_list[i]
            arg_indi['short_name'] = s_name_list[i]
            arg_app = arg_indi.copy()
            arg_list.append(arg_app)

    elif args.subparser_name == 'snp' :
        pass

    return arg_list


def arg() :
    '''
    argument management
    '''

    parse = argparse.ArgumentParser(prog='tools_sunxin_pipeline')

    # main parser

    parse.add_argument('--wdir', help='work directory (base), default=. ', type=str,
                       default=os.getcwd())
    parse.add_argument('-t', '--threads', help='threads, default= 10',
                       type=int, default=10)
    parse.add_argument('-p', '--np', help=' maximum number of pool workers, default= 10',
                       type=int, default=10)

    sub = parse.add_subparsers(help='see subcommand', dest='subparser_name')


    # map parser
    map = sub.add_parser('map', help='raw data process and map')

    map.add_argument('-i', '--indi', help='input indivdual, xxx,xxx', type=str)
    map.add_argument('-s', '--s_name', help='short name, xxx,xxx', type=str)
    map.add_argument('-f', '--file', help='input file list, default=False',
                    type=str, default=False)
    map.add_argument('--se', help='input fq is single end, default=False',
                    action='store_true')
    map.add_argument('--all', help='conduct all analysis, default=True',
                    action='store_false')
    map.add_argument('--cuta', help='cut adaptor and trim sequence, default=False',
                    action='store_true')
    map.add_argument('--trim_ad', help='trim raw read from both end, take as 3,3,'
                                       'default=False', type=str, default=False)
    map.add_argument('--single_strand', help='cut adaptor for ss library, default=False',
                    action='store_true')
    map.add_argument('--mt', help='map to mitochondria sequence, default=False',
                    action='store_true')
    map.add_argument('--wg', help='map to whole genome, default=False',
                    action='store_true')
    map.add_argument('--wg_arg', help='bwa mem extra arguments, default=''',
                    type=str, default='')
    map.add_argument('--clp', help='collapse PE reads, default=False',
                    action='store_true')
    map.add_argument('--stat', help=('calculate statistics for map results,'
                                     'Canel this option if SE read is used'
                                     ' file_name = map_summary, default=False'),
                    action='store_true')
    map.add_argument('--ref_mt', help='Map reference, default=$TIGET_MT',
                     type=str, default='$TIGER_MT')
    map.add_argument('--ref_wg', help='Map reference, default=$TIGER_WG',
                     type=str, default='$TIGER_WG')
    map.add_argument('--clp_only', help='only collapse PE reads, default=False',
                     action='store_true')




    # snp parser
    snp = sub.add_parser('snp', help='snp calling process')



    # p-dis parser
    pdis = sub.add_parser('pdis', help='p-distance calculation')
    pdis.add_argument('-i', '--vcf', help='input vcf file', required=True, type=str)



    # vcf explore parser
    vcf_explore = sub.add_parser('vcf_explore', help='vcf file field exploration')
    vcf_explore.add_argument('-i', '--vcf', help='input vcf file', required=True, type=str)
    vcf_explore.add_argument('-f', '--info_field', help='explore fields in INFO', default=False, type=str)
    vcf_explore.add_argument('-q', '--qual', help='explore QUAL, default=True', default=True,
                        action='store_false')
    vcf_explore.add_argument('-s', '--format_sample', help='sample names to explore in FORMAT, comma separated'
                        , type=str, default=False)
    vcf_explore.add_argument('-sf', '--format_field', help='explore fields in FORMAT, comma separated'
                        , type=str, default=False)


    # diy parser
    diy = sub.add_parser('diy', help='diy command')
    diy.add_argument('-i', '--indi', help='input indivdual, xxx,xxx', type=str)
    diy.add_argument('-s', '--s_name', help='short name, xxx,xxx', type=str)
    diy.add_argument('-f', '--file', help='input file list, default=False', type=str,
                       default=False)
    diy.add_argument('--cmd_str', help='command str, with CMDARG, SHORT_NAME, NAME', type=str)
    diy.add_argument('--work_dir', help='work dir, e.g. /1-map', type=str)
    diy.add_argument('--in_p', help='CMDARG1,,CMDARG2 sep=\',,\'', type=str)
    diy.add_argument('--wait', help='wait for the cmd to finish, default=True',
                    action='store_false')
    diy.add_argument('--out_file', help='out_file name, default=False', type=str,
                    default=False)
    diy.add_argument('--shell', help='shell env, default=False',
                    action='store_true')


    args =parse.parse_args()
    run(args)


if __name__ == '__main__' :
    arg()

