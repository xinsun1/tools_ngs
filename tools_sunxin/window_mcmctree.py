#!/usr/bin/python3

__author__ = 'sunxin'

'''
This script run mcmctree for each windowed based loci.

1) Prepare seq file and tree file
2) Modify the mcmctree control file.
3) Batch run mcmctree
4) Summarize the mcmctree softbond dating result.
5) Generate result for R processing.


File strcuture:

----|----START_DIR
----|----w100s50
        |----scaffold_X_X.philip
        |----RAxML_bestTree.scaffold_X_X
----|----mcmctree_w100s50
        |----mcmctree_all.ctl
        |----scaffold_X_X
            |----scaffold_X_X.philip_mod
            |----scaffold_X_X.tree_mod_no_brlen
            |----mcmctree_scaffold_X_X.ctl
        |---- ...
'''

import sys, os
import argparse
import multiprocessing
from lib_split_file import *
from lib_cmd import *
from ete3 import Tree

def read_nexus(nexus_file) :
    '''
    read in mcmctree generated nexus tree file, return newick tree file
    '''

    nexus_fh = open(nexus_file, 'r')

    n = 0
    while 1 :
        if n!= 4 :
            nexus_lh = nexus_fh.readline().strip()
            n+=1
        else :
            nexus_lh = str(nexus_lh)
            break
    nexus_fh.close()

    out_newick = ''
    is_head = 0         # mark the start of tree
    is_square = 0       # mark the start of square
    for i in nexus_lh[0:30] :
        if i == '(' and is_head == 0 :
            # find the start of the tree
            out_newick += i
            is_head = 1
            continue
        if is_head == 1 :
            # already have a start
            if i == '[' :
                # start of square
                is_square = 1
            elif i == ']' :
                is_square = 0
                is_square_end = 1
            elif is_square == 0 :
                # not in the square
                out_newick += i
    for i in nexus_lh[30:] :
        # already have a start
        if i == '[' :
            # start of square
            is_square = 1
        elif i == ']' :
            is_square = 0
            is_square_end = 1
        elif is_square == 0 :
            # not in the square
            out_newick += i

    return out_newick


def get_most_out_pop_phylo(pop_phylo_str, group_outgroup, group_inner):
    '''
    get most out branch from population phylogeny,
    return the smaller branch with all branch names
    prune the phylogeny with selected population in the beginning.
    '''

    # read population phylogeny
    pop_phylo = Tree(pop_phylo_str)
    pop_phylo.prune([group_outgroup] + group_inner) # prune tree
    pop_phylo.sort_descendants()                  # sort tree
    pop_phylo.set_outgroup(group_outgroup)        # set outgroup

    node_com_an = pop_phylo.get_common_ancestor(group_inner).get_children()
    if len(node_com_an[0].get_leaf_names()) >= \
            len(node_com_an[1].get_leaf_names()) :
        return  node_com_an[1].get_leaf_names()
    else :
        return  node_com_an[0].get_leaf_names()


def run_sum(bed_file, pop_phylo_file, mcmc_dir, group_compare,
            indi_compare, sum_result_file) :
    '''
    summarize mcmctree result.
    read nexus tree file,
    summarize by most out branch of population phylogeny
    '''

    # read bed file and population phylogeny
    bed_fh = open(bed_file, 'r')
    pop_phylo_fh = open(pop_phylo_file, 'r')

    dict_pop_phylo = {}
    while 1 :
        bed_lh = bed_fh.readline().strip().split('\t')
        pop_phylo_lh = pop_phylo_fh.readline().strip()

        if len(bed_lh) != 3 :
            break

        dict_pop_phylo['_'.join(bed_lh)] = str(pop_phylo_lh)

    bed_fh.close()
    pop_phylo_fh.close()

    # read groups for compare
    group_outgroup = group_compare.split(',')[0]   # outgroup, str
    group_inner = group_compare.split(',')[1:]     # except outgroup, list

    # read outgroup sample ID
    sample_outgroup1 = indi_compare.split(',')[0] # outgroup sample 1
    sample_outgroup2 = indi_compare.split(',')[1] # outgroup sample 2
    sample_inner = indi_compare.split(',')[2]     # inner sample

    # open summary result file for write
    sum_fh = open(sum_result_file, 'w')
    print('\t'.join(['Pos', 'First_div_pop', 'TMRCA']),file=sum_fh)

    # move to mcmc directory
    os.chdir(mcmc_dir)
    mcmc_dir_list = os.listdir()

    for i in mcmc_dir_list :
        if i[0:2] != 'sc' :
            continue
        os.chdir(mcmc_dir + '/' + i)

        if 'FigTree.tre' not in os.listdir() :
            continue

        # get common ancestor branch length
        tree_mcmc = Tree(read_nexus('FigTree.tre'))
        tree_mcmc.set_outgroup(tree_mcmc&sample_outgroup1)
        #### look for common ancestor node
        com_an_node = ''
        for k in tree_mcmc.get_children() :
            if len(k.get_leaf_names()) != 1 :
                com_an_node = k
                for j in range(0, len(k.get_children())) :
                    if k.get_children()[j].get_leaf_names() == [sample_outgroup2] :
                        com_an_node = k.get_children()[1 - j]

        len_com_an = tree_mcmc.get_distance(com_an_node, sample_inner)

        # get most out branch
        pop_most_out = get_most_out_pop_phylo(dict_pop_phylo[i],
                                              group_outgroup,
                                              group_inner)

        print('\t'.join([str(i),
                         ','.join(pop_most_out),
                         str(round(len_com_an,2))]),
              file=sum_fh)

    sum_fh.close()


def run_indi(arg_indi) :
    '''
    read Bed file as a list
    '''

    bed_file = arg_indi['child_bed']
    wdir = arg_indi['wdir']
    raxml_dir = arg_indi['raxml_dir']
    mcmctree_dir = arg_indi['mcmctree_dir']
    outgroup_id = arg_indi['outgroup_id']

    os.chdir(wdir + '/' + mcmctree_dir)
    b_fh = open(bed_file, 'r')

    while 1 :
        l_bfh = b_fh.readline().strip().split('\t')

        if len(l_bfh) == 1 :
            break

        bed_name = str('_'.join(l_bfh))

        # change dir
        os.chdir(wdir + '/' + mcmctree_dir + '/' + bed_name)

        # # mkdir
        # # create new directory for each loci
        # run_mkdir = cmd(cmd_line='mkdir -p ./CMDARG1',
        #                 in_p=[bed_name],
        #                 work_dir=wdir + '/' + mcmctree_dir,
        #                 shell=False,
        #                 wait=True,
        #                 )
        # run_mkdir.run()

        # prepare seq file
        # change tab to space for mcmctree
        run_mod_seq_file = cmd_out(
            cmd_line='awk \'{print $1, $2}\' OFS=\'   \' CMDARG1/CMDARG2.philip',
            in_p=[wdir + '/' + mcmctree_dir + '/' + bed_name, bed_name],
            work_dir=wdir + '/' + mcmctree_dir + '/' + bed_name,
            shell=False,
            wait=True,
            outfile=bed_name + '.philip_mod',
        )
        run_mod_seq_file.run()

        # prepare tree file
        # remove phylgenetic tree branch length, add number_of_species and seq_part
        # default seq_part=1, change the script for seq_part
        # reroot the tree
        # add calibration point
        from ete3 import Tree
        loci_tree_raxml_file = wdir + '/' + mcmctree_dir + '/' + bed_name + '/' +  'RAxML_bestTree.' + bed_name
        loci_tree_raxml_fh = open(loci_tree_raxml_file, 'r')

        loci_tree = Tree(loci_tree_raxml_fh.readline().strip())
        loci_tree_raxml_fh.close()

        sample_list = []

        for i in loci_tree.iter_leaves() :
            sample_list.append(i.name)

        loci_tree_file_fh = open(wdir + '/' + mcmctree_dir + '/' + bed_name +
                                 '/' + bed_name + '.tree_mod_no_brlen', 'w')
        print(str(len(sample_list)) + ' ' + '1', file=loci_tree_file_fh)

        loci_tree.set_outgroup(loci_tree&str(outgroup_id))
            # set PUN divergence time
        sample_list.remove(outgroup_id)
        loci_tree.get_common_ancestor(sample_list).name = '\'>18.2<46.2\''
        out_loci_tree = loci_tree.write(format=8)
        out_loci_tree = out_loci_tree.replace('NoName', '')
        print(out_loci_tree, file=loci_tree_file_fh)
        loci_tree_file_fh.close()

        # modify mcmctree control file
        # change mcmctree_all.ctl for each loci
        run_mod_control_file = cmd_out(
            cmd_line='sed \'s/WINDOW_POS/CMDARG1/g\' ../mcmctree_all.ctl',
            in_p=[bed_name],
            work_dir=wdir + '/' + mcmctree_dir + '/' + bed_name,
            shell=True,
            wait=True,
            outfile='mcmctree_' + bed_name + '.ctl',
            )
        run_mod_control_file.run()

        # # run mcmctree for each loci in separate sub-directory
        # run_mcmctree = cmd_out(
        #     cmd_line='mcmctree mcmctree_CMDARG1.ctl',
        #     in_p=[bed_name],
        #     work_dir=wdir + '/' + mcmctree_dir + '/' + bed_name,
        #     shell=False,
        #     wait=True,
        #     outfile='log_' + bed_name,
        # )
        # run_mcmctree.run()


def mcmc_run(np, bed_file, raxml_dir, mcmctree_dir, outgroup_id):

    # read current dir
    # start program directory
    start_dir = os.getcwd()

    # split bed file
    child_bed =  split_file(file=bed_file,
                            np=np,
                            wdir=start_dir)
    child_bed_list = child_bed.split()

    # child file mp run
    # generate argument dict for each child process
    # start child process run
    arg_list = []
    arg_indi = {
        'child_bed'     : '',
        'wdir'          : start_dir,
        'raxml_dir'     : raxml_dir,
        'mcmctree_dir'  : mcmctree_dir,
        'outgroup_id'   : outgroup_id,
    }

    for i in child_bed_list :
        arg_indi['child_bed'] = start_dir + '/tmp_' + bed_file + '/' + i
        arg_app = arg_indi.copy()
        arg_list.append(arg_app)

    mp_pool = multiprocessing.Pool(int(np))
    mp_pool.map(run_indi, arg_list)

    # reduce result
    # remove whole child bed file directory
    child_bed.destory()

def de_arg(args) :

    if args.subparser_name == 'mcmctree' :
        mcmc_run(np=args.np,
                 bed_file=args.bed,
                 raxml_dir=args.dir_raxml,
                 mcmctree_dir=args.dir_mcmc,
                 outgroup_id=args.out_ID
        )
    elif args.subparser_name == 'sum' :
        run_sum(bed_file=args.bed,
                pop_phylo_file=args.pop_phylo,
                mcmc_dir=args.dir_mcmc,
                group_compare=args.group_comp,
                indi_compare=args.indi_comp,
                sum_result_file=args.sum_file
        )


def arg() :
    '''
    argument management
    '''

    parse = argparse.ArgumentParser(prog='window_mcmctree')

    # main parser
    #parse.add_argument('--wdir', help='work directory (base), default=. ', type=str,
    #                   default=os.getcwd())

    sub = parse.add_subparsers(help='see subcommand', dest='subparser_name')

    # mcmctree parser
    mcmctree = sub.add_parser('mcmctree', help='run mcmctree for each window '\
                                               'based loci')

    mcmctree.add_argument('-n', '--np', help='number of pool workers', type=int)
    mcmctree.add_argument('-b', '--bed', help='window bed file', type=str)
    mcmctree.add_argument('--dir_raxml', help='Raxml run directory,'\
                                              'FULL_PATH(/)', type=str)
    mcmctree.add_argument('--dir_mcmc', help='mcmctree run directory,'\
                                             'SHORT_PATH(xxx)', type=str)
    mcmctree.add_argument('--out_ID', help='outgroup sample ID', type=str)

    # summary parser
    sum = sub.add_parser('sum', help='summary mcmctree result')

    sum.add_argument('-b', '--bed', help='window bed file', type=str)
    sum.add_argument('--pop_phylo', help='population phylogeny file', type=str)
    sum.add_argument('--dir_mcmc', help='mcmctree run directory,'\
                                             'FULL_PATH(\.\..)', type=str)
    sum.add_argument('--group_comp', help='group str include populations for '\
                     'compare, with the 1st POP as outgroup. etc, OUT,POP1,POP2',
                     type=str)
    sum.add_argument('--indi_comp', help='individual str include 2 OUT and '\
                     '1 POP. For branch length calculation. etc. OUT1,OUT2,POP',
                     type=str)
    sum.add_argument('--sum_file', help='summary result file name', type=str)

    args = parse.parse_args()
    de_arg(args)

if __name__ == '__main__' :
    arg()


