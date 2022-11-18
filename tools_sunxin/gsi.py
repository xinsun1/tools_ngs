#!/usr/bin/python3

__author__ = 'sunxin'

'''
This script calculates genealogical sorting index (GSI) for a given population.


'''

from ete3 import Tree
import os, sys


COUNT = 0

def get_gs(start_node, pop_list) :
    global COUNT

    node_children = start_node.children
    if len(node_children) == 0 :
        return 0

    COUNT += len(node_children) - 1 # count the current node
    for i in start_node.children :
        node_indi = i.get_leaf_names()
        for j in pop_list :
            if j in node_indi :
                get_gs(i, pop_list)
                break


def get_all_nodes(tree) :
    n_node = 1
    for i in tree.iter_descendants() :
        n_node += len(i.children) - 1

    return n_node


def get_gsi(rooted_tree, pop_list):
    global COUNT

    gs_max = 1
    min_node = len(pop_list) - 1
    common_an = rooted_tree.get_common_ancestor(pop_list)

    #### calculate gs min ####
    gs_all_node_count = get_all_nodes(common_an)
    gs_min = min_node / gs_all_node_count

    #### calculate gs obs ####
    COUNT = 0
    get_gs(common_an, pop_list)
    gs_obs = min_node / COUNT

    #### calculate gsi ####
    gsi = (gs_obs - gs_min) / (gs_max - gs_min)

    return round(gsi, 3)

def run(tree_file, pop_file, out_file) :

    #### read pop file ####
    pop_fh = open(pop_file, 'r')
    pop_dict = {}
    while 1 :
        pop_lh = pop_fh.readline().strip().split(' ')

        if len(pop_lh) != 2 :
            break

        if pop_lh[1] in pop_dict.keys() :
            pop_dict[pop_lh[1]].append(pop_lh[0])
        else :
            pop_dict[pop_lh[1]] = [pop_lh[0]]
    pop_fh.close()
    pop_list = pop_dict.keys()
    pop_list = sorted(pop_list)

    outgroup_list = pop_dict['OUT']

    #### process ####
    tree_fh = open(tree_file, 'r')
    o_fh = open(out_file, 'w')

    header = []
    header_print = 0
    while 1:
        tree_lh = tree_fh.readline().strip()

        if len(str(tree_lh)) <= 1 :
            break

        tree_raw = Tree(tree_lh)
        out_list = []

        for i in range(0, len(outgroup_list)) :
            tree_raw.set_outgroup(outgroup_list[i])
            tree_rooted = tree_raw
            for j in range(0, len(pop_list)) :
                gsi_j = get_gsi(tree_rooted, pop_dict[pop_list[j]])
                out_list.append(str(gsi_j))
                if header_print != 1 :
                    header.append(str(outgroup_list[i])+ "_" + str(pop_list[j]))

        if header_print != 1 :
            print('\t'.join(header), file=o_fh)
            header_print = 1
        print('\t'.join(out_list), file=o_fh)

    tree_fh.close()
    o_fh.close()


if __name__ == '__main__' :
    print('Calculat genealogical sorting index for population.\n'
          '\n'
          'Usage: \n'
          'gsi.py tree_file pop_file out_file_name\n')

    if len(sys.argv) == 1 :
        print('PLEASE INPUT')
    else :
        run(sys.argv[1],       # tree file
            sys.argv[2],       # pop file
            sys.argv[3],       # out file name
            )
