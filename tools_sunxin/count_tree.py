#!/usr/bin/python3

__author__ = 'sunxin'

'''
This script count Tree topology for a given tree file.
Count topology per chromosome / scaffold.

'''


from ete3 import Tree
import sys


def run(tree_file, id_bed, chr_size_file, group_id, out_file_name):


    #### read chr size file ####
    chr_size_fh = open(chr_size_file, 'r')
    chr_dict = {}
    while 1 :
        chr_size_lh = chr_size_fh.readline().strip().split(' ')

        if len(chr_size_lh) != 2 :
            break

        chr_dict[str(chr_size_lh[0])] = int(chr_size_lh[1])
    chr_size_fh.close()

    #### process count ####

    tree_fh = open(tree_file, 'r')
    bed_fh = open(id_bed, 'r')
    group_list = group_id.split(',')

    count_dict = {}
    while 1 :
        bed_lh = bed_fh.readline().strip().split('\t')
        tree_lh = tree_fh.readline().strip()

        if len(bed_lh) != 3 :
            break

        t = Tree(tree_lh)
        t.prune(group_list)
        t.sort_descendants()
        t.set_outgroup(group_list[0])
        ct = str(t.write(format=9))

        chr_bed = str(bed_lh[0])
        if chr_bed in count_dict.keys() :
            if ct in count_dict[chr_bed].keys() :
                count_dict[chr_bed][ct] += 1
            else :
                count_dict[chr_bed][ct] = 1
        else :
            count_dict[chr_bed] = {}
            count_dict[chr_bed][ct] = 1

    tree_fh.close()
    bed_fh.close()

    #### output count result ####
    o_fh = open(out_file_name, 'w')

    o_head = ['chr', 'topology', 'count', 'total_count', 'chr_len', 'comp_group']
    print('\t'.join(o_head), file=o_fh)
    for i in count_dict.keys() :
        ##### calculate total count per chr #####
        chr_total = 0
        for j in count_dict[i].keys() :
            chr_total += count_dict[i][j]

        ##### output per record #####
        for k in count_dict[i].keys() :
            o_list = [i,                           # chr
                      k,                           # topology
                      str(count_dict[i][k]),       # count
                      str(chr_total),              # total count
                      str(chr_dict[i]),            # chr length
                      str(group_id)]               # group for compare
            print('\t'.join(o_list), file=o_fh)

    o_fh.close()

if __name__ == '__main__' :
    print('Count tree topology per chr/scaffold.\n'
          'Take outgroup ids, use the first one as outgroup.'
          '\n'
          'Usage: \n'
          'count_tree.py tree_file tree_ID.bed chr_size_file group_id out_file_name\n')

    if len(sys.argv) == 1 :
        print('PLEASE INPUT')
    else :
        run(sys.argv[1],       # tree file
            sys.argv[2],       # tree ID bed
            sys.argv[3],       # chr size file
            sys.argv[4],       # group ID
            sys.argv[5]        # out file name
            )



