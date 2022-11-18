#!/usr/bin/env python3

__author__ = 'sunxin'

'''

NOT DONE



generate a graph file for qpGraph.

based on the idea of binary tree.

include:
    tree construct
    tree copy
    add new node to tree
    add migration to tree
    tree output
'''

class Node(object) :
    def __init__(self, name=None, outleft=None, outright=None,
                 inleft=None, inright=None, mix_list=None):
        self.name=name                  # Node name
        self.outleft=outleft            # left out node
        self.outright=outright          # right out node
        self.mix_list=mix_list          # mix node list

class Tree(object) :
    def __init__(self, root=None):
        self.root = root

    def clone(self, root):
        if root:
            new_root = Node(root.name)
            new_root.inleft = root.inleft
            new_root.inright = root.inright
            new_root.mix_list = root.mix_list
        else:
            return None

        if root.outleft.inright == root :
            return None
        else :
            new_root.outleft = self.clone(root.outleft)

        if root.outright.inright == root :
            return None
        else :
            new_root.outright = self.clone(root.outright)
        return new_root


