#!/bin/bash
__author__ = 'sunxin'


class db :
    '''
    class for vcf filter database
    '''

    def __init__(self, file_name) :
        '''

        :param file_name:
         db file name
        :return:
        '''

        self.fh = file_name
        self.dict = {}

    def gen_dict(self):
        '''
        generate self.dict
        :return:
        '''

        fh = open(self.fh, 'r')
        while 1 :
            lh = fh.readline().strip().split('\t')

            if len(lh) == 1 :
                break

            if lh[0] in self.dict.keys() :
                if int(int(lh[1])/1000) in self.dict[lh[0]].keys() :
                    self.dict[lh[0]][int(int(lh[1])/1000)].append(str(lh[1]))
                else :
                    self.dict[lh[0]][int(int(lh[1])/1000)] = [str(lh[1])]
            else :
                self.dict[lh[0]] = {}
                self.dict[lh[0]][int(int(lh[1])/1000)] = [str(lh[1])]
        fh.close()

    def in_db(self, chr, pos):
        if str(chr) in self.dict.keys() :
            if int(int(pos)/1000) in self.dict[str(chr)] :
                if (str(pos)) in self.dict[str(chr)][int(int(pos)/1000)] :
                    return 1
                else :
                    return 0
            else :
                return 0
        else:
            return 0