__author__ = 'sunxin'


# all position 0-based

from lib_sam_line import *


class sam_file :

    # class for a sam format file

    def __init__(self, sam_fh, filter_file = False) :

        self.fh = sam_fh                       # file obj
        if filter_file :
            self.filter_fh = filter_file     # filter database file obj
        else :
            self.filter_fh = False
        self.con = {}                          # dict to store consensus STR
        self.pos = []                          # position index for consensus INT
        self.freq = {}                         # frequency dict
        self.filter_db = {}                    # filter database


    def gen_con(self, filter_db = False, filter_ancinet = False) :
        # generate consensus base on the file

        self.fh.seek(0)
        while 1 :
            line = self.fh.readline().strip()

            if len(line) == 0 :
                break

            if self.is_head(line) :
                continue

            sline = sam_line(line)

            if filter_ancinet :
                # filter for this read
                # with ancient signal
                if not self.filt_mut(sline, C2T=True, G2A=True, end=5) :
                    continue
            if filter_db :
                # with mut in filter database
                if self.filt_db(sline) :
                    continue

            self.add2con(sline)


    def add2con(self, sline):
        # add sequence to consensus
        # sam format seq

        if len(sline.con_list) == 0 :
            sline.gen_con()

        poslist = sline.con_list[0]
        seqlist = sline.con_list[1]
        for i in range(0, len(poslist)) :
            if poslist[i][0] != 'I' :
                self.add2con_seg(poslist[i], seqlist[i])
            else :
                self.add2con_seg(poslist[i], seqlist[i], insert=True)


    def add2con_seg(self, pos, seq, insert = False):
        # add sequence segment to consensus
        # pos = [s, e] 0-based, seq = str

        pos = pos
        seq = seq
        if insert :
            # if insertion
            for i in range(1,len(pos)) :
                if pos[i] in self.pos :
                    if str(seq[i - 1]) in self.con[pos[i]].keys() :
                        self.con[pos[i]][str(seq[i - 1])] += 1
                    else :
                        self.con[pos[i]][str(seq[i - 1])] = 1
                else :
                    if i == 1 :
                        # first insertion pos not in self.pos
                        insert_pos = self.pos.index(int(float(pos[i]) - 0.1))
                        new_pos = self.pos[0 : insert_pos + 1]
                        for m in pos[1:] :
                            new_pos.append(m)
                            self.con[m] = {}
                        for m in self.pos[(insert_pos + 1) :] :
                            new_pos.append(m)
                        self.pos = new_pos
                    else :
                        insert_pos = self.pos.index(str(pos[i - 1]))
                        new_pos = self.pos[0 : insert_pos + 1]
                        for m in pos[i :] :
                            new_pos.append(m)
                            self.con[m] = {}
                        for m in self.pos[(insert_pos + 1) :] :
                            new_pos.append(m)
                        self.pos = new_pos

                    if str(seq[i - 1]) in self.con[pos[i]].keys() :
                        self.con[pos[i]][str(seq[i - 1])] += 1
                    else :
                        self.con[pos[i]][str(seq[i - 1])] = 1

        else :
            # not insertion
            if len(self.pos) == 0 :
                for i in range(0, pos[1]) :
                    self.pos.append(i)
                    self.con[i] = {}
            elif float(pos[1]) > float(self.pos[-1]) + 1:
                # take 419.5 for 419
                aps = int(float(self.pos[-1])) + 1
                for i in range(aps, pos[1]) :
                    self.pos.append(i)
                    self.con[i] = {}

            # append seq to consensus
            for i in range(pos[0], pos[1]) :
                base = seq[i - pos[0]]
                if base in self.con[i].keys() :
                    self.con[i][base] += 1
                else :
                    self.con[i][base] = 1


    def cal_freq(self) :
        # calculate allele frequency base on the consensus

        for i in self.pos :
            if len(self.con[i].keys()) == 0 :
                # void position
                self.freq[i] = {}
            else :
                total = 0
                self.freq[i] = {}
                for j in self.con[i].keys() :
                    total += self.con[i][j]
                for j in self.con[i].keys() :
                    self.freq[i][j] = round(float(self.con[i][j] / total), 3)        # round the float value


    def built_filt_db(self) :
        # read filter file when the program start
        # build database based file for read
        # filter file format : pos \t base,base,base

        if len(self.filter_db.keys()) == 0 :
            if self.filter_fh :
                # read in filter file
                fdb_fh = open(self.filter_fh, 'r')
                while 1 :
                    ffl = fdb_fh.readline().strip()

                    if len(ffl) == 0 :
                        break

                    if ffl[0] == '#' :
                        continue
                    ffl = ffl.split('\t')

                    self.filter_db[str(int(ffl[0]) - 1)] = ffl[1].split(',')


    def build_filt_freq(self, low_freq = 0.3) :
        # filter for low frequency
        if len(self.freq.keys()) == 0 :
            self.cal_freq()
        for i in self.pos :
            if len(self.freq[i].keys()) == 0 :
                continue
            for j in self.freq[i].keys() :
                if self.freq[i][j] < low_freq :
                    if i in self.filter_db.keys() :
                        if not j in self.filter_db[i] :
                            self.filter_db[i].append(j)
                    else :
                        self.filter_db[i] = [j]


    def filt_mut(self, sline, C2T = False, G2A = False, end = 5) :
        # mutation type based filter for read
        # currently filter for c2t or g2a in both end(int 5)

        with_filt = False
        if len(sline.snp) == 0 :
            sline.snp()

        for i in sline.snp :
            if (i[0] >= 0 and i[0] <= end - 1) or (i[0] <= sline.seq_len - 1 and i[0] >= sline.seq_len - end) :
                if C2T :
                    if i[1] == 'C' and i[2] == 'T' :
                        with_filt = True
                        break
                if G2A :
                    if i[1] == 'G' and i[2] == 'A' :
                        with_filt = True
                        break

        # return if the read has the filter option
        return with_filt


    def filt_db(self, sline) :
        # filter db based filter for the if the read has N in pos

        with_filt = False
        seqs = float(sline.mpos)
        seqe = float(sline.mpos + sline.seq_len)
        for i in self.filter_db.keys() :
            if with_filt :
                break
            if float(i) >= seqs and float(i) <= seqe :
                for j in self.filter_db[i] :
                    read_base = str(sline.getbase_r(i))
                    if str(j) == read_base :
                        with_filt = True
                        break

        return with_filt


    def iter_con(self, iter_time = 3 , filter_ancient = False, low_freq = 0.3) :
        # iteration generation consensus
        ## build filter from db file

        self.built_filt_db()
        # print(self.filter_db)

        self.con = {}
        if filter_ancient :
            self.gen_con(filter_ancient = True, filter_db = True)
        else :
            self.gen_con(filter_db = True)

        self.freq = {}
        self.cal_freq()
        self.build_filt_freq(low_freq = low_freq)


        for i in range(0, int(iter_time)) :
            # generate consensus from fresh
            self.con = {}
            for j in self.pos :
                self.con[j] = {}

            if filter_ancient :
                self.gen_con(filter_ancient = True, filter_db= True)
            else :
                self.gen_con(filter_db = True)

            # calculate alelle freq from fresh
            self.freq = {}
            self.cal_freq()
            self.build_filt_freq(low_freq = low_freq)


    def conseq(self, iupac = False, filter_con = False) :
        # generate consensus sequence base on the current con dict
        # return a seq str

        if len(self.pos) == 0 :
            self.gen_con()

        con = ''
        ambi = 0
        for i in self.pos :
            if len(self.con[i].keys()) == 0 :
                if float(i) == int(float(i)) :
                    con += 'N'
            elif len(self.con[i].keys()) == 1 :
                con += '/'.join(self.con[i].keys())
            else :
                con += '(' + '/'.join(self.con[i].keys()) + ')'
                ambi += 1

        return [con, ambi]


    def out_ambi(self) :
        # print ambiguous site allele frequency

        for i in self.pos :
            if len(self.con[i].keys()) > 1 :
                print(i ,self.freq[i])


    def is_head(self,line) :
        # determine if the sequence is sam head

        if line[0] == '@' :
            return True
        else :
            return False


    def snp_stat(self, read_length = 150) :
        # generate snp stats for this sam file
        # position

        self.snpstat_L = {}
        self.snpstat_R = {}

        for i in range(0,150) :
            self.snpstat_L[str(i)] = {'AT' : 0, 'AC' : 0, 'AG' : 0, 'AN' : 0,
                                      'TA' : 0, 'TC' : 0, 'TG' : 0, 'TN' : 0,
                                      'CA' : 0, 'CT' : 0, 'CG' : 0, 'CN' : 0,
                                      'GA' : 0, 'GT' : 0, 'GC' : 0, 'GN' : 0}
            self.snpstat_R[str(i)] = {'AT' : 0, 'AC' : 0, 'AG' : 0, 'AN' : 0,
                                      'TA' : 0, 'TC' : 0, 'TG' : 0, 'TN' : 0,
                                      'CA' : 0, 'CT' : 0, 'CG' : 0, 'CN' : 0,
                                      'GA' : 0, 'GT' : 0, 'GC' : 0, 'GN' : 0}
        self.fh.seek(0)

        file_length = 0
        while 1 :
            line = self.fh.readline().strip()

            if len(line) == 0 :
                break

            if self.is_head(line) :
                continue

            sline = sam_line(line)
            sline_length = sline.seq_len
            file_length += 1

            sline.gen_snp()

            for i in sline.snp :
                self.snpstat_L[str(i[0])][str(i[1]) + str(i[2])] += 1
                self.snpstat_R[str(sline_length - int(i[0]) - 1)][str(i[1]) + str(i[2])] += 1
        return(self.snpstat_L, self.snpstat_R, file_length)





















