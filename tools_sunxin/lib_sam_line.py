__author__ = 'sunxin'

# all position 0-based

class sam_line :

    # class for sam format line

    def __init__(self, line) :
        self.line = str(line)                                       # read line as str
        self.list = self.line.split('\t')                           # read line list
        self.cigar = str(self.list[5].strip())                      # cigar str
        self.cigar_list = []                                        # cigar list
        self.seq = str(self.list[9])                                # DNA sequence str
        self.seq_len = int(len(self.seq))                           # sequence str length
        self.mpos = int(self.list[3]) - 1

        for i in range(11, len(self.list)) :
            if str(self.list[i])[0:2] == 'NM' :
                self.NM = int(str(self.list[i]).split(':')[2])      # NM number of mismatch int
            if str(self.list[i])[0:2] == 'MD' :
                self.MD = str(str(self.list[i]).split(':')[2])      # MD tag str

        self.MD_list = []                                           # MD list
        self.con_list = []                                          # consensus list
        self.snp = []                                               # snp list


    def is_head(self) :
        # determine if the sequence is sam head
        #WARNING : already move to lib_sam_file

        if self.line[0] == '@' :
            return True
        else :
            return False


    def cigar_cal(self) :
        # decode cigar string, return a list, only MID calculated
        # [[M, I , D, M], [10, 1, 2, 10]]

        pre = 0
        cigar_list1 = []
        cigar_list2 = []
        for i in range(0,len(self.cigar)) :
            for j in ['M', 'I', 'D'] :
                if j == self.cigar[i] :
                    cigar_list1.append(j)
                    cigar_list2.append(int(self.cigar[pre : i]))
                    pre = i + 1
                    break
        cigar_list = [cigar_list1, cigar_list2]
        self.cigar_list = cigar_list


    def md_decode(self) :
        # decode MD tag, no information for I, only for M , MIS, D
        # [[M, MIS, D, MIS, M], [10, A, AT, G, 10]]

        MD_list1 = []
        MD_list2 = []
        # 'reference based position'
        pre = [-1, -1] # list for Match string
        de = [-1, -1] # list for Deletion string
        for i in range(0, len(self.MD)) :
            if self.MD[i].isnumeric() :
                if pre[0] == -1 and de[0] != -1 :
                    # ^AT 123

                    pre[0] = i
                    MD_list1.append('D')
                    if de[1] == -1 :
                        MD_list2.append(self.MD[de[0]])
                    else :
                        MD_list2.append(self.MD[de[0] : de[1]])
                    de = [-1 , -1]

                elif pre[0] == -1 and de[0] == -1 and i == 0 :
                    # start
                    pre[0] = i

                elif pre[0] == -1 and de[0] == -1 and i != 0 :
                    # T 123
                    pre[0] = i

                else :
                    # 111  123
                    pre[1] = i + 1
            else :
                if pre[0] != -1 :
                    # 123 A/^
                    MD_list1.append('M')
                    if pre[1] == -1 :
                        MD_list2.append(self.MD[pre[0]])
                    else :
                        MD_list2.append(self.MD[pre[0] : pre[1]])

                    pre = [-1, -1]

                    if self.MD[i] == '^' :
                        # 123 ^
                        de[0] = i + 1
                    else :
                        # 123 A
                        if de[0] == -1 :
                            MD_list1.append('MIS')
                            MD_list2.append(self.MD[i])
                        else :
                            # 123 ^AG
                            if i != de[0] :
                                de[1] = i + 1
                else :
                    if self.MD[i] == '^' :
                        # A ^
                        de[0] = i + 1
                    else :
                        if de[0] == -1 :
                            # A A
                            MD_list1.append('MIS')
                            MD_list2.append(self.MD[i])
                        else :
                            # ^ AA
                            if i != de[0] :
                                de[1] = i + 1
            if i == len(self.MD) - 1 :
                MD_list1.append('M')
                if pre[1] == -1 :
                    MD_list2.append(self.MD[pre[0]])
                else :
                    MD_list2.append(self.MD[pre[0] : pre[1]])


        self.MD_list = [MD_list1, MD_list2]


    def gen_con(self):
        # generate consensus build information for this seq , include M, I, D
        # 0 based
        # [[10,30],[30,33],[I, 32.1, 32.2, 32.3],[33,35]]
        # ['AGGAGAGA...', '---', 'AAA', 'AAA...']
        # selfbased position cigar_list2

        con_list1 = []
        con_list2 = []

        if len(self.cigar_list) == 0 :
            self.cigar_cal()

        pos = 0
        mpos = 0 + int(self.mpos)
        for i in range(0, len(self.cigar_list[0])) :
            if self.cigar_list[0][i] == 'M' :
                con_list1.append([mpos, mpos + self.cigar_list[1][i]])
                con_list2.append(self.seq[pos : pos + self.cigar_list[1][i]])
                mpos += self.cigar_list[1][i]
                pos += self.cigar_list[1][i]
            elif self.cigar_list[0][i] == 'D' :
                con_list1.append([mpos, mpos + self.cigar_list[1][i]])
                con_list2.append('-' *  self.cigar_list[1][i])
                mpos += self.cigar_list[1][i]
            else :
                insert_list = ['I']
                for j in range(0, self.cigar_list[1][i]) :
                    insert_list.append(str(mpos - 1) + '.' + str(j + 1))
                con_list1.append(insert_list)
                con_list2.append(self.seq[pos : pos + self.cigar_list[1][i]])
                pos += self.cigar_list[1][i]

        self.con_list = [con_list1, con_list2]


    def getbase_r(self, rpos) :
        # get base by reference based position (0-based)
        # include deletion, insertion
        # rpos float

        base = ''
        if float(rpos) != float(int(float(rpos))) :
            return ''


        if len(self.con_list) == 0 :
            self.gen_con()

        for i in range(0, len(self.con_list[0])) :
            if base != '' :
                break
            if self.con_list[0][i][0] == 'I' :
                for j in range(1, len(self.con_list[0][i])) :
                    if float(self.con_list[0][i][j]) == float(rpos) :
                        base = self.con_list[1][i][j - 1].upper()
                        break

            elif float(rpos) >= float(self.con_list[0][i][0]) and float(rpos) < float(self.con_list[0][i][1]) :
                base = self.con_list[1][i][int(rpos) - int(self.con_list[0][i][0])].upper()
                break

        # return base in the read
        return base


    def getbase_md(self, mdpos) :
        # get base from md pos (0-based)
        # md pos include M, MIS, D not include I
        if len(self.cigar_list) == 0 :
            self.cigar_cal()
        if len(self.MD_list) == 0 :
            self.md_decode()

        pos = -1
        spos = -1
        for i in range(0, len(self.cigar_list[0])) :
            if self.cigar_list[0][i] == 'M' :
                pos += self.cigar_list[1][i]
                spos += self.cigar_list[1][i]
            elif self.cigar_list[0][i] == 'D' :
                pos += self.cigar_list[1][i]
            else :
                spos += self.cigar_list[1][i]

            if pos >= mdpos :
                base = self.seq[spos - pos + mdpos].upper()
                break

        # return  [read relative pos, base in the read]
        return [spos - pos + mdpos, base]


    def gen_snp(self) :
        # generate snp information, mismatch
        # self.snp = [read relative pos, base in the reference, base in the read]

        if len(self.cigar_list) == 0 :
            self.cigar_cal()
        self.md_decode()
        mdpos = -1
        for i in range(0, len(self.MD_list[0])) :
            if self.MD_list[0][i] == 'MIS' :
                mdpos += 1
                gmd = self.getbase_md(mdpos)
                # store snp information [read pos, ref base, base in the read]

                self.snp.append([gmd[0], self.MD_list[1][i], gmd[1]])
            elif self.MD_list[0][i] == 'M' :
                mdpos += int(self.MD_list[1][i])
            else :
                mdpos += len(self.MD_list[1][i])


    def map_tg(self, pos) :
        # add mapping position tag to the sam file
        # nuclear or mitochondria
        # YM stand for whether the seq is mapping to mitochondria: 1:Yes 0:No

        if str(pos) == 'MT' :
            tag = 'YM:i:1'
        elif str(pos) == 'WG' :
            tag = 'YM:i:0'

