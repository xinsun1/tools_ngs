__author__ = 'sunxin'


class fast_line :

    # class for fasta/fastq line
    def __init__(self, name, line, pair=False, qual=False) :
        # pair info get from the name line[1]

        self.name = name
        self.line = line
        self.pair_line = pair
        if pair :
            self.pair = pair.split(':')[0]
        if qual :
            self.qual = qual

    def print(self, ofh) :

        if self.qual :
            print(self.name + ' ' + self.pair_line, file=ofh)
            print(self.line, file=ofh)
            print('+', file=ofh)
            print(self.qual, file=ofh)
        else :
            print(self.name, file=ofh)
            print(self.line, file=ofh)




class fasta :

    # class for fasta file

    def __init__(self, file, is_fastq=False, is_pair=False) :

        self.fh = open(file, 'r')              # file obj
        self.is_fastq = is_fastq    # if the file is fastq file
        self.is_pair = is_pair

    def next(self) :

        name_line = self.fh.readline().strip().split(' ')

        if len(name_line) == 1 :            # file ended
            return 0

        name = name_line[0]
        pair = False
        if self.is_pair :
            pair = name_line[1]
        line = self.fh.readline().strip()
        qual = False

        if self.is_fastq :
            self.fh.readline()
            qual = self.fh.readline().strip()

        fline = fast_line(name = name, line = line, pair = pair, qual = qual)

        return fline

    def reset(self) :

        # return to the topline of the file

        self.fh.seek(0)

