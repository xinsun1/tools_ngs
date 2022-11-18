__author__ = 'sunxin'

from lib_cmd import *
import subprocess, os, time

class split_file() :
    '''
    class for file object that need to be split
    '''

    def __init__(self, file, np, wdir, chunksize=1) :
        '''
        file
            in putfile name

        np
            number of processes
        '''

        self.file = file
        self.np = int(np)
        self.wdir = wdir
        self.chunk = int(chunksize)
        self.out = []

    def getlen_wc(self, file_name, wdir) :
        '''
        calculate the total line number of the file by 'wc -l'
        '''

        get_len = cmd_sh(cmd_line='wc -l CMDARG1',
                         in_p=[file_name],
                         work_dir=wdir,
                         shell=True)
        file_len = int(get_len.run().split(' ')[0])
        return file_len

    def split(self) :
        '''
        split file for tmp files
        '''

        os.chdir(self.wdir)
        subprocess.Popen(['mkdir', 'tmp_' + self.file])
        time.sleep(1)

        f_size = self.getlen_wc(file_name=self.file,
                                wdir=self.wdir)
        tmp_sp_size = f_size / self.chunk / self.np
        if float(tmp_sp_size) == float(int(tmp_sp_size)) :
            sp_size = int(tmp_sp_size)
        else :
            sp_size = int(tmp_sp_size) + 1

        self.reset()
        for i in range(self.np) :
            if i < 10 :
                self.out.append(str('tmp_0' + str(i) + '_' + self.file))
            else :
                self.out.append(str('tmp_' + str(i) + '_' + self.file))

        fh = open(self.file, 'r')

        os.chdir(self.wdir + '/tmp_' + self.file)

        ln = 0
        while 1 :
            a = fh.readline()
            if len(a) == 0 :
                ofh.close()
                break
            ln += 1

            # file control part
            if int(int(ln) % (self.chunk * sp_size)) == 1 :
                index_file = int(int(ln) / (self.chunk * sp_size))
                if index_file != 0 :
                    ofh.close()
                ofh = open(self.out[index_file], 'w')

            print(a, file=ofh, end='')
        fh.close()
        return self.out

    def reset(self) :
        self.out = []

    def destory(self) :
        '''
        remove whole sub directory
        '''

        rm_dir = cmd(cmd_line='rm -fr CMDARG1',
                     in_p=['tmp_' + self.file],
                     work_dir=self.wdir,
                     shell=False,
                     wait=False)
        rm_dir.run()

    def split_head(self) :
        '''
        split file with head
        '''

        # write head
        # head start with '#'

        self.reset()
        os.chdir(self.wdir)
        subprocess.Popen(['mkdir', 'tmp_' + self.file])
        time.sleep(1)

        get_head = cmd_out(cmd_line='grep ^# CMDARG1',
                           in_p=[self.file],
                           work_dir=self.wdir,
                           shell=True,
                           wait=True,
                           outfile='./tmp_' + self.file + '/tmp_head_' + self.file)
        get_head.run()

        # write body
        get_body = cmd_out(cmd_line='grep -v ^# CMDARG1',
                           in_p=[self.file],
                           work_dir=self.wdir,
                           shell=True,
                           wait=True,
                           outfile='./tmp_' + self.file + '/tmp_body_' + self.file)
        get_body.run()

        # change to subdir
        os.chdir(self.wdir + '/tmp_' + self.file)
        body_len  = self.getlen_wc(file_name='tmp_body_' + self.file,
                                   wdir=self.wdir + '/tmp_' + self.file)

        child_file_len = (int(body_len / self.chunk / self.np) + 1) * self.chunk

        gen_child = cmd(cmd_line='split -l CMDARG1 -d CMDARG2 CMDARG3',
                        in_p=[child_file_len, 'tmp_body_' + self.file, 'tmp_' + self.file + '_'],
                        work_dir=self.wdir + '/tmp_' + self.file,
                        shell=False)
        gen_child.run()

        num_file = body_len / child_file_len
        if float(num_file) == float(int(num_file)) :
            num_file = int(num_file)
        else :
            num_file = int(num_file) + 1
        for i in range(num_file) :
            if i < 10 :
                child_file_head = cmd_out(cmd_line='cat CMDARG1 CMDARG2',
                                          in_p=['tmp_head_' + self.file,
                                                'tmp_' + self.file + '_0' + str(i)],
                                          work_dir=self.wdir + '/tmp_' + self.file,
                                          outfile='tmp_0' + str(i) + '_' + self.file,
                                          shell=True,
                                          wait=True)
                child_file_head.run()
                self.out.append('tmp_0' + str(i) + '_' + self.file)
                rm_tmp = cmd(cmd_line='rm CMDARG1',
                             in_p=['tmp_' + self.file + '_0' + str(i)],
                             work_dir=self.wdir + '/tmp_' + self.file,
                             shell=False,
                             wait=False)
                rm_tmp.run()
            else :
                child_file_head = cmd_out(cmd_line='cat CMDARG1 CMDARG2',
                                          in_p=['tmp_head_' + self.file,
                                                'tmp_' + self.file + '_' + str(i)],
                                          work_dir=self.wdir + '/tmp_' + self.file,
                                          outfile='tmp_' + str(i) + '_' + self.file,
                                          shell=True,
                                          wait=True)
                child_file_head.run()
                self.out.append('tmp_' + str(i) + '_' + self.file)
                rm_tmp = cmd(cmd_line='rm CMDARG1',
                             in_p=['tmp_' + self.file + '_' + str(i)],
                             work_dir=self.wdir + '/tmp_' + self.file,
                             shell=False,
                             wait=False)
                rm_tmp.run()

        return self.out















































