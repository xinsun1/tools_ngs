__author__ = 'sunxin'



class cmd :
    def __init__(self, cmd_line, in_p, work_dir, shell=False, wait=True) :
        '''

        :param cmd_line:
         ful cmd line
        :param in_p:
         input, list of str consider as input
        :param work_dir:
         cmd working dir
        :param wait:
         if the cmd need wait
        :return:
         re
        '''

        self.cmd = str(cmd_line)
        self.input = in_p
        self.in_len = len(self.input)
        self.dir = str(work_dir)
        self.shell = shell
        self.wait = wait

    def run(self):
        '''
        replace cmd with args in input
        subprocess.Popen
        :return:
        '''

        import subprocess, os

        run_str = self.cmd
        for i in range(0, self.in_len) :
            run_str = run_str.replace("CMDARG" + str(int(i) + 1), str(self.input[i]))

        run_arg = run_str.split(" ")            # transform to list for subprocess
        os.chdir(self.dir)          # change the working dir

        # write run log
        fh = open('log_run', 'a')
        print('Run the command:', file=fh)
        print(run_str, file=fh)
        fh.close()
        run_subp  = subprocess.Popen(run_arg, shell=self.shell)

        if self.wait :
            run_subp.wait()

class cmd_sh(cmd) :

    def run(self) :
        '''
        run in shell communication mode
        return shell out for the first line
        '''

        import subprocess, os

        run_str = self.cmd
        for i in range(0, self.in_len) :
            run_str = run_str.replace("CMDARG" + str(int(i) + 1), str(self.input[i]))

        os.chdir(self.dir)          # change the working dir

        # write run log
        fh = open('log_run', 'a')
        print('Run the command:', file=fh)
        print(run_str, file=fh)
        fh.close()

        run_subp = subprocess.Popen(run_str, stdout=subprocess.PIPE, shell=self.shell)
        run_subp.wait()
        out = run_subp.stdout.read().decode().strip()

        return out

class cmd_out(cmd) :
    def __init__(self, cmd_line, in_p, work_dir, outfile, shell=False, wait=True) :
        cmd.__init__(self,
                     cmd_line=cmd_line,
                     in_p=in_p,
                     work_dir=work_dir,
                     wait=wait,
                     shell=shell
                     )
        self.outfile = outfile

    def run(self) :

        import subprocess, os

        run_str = self.cmd
        for i in range(0, self.in_len) :
            run_str = run_str.replace("CMDARG" + str(int(i) + 1), str(self.input[i]))

        os.chdir(self.dir)          # change the working dir

        # write run log
        fh = open('log_run', 'a')
        print('Run the command:', file=fh)
        print(run_str, file=fh)
        fh.close()

        with open(self.outfile, 'w') as ofh :
            subprocess.call(run_str, stdout=ofh, shell=True)
