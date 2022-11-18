__author__ = 'sunxin'

from lib_cmd import *

'''
file lib for pipeline

consider each sample as a file object.

Last edit : 20190321

'''


class file:
    def __init__(self, name, short_name, work_dir, cwdir=True) :
        '''
        name
           file name, no suffix, FILE1.fq.gz, FILE2.fq.gz
        short_name
        work_dir
        cwdir = True
        '''

        self.name = name
        self.short_name = short_name
        self.wdir = work_dir

        if cwdir :
            self.init_dir()

    def init_dir(self) :
        '''
        work_dir
            WORK_DIR(full dir)
            |
            |-file_name_list
            |-0-raw
            |-1-trim_fq
            |-2-map
                |
                |-1-mt
                |-2-wg
            |-3-vcf
            # |-4-population_genetics
            |
        '''

        import subprocess, os, time
        os.chdir(self.wdir)

        if not  '0-raw' in os.listdir() :
            subprocess.Popen(['mkdir', '0-raw'])
            time.sleep(1)
        if not  '1-trim_fq' in os.listdir() :
            subprocess.Popen(['mkdir', '1-trim_fq'])
            time.sleep(1)
        if not  '2-map' in os.listdir() :
            subprocess.Popen(['mkdir', '2-map'])
            time.sleep(1)
        if not  '3-vcf' in os.listdir() :
            subprocess.Popen(['mkdir', '3-vcf'])
            time.sleep(1)
        os.chdir(self.wdir + '/2-map')
        if not  '1-mt' in os.listdir() :
            subprocess.Popen(['mkdir', '1-mt'])
            time.sleep(1)
        if not  '2-wg' in os.listdir() :
            subprocess.Popen(['mkdir', '2-wg'])
            time.sleep(1)


    def cuta(self, t=10, qual=30, min=25, trim_start=False) :
        '''
        use cutadapt to trim the raw file,
        creat self.trim_file
        t
            threads to use
        qual
            trim quality, default=30
        min
            trim minimum length, default=25
        trim_start
            take as 3,-3 trim 3bp from both end
        '''

        if trim_start :
            trim_s = trim_start.split(',')[0]   #  trim from left
            trim_e = trim_start.split(',')[1]   #  trim from end
            cuta_O = 3 + int(trim_e)

            cuta = cmd(cmd_line=('cutadapt -j CMDARG1 -q CMDARG2,CMDARG2 -m CMDARG3 '
                                 '-a N{CMDARG7}AGATCGGAAGAGC -A N{CMDARG7}AGATCGGAAGAGC '
                                 '-u CMDARG6 '
                                 '-O CMDARG8 '
                                 '-o trim_CMDARG4_1.fq -p trim_CMDARG4_2.fq '
                                 '../0-raw/CMDARG5_1.fq.gz ../0-raw/CMDARG5_2.fq.gz'),
                       in_p=[t, qual, min, self.short_name, self.name,
                             trim_s, trim_e, cuta_O],
                       work_dir=self.wdir + '/1-trim_fq',
                       shell=False)
        else :
            cuta = cmd(cmd_line=('cutadapt -j CMDARG1 -q CMDARG2,CMDARG2 -m CMDARG3 '
                                 '-a AGATCGGAAGAGC -A AGATCGGAAGAGC '
                                 '-o trim_CMDARG4_1.fq -p trim_CMDARG4_2.fq '
                                 '../0-raw/CMDARG5_1.fq.gz ../0-raw/CMDARG5_2.fq.gz'),
                       in_p=[t, qual, min, self.short_name, self.name],
                       work_dir=self.wdir + '/1-trim_fq',
                       shell=False)
        cuta.run()


    def cuta_ss(self, t=10, qual=30, min=20, trim_start=False) :
        '''
        use cutadapt to trim the raw file,
        creat self.trim_file
        qual
            trim quality, default=30
        min
            trim minimum length, default=20
        trim_start
            take as 3,-3 trim 3bp from both end
        '''

        rmc = cmd(cmd_line=('cutadapt -j 20 -u 1 -o trim_C_CMDARG1_C_1.fq '
                            '../0-raw/CMDARG2_1.fq.gz'),
                  in_p=[self.short_name, self.name],
                  work_dir=self.wdir + '/1-trim_fq',
                  shell=False)
        rmc.run()

        if trim_start:
            trim_s = trim_start.split(',')[0]  # trim from left
            trim_e = trim_start.split(',')[1]  # trim from end
            cuta_O = 3 + int(trim_e)

            cuta = cmd(cmd_line=('cutadapt -j CMDARG1 -q CMDARG2,CMDARG2 -m CMDARG3 '
                                 '-a N{CMDARG7}AGATCGGAAGAGC -A N{CMDARG7}GAGATCGGAAGAGC '
                                 '-u CMDARG6 '
                                 '-O CMDARG8 '
                                 '-o trim_CMDARG4_1.fq -p trim_CMDARG4_2.fq '
                                 'trim_C_CMDARG4_C_1.fq ../0-raw/CMDARG5_2.fq.gz'),
                       in_p=[t, qual, min, self.short_name, self.name,
                             trim_s, trim_e, cuta_O],
                       work_dir=self.wdir + '/1-trim_fq',
                       shell=False)
        else :
            cuta = cmd(cmd_line=('cutadapt -j CMDARG 5 -q CMDARG1,CMDARG1 -m CMDARG2 '
                                 '-a AGATCGGAAGAGC -A GAGATCGGAAGAGC '
                                 '-o trim_CMDARG3_1.fq -p trim_CMDARG3_2.fq '
                                 'trim_C_CMDARG3_C_1.fq ../0-raw/CMDARG4_2.fq.gz'),
                       in_p=[qual, min, self.short_name, self.name, t],
                       work_dir=self.wdir + '/1-trim_fq',
                       shell=False)
        cuta.run()


    def adar(self, threads=10) :
        '''
        collapse pair-end reads to single end

        out files
            SHORT_NAME.collapsed
            SHORT_NAME.pair1.truncated
            SHORT_NAME.pair2.truncated

        :return:
        '''

        adar = cmd(cmd_line=('AdapterRemoval --file1 trim_CMDARG1_1.fq '
                             '--file2 trim_CMDARG1_2.fq --minalignmentlength 20 '
                             '--threads CMDARG2 --collapse --basename CMDARG1'),
                   in_p=[self.short_name, str(threads)],
                   work_dir=self.wdir + '/1-trim_fq',
                   shell=False)
        adar.run()


    def map(self, map_tool, ref, dir_str, fq, type, t, sam_f, mem_arg=''):
        '''
        mapping function
        :param ref:
         reference
        :param dir:
         working dir
        :param fq:
         fq for map
        :param type:
         mt or wg or clp_mt or clp_wg
        :param t:
         threads for mapping
        :param sam_f:
         -F 3852 -f 3
        :param mem_arg:
         args for bwa mem
        :return:
        '''

        if map_tool == 'mem' :
            map = cmd_out(cmd_line=('bwa mem -t CMDARG1 CMDARG2 CMDARG3 CMDARG4'),  #[threads, ref, fq]
                          in_p=[t, mem_arg, ref, fq],
                          work_dir=dir_str,
                          outfile=self.short_name + '_' + str(type) + '.sam',
                          shell=True)
            map.run()
        elif map_tool == 'aln' :
            map = cmd_out(cmd_line='bwa aln -l 16500 -t CMDARG1 CMDARG2 CMDARG3',
                          in_p=[t, ref, fq],
                          work_dir=dir_str,
                          outfile=self.short_name + '_' + str(type) + '.sai',
                          shell=True)
            map.run()
            map2 = cmd_out(cmd_line='bwa samse CMDARG1 CMDARG2 CMDARG3',
                           in_p=[ref, self.short_name + '_' + str(type) + '.sai', fq],
                           work_dir=dir_str,
                           outfile=self.short_name + '_' + str(type) + '.sam',
                           shell=True)
            map2.run()

        s2b = cmd(cmd_line=('samtools view -b -h -q 20 CMDARG2 ' #[short_name, sam_f, type]
                            'CMDARG1_CMDARG3.sam -o CMDARG1_CMDARG3_f.bam'),
                  in_p=[self.short_name, sam_f, str(type)],
                  work_dir=dir_str,
                  shell=False)
        s2b.run()

        b_nsort = cmd(cmd_line=('samtools sort -n -@ 20 -o CMDARG1_CMDARG2_f_nsort.bam ' #[short_name, type]
                               'CMDARG1_CMDARG2_f.bam'),
                     in_p=[self.short_name, str(type)],
                     work_dir=dir_str,
                     shell=False)
        b_nsort.run()


        b_fixmate = cmd(cmd_line=('samtools fixmate -m CMDARG1_CMDARG2_f_nsort.bam ' #[short_name, type]
                               'CMDARG1_CMDARG2_f_nsort_fix.bam'),
                     in_p=[self.short_name, str(type)],
                     work_dir=dir_str,
                     shell=False)
        b_fixmate.run()

        b_sort = cmd(cmd_line=('samtools sort -@ 20 -o CMDARG1_CMDARG2_f_sort.bam ' #[short_name, type]
                               'CMDARG1_CMDARG2_f_nsort_fix.bam'),
                     in_p=[self.short_name, str(type)],
                     work_dir=dir_str,
                     shell=False)
        b_sort.run()

        b_stat1 = cmd_out(cmd_line=('samtools stats CMDARG1_CMDARG2_f_sort.bam'), #[short_name, type]
                     in_p=[self.short_name, str(type)],
                     work_dir=dir_str,
                     wait=True,
                     outfile='stat_' + self.short_name + '_' + str(type) + '_f_sort.bam',
                     shell=False)
        b_stat1.run()

        rmdup = cmd(cmd_line=('samtools markdup -r CMDARG1_CMDARG2_f_sort.bam ' #[short_name, type]
                              'CMDARG1_CMDARG2_f_sort_rmdup.bam'),
                    in_p=[self.short_name, str(type)],
                    work_dir=dir_str,
                    shell=False)
        rmdup.run()

        b_index = cmd(cmd_line=('samtools index CMDARG1_CMDARG2_f_sort_rmdup.bam'), #[short_name, type]
                      in_p=[self.short_name, str(type)],
                      work_dir=dir_str,
                      wait=False,
                      shell=False)
        b_index.run()

        rm_1 = cmd(cmd_line=('rm CMDARG1_CMDARG2_f.bam'), #[short_name, type]
                   in_p=[self.short_name, str(type)],
                   work_dir=dir_str,
                   wait=False,
                   shell=False)
        rm_1.run()

        rm_2 = cmd(cmd_line=('rm CMDARG1_CMDARG2_f_sort.bam'), #[short_name, type]
                   in_p=[self.short_name, str(type)],
                   work_dir=dir_str,
                   wait=False,
                   shell=False)
        rm_2.run()

        rm_3 = cmd(cmd_line=('rm CMDARG1_CMDARG2_f_nsort.bam'), #[short_name, type]
                   in_p=[self.short_name, str(type)],
                   work_dir=dir_str,
                   wait=False,
                   shell=False)
        rm_3.run()

        rm_4 = cmd(cmd_line=('rm CMDARG1_CMDARG2_f_nsort_fix.bam'), #[short_name, type]
                   in_p=[self.short_name, str(type)],
                   work_dir=dir_str,
                   wait=False,
                   shell=False)
        rm_4.run()


        b_stat = cmd_out(cmd_line=('samtools stats CMDARG1_CMDARG2_f_sort_rmdup.bam'), #[short_name, type]
                     in_p=[self.short_name, str(type)],
                     work_dir=dir_str,
                     wait=False,
                     outfile='stat_' + self.short_name + '_' + str(type) + '_f_sort_rmdup.bam',
                     shell=False)
        b_stat.run()


    def rmdup(self) :
        # duplication remove using samtools version 1.7
        pass


    def cuta_se(self, t=10, qual=30, min=20, trim_start=False) :
        '''
        use cutadapt to trim the raw file,
        creat self.trim_file
        t
            threads to use
        qual
            trim quality, default=30
        min
            trim minimum length, default=20
        trim_start
            trim sequence from both end
            take 3,3 for 3bp from both end
        '''


        if trim_start :
            trim_s = trim_start.split(',')[0]  # trim from left
            trim_e = trim_start.split(',')[1]  # trim from end
            cuta_O = 3 + int(trim_e)

            cuta = cmd(cmd_line=('cutadapt -j CMDARG1 -q CMDARG2 -m CMDARG3 '
                                 '-a N{CMDARG7}AGATCGGAAGAGC '
                                 '-u CMDARG6 '
                                 '-O CMDARG8 '
                                 '-o trim_CMDARG4.fq '
                                 '../0-raw/CMDARG5.fq.gz'),
                       in_p=[t, qual, min, self.short_name, self.name,
                             trim_s, trim_e, cuta_O],
                       work_dir=self.wdir + '/1-trim_fq',
                       shell=False)
        else :
            cuta = cmd(cmd_line=('cutadapt -j CMDARG1 -q CMDARG2 -m CMDARG3 '
                                 '-a AGATCGGAAGAGC '
                                 '-o trim_CMDARG4.fq '
                                 '../0-raw/CMDARG5.fq.gz'),
                       in_p=[t, qual, min, self.short_name, self.name],
                       work_dir=self.wdir + '/1-trim_fq',
                       shell=False)
        cuta.run()

    def map_se(self, mt_ref, wg_ref, mem_arg, t=10):
        '''
        map pipeline for single end sequencing
        '''

        self.map(ref=mt_ref,
                 map_tool='mem',
                 dir_str=self.wdir + '/2-map/1-mt',
                 fq='../../1-trim_fq/trim_' + self.short_name + '.fq ',
                 type='mt',
                 t=t,
                 sam_f='-F 3852')

        self.map(ref=wg_ref,
                 map_tool='mem',
                 dir_str=self.wdir + '/2-map/2-wg',
                 fq='../../1-trim_fq/trim_' + self.short_name + '.fq ',
                 type='wg',
                 t=t,
                 sam_f='-F 3852',
                 mem_arg=mem_arg)



    def mt_map(self, ref, t=10) :
        '''
        map sequence to mitochondria
        map -> sam 2 bam -> sort -> rmdup

        :param t:
        threads use for map
        :return:
        '''

        self.map(ref=ref,
                 map_tool='mem',
                 dir_str=self.wdir + '/2-map/1-mt',
                 fq='../../1-trim_fq/trim_' + self.short_name + '_1.fq ' +
                    '../../1-trim_fq/trim_' + self.short_name + '_2.fq',
                 type='mt',
                 t=t,
                 sam_f='-F 3852 -f 3')


    def wg_map(self, ref, mem_arg, t=10) :
        '''
        map sequence to whole genome
        map -> sam 2 bam -> sort -> rmdup

        :param t:
         threads to use for map
        :return:
        '''

        self.map(ref=ref,
                 map_tool='mem',
                 dir_str=self.wdir + '/2-map/2-wg',
                 fq='../../1-trim_fq/trim_' + self.short_name + '_1.fq ' +
                    '../../1-trim_fq/trim_' + self.short_name + '_2.fq',
                 type='wg',
                 t=t,
                 sam_f='-F 3852 -f 3',
                 mem_arg=mem_arg)


    def collapse_mt(self, ref, t=10) :
        '''
        map collapsed sequence to mitochondria
        map -> sam 2 bam -> sort -> rmdup

        :param t:
        threads use for map
        :return:
        '''

        self.map(ref=ref,
                 map_tool='aln',
                 dir_str=self.wdir + '/2-map/1-mt',
                 fq='../../1-trim_fq/' + self.short_name + '.collapsed',
                 type='clp_mt',
                 t=t,
                 sam_f='-F 3852')


    def collapse_wg(self, ref, t=10) :
        '''
        map collapsed sequence to whole genome
        map -> sam 2 bam -> sort -> rmdup

        :param t:
         threads to use for map
        :return:
        '''

        self.map(ref=ref,
                 map_tool='aln',
                 dir_str=self.wdir + '/2-map/2-wg',
                 fq='../../1-trim_fq/' + self.short_name + '.collapsed',
                 type='clp_wg',
                 t=t,
                 sam_f='-F 3852')


    def read_stat(self, wdir, type):
        '''
        function for read samtools stats file
        :param wdir:
         work dir
        :param type:
         mt, wg, clp_mt, clp_wg
        :return:
        [read count, base count, error rate]
        '''

        run_read_c = cmd_sh(cmd_line='grep \'reads mapped:\' stat_CMDARG1_CMDARG2_f_sort_rmdup.bam',
                                                #[short_name, type]
                        in_p=[self.short_name, type],
                        work_dir=wdir,
                        shell=True)
        try :
            read_c_tmp = run_read_c.run()
            read_c = int(read_c_tmp.split('\t')[2])
        except IndexError :
            read_c = 0
        run_base_c = cmd_sh(cmd_line='grep \'bases mapped (cigar):\' stat_CMDARG1_CMDARG2_f_sort_rmdup.bam',
                                                #[short_name, type]
                        in_p=[self.short_name, type],
                        work_dir=wdir,
                        shell=True)
        try :
            base_c_tmp = run_base_c.run()
            base_c = int(base_c_tmp.split('\t')[2])
        except IndexError :
            base_c = 0

        run_err_r = cmd_sh(cmd_line='grep \'error rate:\' stat_CMDARG1_CMDARG2_f_sort_rmdup.bam',
                                                #[short_name, type]
                        in_p=[self.short_name, type],
                        work_dir=wdir,
                        shell=True)
        try :
            err_r_tmp = run_err_r.run()
            err_r = str(err_r_tmp.split('\t')[2])
        except IndexError :
            err_r = '0'

        run_base_c1 = cmd_sh(cmd_line='grep \'bases mapped (cigar):\' stat_CMDARG1_CMDARG2_f_sort.bam',
                                                #[short_name, type]
                        in_p=[self.short_name, type],
                        work_dir=wdir,
                        shell=True)
        try :
            base_c1_tmp = run_base_c1.run()
            base_c1 = int(base_c1_tmp.split('\t')[2])
        except IndexError :
            base_c1 = 0


        return [read_c, base_c, err_r, base_c1]


    def map_stat(self, clp=False, mt_len=16793, wg_len=2332967107, se=False) :
        '''
        mapping result statistics

        :return:
        list[sample ID,
         Raw Read count, Raw read base,
         ada cut base, ada cut %,
         mt_reads, mt base, error rate, aver_mt_depth, mt_raw %,
         wg_reads, wg base, aver_len, error_rate, aver_wg_depth, wg_raw %,
         mt/wg, wg_ada %]

        '''

        map_out = [self.short_name]


        #### read raw data info ####
        c_raw_c = cmd_sh(cmd_line='zcat CMDARG1_1.fq.gz | wc -l',
                         in_p=[self.name],
                         work_dir=self.wdir + '/0-raw',
                         shell=True)
        try :
            raw_read_c = int(c_raw_c.run())
        except FileNotFoundError :
            print(self.short_name + ' raw fq not found use 0 ')
            raw_read_c = 0

        raw_read_base = int(raw_read_c / 4 * 2 * 150)
        map_out += [raw_read_c, raw_read_base]


        #### read after ada cut info ####
        ada_cut_c1 = cmd_sh(cmd_line=('awk \'BEGIN{a=0}{if(NR%4==2){a+=length($0)}}END{print a}\' '
                                     'trim_CMDARG1_1.fq'),      # [short_name]
                           in_p=[self.short_name],
                           work_dir=self.wdir + '/1-trim_fq',
                           shell=True)
        ada_cut_c2 = cmd_sh(cmd_line=('awk \'BEGIN{a=0}{if(NR%4==2){a+=length($0)}}END{print a}\' '
                                     'trim_CMDARG1_2.fq'),      # [short_name]
                           in_p=[self.short_name],
                           work_dir=self.wdir + '/1-trim_fq',
                           shell=True)
        try :
            ada_cut_c1_tmp = ada_cut_c1.run()
            ada_cut_c2_tmp = ada_cut_c2.run()
            ada_cut_c = int(ada_cut_c1_tmp) + int(ada_cut_c2_tmp)
        except ValueError :
            print(self.short_name + ' adaptor cut fq not found use 0 ')
            ada_cut_c = 0

        try :
            ada_cut_raw = round(float(ada_cut_c) / float(raw_read_base), 5)
        except ZeroDivisionError :
            print(self.short_name + ' divided adaptor proportion by 0')
            ada_cut_raw = 0
        map_out += [ada_cut_c, ada_cut_raw]


        #### read mappped info ####
        if clp :
            mt = self.read_stat(wdir=self.wdir + '/2-map/1-mt',
                                 type='clp_mt')
            wg = self.read_stat(wdir=self.wdir + '/2-map/2-wg',
                                 type='clp_wg')
        else :
            mt = self.read_stat(wdir=self.wdir + '/2-map/1-mt',
                                type='mt')
            wg = self.read_stat(wdir=self.wdir + '/2-map/2-wg',
                                type='wg')

        mt_dp = round(float(mt[1])/float(mt_len), 5)
        wg_dp = round(float(wg[1])/float(wg_len), 5)

        # [mt_reads, mt base, error rate, aver_mt_depth, mt_raw %,] #
        try :
            map_out += [mt[0], mt[1], mt[2], mt_dp,
                        round(float(mt[1])/float(raw_read_base), 5)]
        except ZeroDivisionError :
            map_out += [mt[0], mt[1], mt[2], mt_dp, 0]

        # [wg_reads, wg base, aver_len, error_rate, aver_wg_depth, wg_raw %,] #
        try :
            map_out += [wg[0], wg[1],
                        round(float(wg[1])/float(wg[0]), 3) if float(wg[0]) != 0 else 0,
                        wg[2], wg_dp,
                        round(float(wg[1])/float(raw_read_base), 5)]
        except ZeroDivisionError :
            map_out += [wg[0], wg[1],
                        round(float(wg[1])/float(wg[0]), 3) if float(wg[0]) != 0 else 0,
                        wg[2], wg_dp,
                        0]

        # [mt/wg, wg_ada %] #
        try :
            map_out += [round(mt_dp/wg_dp, 3)]
        except ZeroDivisionError :
            map_out += [0]
        try :
            map_out += [round(float(wg[1])/float(ada_cut_c), 5)]
        except ZeroDivisionError :
            map_out += [0]

        try :
            map_out += [round(float(wg[3])/float(wg[1]), 5)] #[wg_before_rmdup/wg_base %]
        except ZeroDivisionError :
            map_out += [0]



        import os
        os.chdir(self.wdir)
        fh = open('map_summary', 'a')
        print('\t'.join(str(i) for i in map_out), file=fh)
        fh.close()


    def diy(self, work_dir, cmd_str, in_p, wait=True, out_file=False, shell=False) :
        '''

        :param work_dir:
        cmd work directory
        work_directory : self.wdir + work dir '/1-map'
        :param cmd_str:
        cmd : xxxx CMDARG1 CMDARG2
        str within cmd SHORT_NAME NAME
        :param in_p:
        str sep=',,'
        CMDARG1,,CMDARG2
        :param out_file:
        default FALSE, out_file_name
        :param shell:
        if run with shell mode
        :return:
        '''

        run_str = cmd_str
        run_str = run_str.replace("SHORT_NAME", self.short_name)
        run_str = run_str.replace("NAME", self.name)

        if out_file :
            out_f_str = out_file
            out_f_str = out_f_str.replace("SHORT_NAME", self.short_name)
            out_f_str = out_f_str.replace("NAME", self.name)

            diy_r = cmd_out(cmd_line=run_str,
                            in_p=str(in_p).split(',,'),
                            work_dir=self.wdir + work_dir,
                            outfile=out_f_str,
                            wait=wait,
                            shell=shell)
            diy_r.run()
        else :
            diy_r = cmd(cmd_line=run_str,
                        in_p=str(in_p).split(',,'),
                        work_dir=self.wdir + work_dir,
                        wait=wait,
                        shell=shell)
            diy_r.run()



