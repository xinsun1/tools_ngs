__author__ = 'sunxin'

import os
import time
import subprocess

class file :
    def __init__(self, name, short_name, work_dir) :
        '''
        name
           file name, no suffix, FILE_1.fq
        short_name
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
            |-4-population_genetics
            |
        '''

        self.name = name
        self.short_name = short_name
        self.wdir = work_dir

    def cutada(self, qual=30, min=40, trim_start=False) :
        '''
        use cutadapt to trim the raw file,
        creat self.trim_file
        qual
            trim quality, default=30
        min
            trim minimum length, default=40
        '''

        os.chdir(self.wdir)
        subprocess.Popen(['mkdir', '1-trim_fq'])
        time.sleep(3)
        os.chdir(self.wdir + '/1-trim_fq')

        if trim_start :
            trim_arg = ['cutadapt',
                        '-q', str(qual) + ',' + str(qual),
                        '-m', str(min),
                        '-a', 'AGATCGGAAGAGC',
                        '-A', 'AGATCGGAAGAGC',
                        '-u', '10',
                        '-o', 'trim_' + self.short_name + '_1.fq',
                        '-p', 'trim_' + self.short_name + '_2.fq',
                        self.wdir + '/0-raw/' + self.name + '_1.fq.gz',
                        self.wdir + '/0-raw/' + self.name + '_2.fq.gz'
                        ]
        else :
            trim_arg = ['cutadapt',
                        '-q', str(qual) + ',' + str(qual),
                        '-m', str(min),
                        '-a', 'AGATCGGAAGAGC',
                        '-A', 'AGATCGGAAGAGC',
                        '-o', 'trim_' + self.short_name + '_1.fq',
                        '-p', 'trim_' + self.short_name + '_2.fq',
                        self.wdir + '/0-raw/' + self.name + '_1.fq.gz',
                        self.wdir + '/0-raw/' + self.name + '_2.fq.gz'
                        ]

        cuta = subprocess.Popen(trim_arg)
        cuta.wait()

    def cutada_ss(self, qual=30, min=40) :
        '''
        use cutadapt to trim the raw file,
        ONLY cut adapter without quality cut
        creat self.trim_ada_file
        min
            trim minimum length, default=40
        '''

        os.chdir(self.wdir)
        subprocess.Popen(['mkdir', '1-trim_fq'])
        time.sleep(3)
        os.chdir(self.wdir + '/1-trim_fq')

        trim_arg1 = ['cutadapt', '-u', '1',
                     '-o', 'trim_C_' + self.short_name + '_C_1.fq',
                     self.wdir + '/0-raw/' + self.name + '_1.fq.gz']

        cuta1 = subprocess.Popen(trim_arg1)
        cuta1.wait()

        trim_arg2 = ['cutadapt',
                    '-q', str(qual) + ',' + str(qual),
                    '-m', str(min),
                    '-a', 'AGATCGGAAGAGC',
                    '-A', 'GAGATCGGAAGAGC',
                    '-o', 'trim_' + self.short_name + '_1.fq',
                    '-p', 'trim_' + self.short_name + '_2.fq',
                    'trim_C_' + self.short_name + '_C_1.fq',
                    self.wdir + '/0-raw/' + self.name + '_2.fq.gz'
                    ]

        cuta2 = subprocess.Popen(trim_arg2)
        cuta2.wait()


    def fqc(self) :
        '''
        run fastqc for trimmed files
        '''

        os.chdir(self.wdir + '/1-trim_fq')

        subprocess.Popen(['fastqc', '-t', '20',
                          '--extract',
                          '-o', './',
                          'trim_' + self.short_name + '_1.fq',
                          'trim_' + self.short_name + '_2.fq',
                          '&'])
        # subprocess.Popen(['rm',
        #                   'trim_' + self.short_name + '_1_fastqc.zip',
        #                   'trim_' + self.short_name + '_2_fastqc.zip'])'


    def mt_map(self, t=10) :
        '''
        map sequence to mitochondria
        map -> sam 2 bam -> sort -> rmdup

        t
            threads to use for bwa mem
        '''

        os.chdir(self.wdir)
        subprocess.Popen(['mkdir', '2-map'])
        time.sleep(3)
        os.chdir(self.wdir + '/2-map')
        subprocess.Popen(['mkdir', '1-mt'])
        time.sleep(3)
        os.chdir(self.wdir + '/2-map/1-mt')

        with open(self.short_name + '_mt.sam', 'w') as ofh :
            a = ['bwa', 'mem',
                 '-t', str(t),
                 '$TIGER_MT',
                 '../../1-trim_fq/trim_' + self.short_name + '_1.fq',
                 '../../1-trim_fq/trim_' + self.short_name + '_2.fq']
            subprocess.call(' '.join(a), stdout=ofh, shell=True)

        s2b = subprocess.Popen(['samtools', 'view', '-b', '-h',
                                '-q', '1',
                                '-F', '3852',
                                '-f', '3',
                                self.short_name + '_mt.sam',
                                '-o', self.short_name + '_mt_f.bam'])
        s2b.wait()

        b_sort = subprocess.Popen(['samtools', 'sort',
                                   '-@', '20',
                                   '-o', self.short_name + '_mt_f_sort.bam',
                                   self.short_name + '_mt_f.bam'
                                   ])
        b_sort.wait()

        rmdup = subprocess.Popen(['samtools', 'rmdup',
                                  '-s', self.short_name + '_mt_f_sort.bam',
                                  self.short_name + '_mt_f_sort_rmdup.bam'
                                  ])
        rmdup.wait()

        subprocess.Popen(['rm', self.short_name + '_mt_f.bam'])
        subprocess.Popen(['rm', self.short_name + '_mt_f_sort.bam'])

    def wg_map(self, t=10) :
        '''
        map sequence to whole genome
        map -> sam 2 bam -> sort -> rmdup

        t
            threads to use for bwa mem
        '''

        os.chdir(self.wdir)
        subprocess.Popen(['mkdir', '2-map'])
        time.sleep(3)
        os.chdir(self.wdir + '/2-map')
        subprocess.Popen(['mkdir', '2-wg'])
        time.sleep(3)
        os.chdir(self.wdir + '/2-map/2-wg')

        with open(self.short_name + '_wg.sam', 'w') as ofh :
            a = ['bwa', 'mem', '-t', str(t),
                '$TIGER_WG',
                '../../1-trim_fq/trim_' + self.short_name + '_1.fq',
                '../../1-trim_fq/trim_' + self.short_name + '_2.fq']
            subprocess.call(' '.join(a), stdout=ofh, shell=True)

        s2b = subprocess.Popen(['samtools', 'view', '-b', '-h',
                                '-q', '1',
                                '-F', '3852',
                                '-f', '3',
                                self.short_name + '_wg.sam',
                                '-o', self.short_name + '_wg_f.bam'])
        s2b.wait()

        b_sort = subprocess.Popen(['samtools', 'sort',
                                   '-@', '20',
                                   '-o', self.short_name + '_wg_f_sort.bam',
                                   self.short_name + '_wg_f.bam'
                                   ])
        b_sort.wait()

        rmdup = subprocess.Popen(['samtools', 'rmdup',
                                  '-s', self.short_name + '_wg_f_sort.bam',
                                  self.short_name + '_wg_f_sort_rmdup.bam'
                                  ])
        rmdup.wait()

        subprocess.Popen(['rm', self.short_name + '_wg_f.bam'])
        subprocess.Popen(['rm', self.short_name + '_wg_f_sort.bam'])
















