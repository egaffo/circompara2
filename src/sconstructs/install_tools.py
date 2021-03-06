import os, sys

def SymLink(target, source, env):
    os.symlink(os.path.abspath(str(source[0])), os.path.abspath(str(target[0])))

env = Environment(ENV=os.environ, SHELL = '/bin/bash')
tools_dir = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'tools')
bin_dir = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'bin')

#env.PrependENVPath('PATH', os.path.join(tools_dir, 'bin'))

PYTHON_VERSION = str(sys.version_info.major) + '.' + str(sys.version_info.minor)

python_lib_dir = os.path.join(tools_dir, 'lib', 'python'+PYTHON_VERSION,'site-packages')
SET_PIP_USER = '--user'

## initialize VIRTUAL_ENV
if not env['ENV'].has_key('VIRTUAL_ENV'):
    env['ENV']['VIRTUAL_ENV'] = ''

## pip does not work with --user if run inside a virtualenv
if env['ENV']['VIRTUAL_ENV']:
    print 'Running in virtualenv. Python packages will be installed in the virtualenv directories'
    SET_PIP_USER = ''
    python_lib_dir = os.path.join(env['ENV']['VIRTUAL_ENV'], 'lib', 
                                  'python' + PYTHON_VERSION, 'site-packages')
else:
    env.PrependENVPath('PYTHONUSERBASE', tools_dir)

env.PrependENVPath('PYTHONPATH', python_lib_dir)

## PIP
## use virtualenv pip if run inside virtualenv
if env['ENV']['VIRTUAL_ENV']:
    pip_targets = os.path.join(env['ENV']['VIRTUAL_ENV'], 'bin', 'pip2')
    pip = env.Command(os.path.join(bin_dir, 'pip'), pip_targets, SymLink)
else:
    pip_file = 'get-pip.py'
    pip_url = 'https://bootstrap.pypa.io/pip/2.7/' + pip_file
    pip_targets = [os.path.join(tools_dir, pip_file),
                   os.path.join(tools_dir, 'bin', 'pip')]
    pip_cmd = ' && '.join(['wget -O ${TARGETS[0]} ' + pip_url, 
                           'python ${TARGETS[0]} ' + SET_PIP_USER + ' pip==20.1.1'])
    pip = env.Command(pip_targets, 
                      [], 
                      pip_cmd)
    pip = env.Command(os.path.join(bin_dir, 'pip'), pip[1], SymLink)

# PYSAM
## freeze pysam to v0.15.4 since
## v0.16 does not read gzip'ed files
PYSAM_dir = os.path.join(python_lib_dir, 'pysam')
PYSAM_target = [os.path.join(PYSAM_dir, 'samtools.py')]
PYSAM = env.Command(PYSAM_target, [pip], 
                        ['pip install --ignore-installed ' + SET_PIP_USER + ' pysam==0.15.4'])


# BIOPYTHON
## N.B: if the virtualenv is not a subdirectory of circompara (installation),
## i.e. not under the Scons directory scope, then the installation could miss
## the targets in the virtualenv directories
BIOPYTHON_dir = os.path.join(python_lib_dir, 'Bio')
BIOPYTHON_target = [os.path.join(BIOPYTHON_dir, 'SeqIO', 'FastaIO.py')]
BIOPYTHON = env.Command(BIOPYTHON_target, [pip, PYSAM], 
                        ['pip install ' + SET_PIP_USER + ' biopython==1.76']) #--ignore-installed 

## CYTHON
#cython_dir = os.path.join(python_lib_dir, 'cython')
if env['ENV']['VIRTUAL_ENV']:
    cython_target = [os.path.join(env['ENV']['VIRTUAL_ENV'], 'bin', 'cython')]
else:
    cython_target = [os.path.join(tools_dir, 'bin', 'cython')]
cython_cmd = 'pip install ' + SET_PIP_USER + ' Cython==0.29.19 --install-option="--no-cython-compile"' #--ignore-installed 
cython = env.Command(cython_target,
                     [pip, PYSAM],
                     cython_cmd)

## HTSeq
HTSeq_dir = os.path.join(python_lib_dir, 'HTSeq')
#HTSeq_target = [os.path.join(HTSeq_dir, 'lib','python2.7','site-packages',
#                'HTSeq', '__init__.py'),
HTSeq_target = [os.path.join(HTSeq_dir, '__init__.py')]
if env['ENV']['VIRTUAL_ENV']:
    HTSeq_target.append(os.path.join(env['ENV']['VIRTUAL_ENV'], 'bin', 'htseq-count'))
else:
    HTSeq_target.append(os.path.join(tools_dir, 'bin', 'htseq-count'))
HTSeq = env.Command(HTSeq_target, [pip, cython, PYSAM], 
                    ['pip install ' + SET_PIP_USER + ' HTSeq==0.12.4']) # --ignore-installed
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), HTSeq[1], SymLink)

# CIRCEXPLORER2
if env['ENV']['VIRTUAL_ENV']:
    CIRCEXPLORER2_target = [os.path.join(env['ENV']['VIRTUAL_ENV'], 'bin', 'CIRCexplorer2')]#, 
else:
    CIRCEXPLORER2_target = [os.path.join(tools_dir, 'bin', 'CIRCexplorer2')]#, 
CIRCEXPLORER2 = env.Command(CIRCEXPLORER2_target, [pip, HTSeq, PYSAM],
                           ['pip install ' + SET_PIP_USER + ' circexplorer2==2.3.8']) #--ignore-installed not added because it will install pysam > v0.15.4 that does not support reading gzip'd files
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), CIRCEXPLORER2[0], SymLink)

# DCC
## PANDAS (DCC dependency)
pandas_dir = os.path.join(python_lib_dir, 'pandas')
pandas = env.Command([os.path.join(pandas_dir, '__init__.py')], 
                         [pip, HTSeq, PYSAM], 
                         'pip install ' + SET_PIP_USER + ' pandas==0.23') #--ignore-installed 

dcc_dir = os.path.join(tools_dir, 'DCC-0.4.8')
dcc_tar = 'v0.4.8.tar.gz'
dcc_url = 'https://github.com/dieterich-lab/DCC/archive/' + dcc_tar
dcc_target = [os.path.join(tools_dir, dcc_tar)]
if env['ENV']['VIRTUAL_ENV']:
    dcc_target.append(os.path.join(env['ENV']['VIRTUAL_ENV'], 'bin', 'DCC'))
else:
    dcc_target.append(os.path.join(tools_dir, 'bin', 'DCC'))
dcc_target.append(os.path.join(python_lib_dir, 'DCC-0.4.8-py'+ PYTHON_VERSION +'.egg'))
dcc = env.Command(dcc_target, 
                  [pandas, PYSAM], 
                  ['wget -O ${TARGETS[0]} ' + dcc_url,
                   'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                   'cd ' + dcc_dir + ' && '\
                   'cp ' + os.path.join(env['ENV']['CIRCOMPARA_HOME'], 
					'src', 'utils', 'python', 
					'DCC_patch_CombineCounts.py') +\
                   ' ' + os.path.join('DCC', 'CombineCounts.py') + ' && '\
                   'python setup.py install ' + SET_PIP_USER + '',
                   'cd ' + Dir('#').abspath]
                  )

dcc_link = env.Command(os.path.join(bin_dir, "${SOURCE.file}"), dcc[1], SymLink)

# HISAT2
## NB: HISAT2 v2.2.0 is not compatible because of issue #255
HISAT2_zip = 'hisat2-2.1.0-Linux_x86_64.zip'
HISAT2_link = 'https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download'
HISAT2_target = [os.path.join(tools_dir, HISAT2_zip), 
                 os.path.join(tools_dir, 'hisat2-2.1.0', 'hisat2'),
                 os.path.join(tools_dir, 'hisat2-2.1.0', 'hisat2-build')]
#HISAT2_source = HISAT2_link
HISAT2 = env.Command(HISAT2_target, [], ['wget -O ${TARGETS[0]} ' + HISAT2_link,
                                         'unzip -o -u -d ${TARGETS[0].dir} ${TARGETS[0]}'])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), HISAT2[1], SymLink)
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), HISAT2[2], SymLink)

# BOWTIE2
BOWTIE2_zip = 'bowtie2-2.4.1-linux-x86_64.zip'
BOWTIE2_link = 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.4.1/' + BOWTIE2_zip
BOWTIE2_target = [os.path.join(tools_dir, BOWTIE2_zip), 
                  os.path.join(tools_dir, 'bowtie2-2.4.1-linux-x86_64', 'bowtie2'),
                  os.path.join(tools_dir, 'bowtie2-2.4.1-linux-x86_64', 'bowtie2-inspect'),
                  os.path.join(tools_dir, 'bowtie2-2.4.1-linux-x86_64', 'bowtie2-build')]
#BOWTIE2_source = 
BOWTIE2 = env.Command(BOWTIE2_target, [], ['wget -O ${TARGETS[0]} ' + BOWTIE2_link,
                                           'unzip -o -u -d ${TARGETS[0].dir} ${TARGETS[0]}'])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), BOWTIE2[1], SymLink)
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), BOWTIE2[2], SymLink)
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), BOWTIE2[3], SymLink)

## BWA-MEM
## no binary package available for BWA > v0.7.15 :(
BWAMEM_tar = 'bwakit-0.7.15_x64-linux.tar.bz2'
BWAMEM_link = 'https://sourceforge.net/projects/bio-bwa/files/bwakit/' + BWAMEM_tar
BWAMEM_target = [os.path.join(tools_dir, BWAMEM_tar),
                 os.path.join(tools_dir, 'bwa.kit', 'bwa')]
#BWAMEM_source = 
BWAMEM = env.Command(BWAMEM_target, [], ['wget -O ${TARGETS[0]}  ' + BWAMEM_link,
                                         'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}'])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), BWAMEM[1], SymLink)

# STAR
STAR_tar = '2.6.1e.tar.gz'
STAR_link = 'https://github.com/alexdobin/STAR/archive/' + STAR_tar
STAR_target = [os.path.join(tools_dir, 'STAR_' + STAR_tar), 
               os.path.join(tools_dir, 'STAR-2.6.1e', 'bin', 'Linux_x86_64_static', 'STAR')]
#STAR_source = 
STAR = env.Command(STAR_target, [], ['wget -O ${TARGETS[0]} ' + STAR_link,
                                     'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}'])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), STAR[1], SymLink)

# HTSLIB (required by Segemehl)
HTSLIB_dir = os.path.join(tools_dir, 'htslib-1.10.2')
HTSLIB_tar = 'htslib-1.10.2.tar.bz2'
HTSLIB_target = [os.path.join(tools_dir, HTSLIB_tar),
                 os.path.join(HTSLIB_dir, 'lib', 'pkgconfig', 'htslib.pc')]
HTSLIB_link = 'https://github.com/samtools/htslib/releases/download/1.10.2/' +\
              HTSLIB_tar
HTSLIB = env.Command(HTSLIB_target,
                     [],
                     ['wget -O $TARGET  ' + HTSLIB_link,
                      'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                      ' && '.join(['cd ' + HTSLIB_dir, 
                                   './configure --prefix=`pwd`',
                                   'make', 'make install',
                                   'cd ' + Dir('.').abspath])
                      ])

# SEGEMEHL
SEGEMEHL_tar = 'segemehl-0.3.4.tar.gz'
SEGEMEHL_dir = os.path.join(tools_dir, 'segemehl-0.3.4')
SEGEMEHL_link = 'http://www.bioinf.uni-leipzig.de/Software/segemehl/downloads/' + SEGEMEHL_tar
SEGEMEHL_target = [os.path.join(tools_dir, SEGEMEHL_tar),
                   os.path.join(SEGEMEHL_dir, 'segemehl.x'),
                   os.path.join(SEGEMEHL_dir, 'haarz.x')]

env.PrependENVPath('PKG_CONFIG_PATH',
                   os.path.dirname(HTSLIB[1].abspath))
SEGEMEHL = env.Command(SEGEMEHL_target, 
                       HTSLIB, 
                       ['wget -O $TARGET  ' + SEGEMEHL_link,
                       'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                       'cd ${TARGETS[1].dir} && make all',
                       'cd ' + Dir('.').abspath])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), SEGEMEHL[1], SymLink)
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), SEGEMEHL[2], SymLink)

# TRIMMOMATIC
TRIMMOMATIC_zip = 'Trimmomatic-0.39.zip'
TRIMMOMATIC_link = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/' + TRIMMOMATIC_zip
TRIMMOMATIC_target = [os.path.join(tools_dir, TRIMMOMATIC_zip), 
                      os.path.join(tools_dir, 'Trimmomatic-0.39', 'trimmomatic-0.39.jar')]
#TRIMMOMATIC_source = 
TRIMMOMATIC = env.Command(TRIMMOMATIC_target, [], ['wget -O $TARGET  ' + TRIMMOMATIC_link,
                                                   'unzip -o -u -d ${TARGETS[0].dir} ${TARGETS[0]}'])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), TRIMMOMATIC[1], SymLink)


# FASTQC
FASTQC_zip = 'fastqc_v0.11.9.zip'
FASTQC_link = 'http://www.bioinformatics.babraham.ac.uk/projects/fastqc/' + FASTQC_zip
FASTQC_target = [os.path.join(tools_dir, FASTQC_zip), 
                 os.path.join(tools_dir, 'FastQC', 'fastqc')]
#FASTQC_source = 
FASTQC = env.Command(FASTQC_target, [], ['wget -O $TARGET  ' + FASTQC_link,
                                         'unzip  -o -u -d ${TARGETS[0].dir} ${TARGETS[0]} '\
                                         '&& chmod +x ${TARGETS[1]}'])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), FASTQC[1], SymLink)


# BEDTOOLS
BEDTOOLS_tar = 'bedtools-2.29.2.tar.gz'
BEDTOOLS_dir = os.path.join(tools_dir, 'bedtools2')
BEDTOOLS_link = 'https://github.com/arq5x/bedtools2/releases/download/v2.29.2/' + BEDTOOLS_tar
BEDTOOLS_target = [os.path.join(tools_dir, BEDTOOLS_tar), 
                   os.path.join(BEDTOOLS_dir, 'bin', 'bedtools')]
#BEDTOOLS_source = 
BEDTOOLS = env.Command(BEDTOOLS_target, [], 
                       ['wget -O $TARGET ' + BEDTOOLS_link, 
                        'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}', 
                        'cd ' + BEDTOOLS_dir + ' && make', 
                        'cd ' + Dir('.').abspath])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), BEDTOOLS[1], SymLink)

# SAMTOOLS
SAMTOOLS_tar = 'samtools-1.10.tar.bz2'
SAMTOOLS_dir = os.path.join(tools_dir, 'samtools-1.10')
SAMTOOLS_link = 'https://github.com/samtools/samtools/releases/download/1.10/' + SAMTOOLS_tar
SAMTOOLS_target = [os.path.join(tools_dir, 'samtools-1.10.tar.bz2'),
                   os.path.join(SAMTOOLS_dir, 'samtools')]
#SAMTOOLS_source = 
SAMTOOLS = env.Command(SAMTOOLS_target, [], ['wget -O $TARGET  ' + SAMTOOLS_link,
                                             'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                                             'cd ' + SAMTOOLS_dir + ' && '\
                                             'make prefix=' + SAMTOOLS_dir + ' install',
                                             'cd ' + Dir('.').abspath])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), SAMTOOLS[1], SymLink)

# CUFFLINKS
CUFFLINKS_tar = 'cufflinks-2.2.1.Linux_x86_64.tar.gz'
CUFFLINKS_link = 'http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/' + CUFFLINKS_tar
CUFFLINKS_target = [os.path.join(tools_dir, CUFFLINKS_tar), 
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cufflinks'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffcompare'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffdiff'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffmerge'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffnorm'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'cuffquant'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'gffread'),
                    os.path.join(tools_dir, 'cufflinks-2.2.1.Linux_x86_64', 'gtf_to_sam')]
#CUFFLINKS_source = 
CUFFLINKS = env.Command(CUFFLINKS_target, [], ['wget -O $TARGET  ' + CUFFLINKS_link,
                                               'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}'])
for t2link in CUFFLINKS[1:]:
    env.Command(os.path.join(bin_dir, "${SOURCE.file}"), t2link, SymLink)


# CIRI
CIRI_zip = 'CIRI_v2.0.6.zip'
CIRI_link = 'http://downloads.sourceforge.net/project/ciri/CIRI2/' + CIRI_zip
CIRI_target = [os.path.join(tools_dir, CIRI_zip), 
               os.path.join(tools_dir, 'CIRI_v2.0.6', 'CIRI2.pl')]
#CIRI_source = 
CIRI = env.Command(CIRI_target, [], ['wget -O ${TARGETS[0]} ' + CIRI_link,
                                     'unzip -o -u -d ${TARGETS[0].dir} ${TARGETS[0]}'])
env.Command(os.path.join(bin_dir, "CIRI.pl"), CIRI[1], SymLink)


# FIND-CIRC
FINDCIRC_tar = 'find_circ.zip'
#FINDCIRC_link = 'http://www.circbase.org/download/' + FINDCIRC_tar
FINDCIRC_link = 'https://github.com/marvin-jens/find_circ/archive/master.zip'
FINDCIRC_target = [os.path.join(tools_dir, FINDCIRC_tar),
                   os.path.join(tools_dir, 'find_circ-master', 'find_circ.py'),
		   os.path.join(tools_dir, 'find_circ-master', 'unmapped2anchors.py')]
#FINDCIRC_source = 
FINDCIRC = env.Command(FINDCIRC_target, [], ['wget -O $TARGET  ' + FINDCIRC_link,
                                             'unzip -o -u -d ' + tools_dir + ' ${TARGETS[0]}'])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), FINDCIRC[1], SymLink)
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), FINDCIRC[2], SymLink)

# gtfToGenePred
gtfToGenePred_link = 'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred'
gtfToGenePred_target = [os.path.join(tools_dir, 'gtfToGenePred')]
gtfToGenePred = env.Command(gtfToGenePred_target, [], 
                            ['wget -O $TARGET ' + gtfToGenePred_link, 
                            Chmod('$TARGET', 0775)])
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), gtfToGenePred, SymLink)

# optparse, ggplot2, DATA.TABLE, plyr, scales, reshape2, ggthemes, RSvgDevice
# Bioconductor: ReportingTools, ballgown
env['ENV']['R_LIBS'] = os.path.join(tools_dir, "R_libs")
#R_libs_targets = [os.path.join(tools_dir, 'R_libs', 'DESeq2', 'R', 'DESeq2')]
R_libs_targets = [os.path.join(tools_dir, 'R_libs', 'data.table', 'R', 'data.table')]
R_libs = env.Command(R_libs_targets, [], 'install_R_libs.R')

# BOWTIE v1
# v1.2.1.1, v1.2.1, and v1.2 do not work!
#BOWTIE1_zip = 'bowtie-1.2.1.1-linux-x86_64.zip'
#BOWTIE1_link = 'https://sourceforge.net/projects/bowtie-bio/files/bowtie/'\
#		'1.2.1.1/' + BOWTIE1_zip
#BOWTIE1_dir = 'bowtie-1.2.1.1'
#BOWTIE1_zip = 'bowtie-1.2.1-linux-x86_64.zip'
#BOWTIE1_link = 'https://sourceforge.net/projects/bowtie-bio/files/bowtie/'\
#               '1.2.1/' + BOWTIE1_zip
#BOWTIE1_dir = 'bowtie-1.2.1'
BOWTIE1_zip = 'bowtie-1.1.2-linux-x86_64.zip'
BOWTIE1_link = 'https://sourceforge.net/projects/bowtie-bio/files/bowtie/'\
               '1.1.2/' + BOWTIE1_zip
BOWTIE1_dir = 'bowtie-1.1.2'

BOWTIE1_target = [os.path.join(tools_dir, BOWTIE1_zip), 
                  os.path.join(tools_dir, BOWTIE1_dir, 'bowtie'),
                  os.path.join(tools_dir, BOWTIE1_dir, 'bowtie-inspect'),
                  os.path.join(tools_dir, BOWTIE1_dir, 'bowtie-build')]
BOWTIE1 = env.Command(BOWTIE1_target, 
                      [], 
                      ['wget -O ${TARGETS[0]} ' + BOWTIE1_link,
                       'unzip -o -u -d ${TARGETS[0].dir} ${TARGETS[0]}']
                     )
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), BOWTIE1[1], SymLink)
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), BOWTIE1[2], SymLink)
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), BOWTIE1[3], SymLink)

# TOPHAT2
tophat2_dir = 'tophat-2.1.0.Linux_x86_64' #'tophat-2.1.1.Linux_x86_64'
tophat2_tar = 'tophat-2.1.0.Linux_x86_64.tar.gz' #'tophat-2.1.1.Linux_x86_64.tar.gz' 
tophat2_url = 'http://ccb.jhu.edu/software/tophat/downloads/' + tophat2_tar

tophat2_target = [os.path.join(tools_dir, tophat2_tar), 
                  os.path.join(tools_dir, tophat2_dir, 'tophat2')]

tophat2 = env.Command(tophat2_target, 
                      [], 
                      ['wget -O ${TARGETS[0]} ' + tophat2_url, 
                       'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}']
                     )
env.Command(os.path.join(bin_dir, "${SOURCE.file}"), tophat2[1], SymLink)

# STRINGTIE
stringtie_dir =	'stringtie-2.1.4.Linux_x86_64'
stringtie_tar = 'stringtie-2.1.4.Linux_x86_64.tar.gz'
stringtie_url = 'http://ccb.jhu.edu/software/stringtie/dl/' + stringtie_tar
stringtie_target = [os.path.join(tools_dir, stringtie_tar),
		    os.path.join(tools_dir, stringtie_dir, 'stringtie')]
stringtie = env.Command(stringtie_target, 
			[], 
			['wget -O ${TARGETS[0]} ' + stringtie_url,
			 'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}']
			)

stringtie_link = env.Command(os.path.join(bin_dir, "${SOURCE.file}"), stringtie[1], SymLink)

# CIRCRNA_FINDER
cfinder_dir = os.path.join(tools_dir, 'circRNA_finder-1.1')
cfinder_tar = 'v1.1.tar.gz'
cfinder_url = 'https://github.com/orzechoj/circRNA_finder/archive/' + cfinder_tar
cfinder_target = [os.path.join(tools_dir, cfinder_tar), 
                  [os.path.join(cfinder_dir, f) for f in
                                                ['postProcessStarAlignment.pl',
                                                  'filterCirc.awk', 
                                                  'filterSpliceSiteCircles.pl', 
                                                  'nrForwardSplicedReads.pl', 
                                                  'starCirclesToBed.pl']]
                 ]
cfinder = env.Command(cfinder_target,
                      [],
                      ['wget -O ${TARGETS[0]} ' + cfinder_url,
                       'tar -xzf ${TARGETS[0]} -C ${TARGETS[0].dir}'])
for t in cfinder[1:]:
    env.Command(os.path.join(bin_dir, "${SOURCE.file}"), t, SymLink)

# SAMTOOLS <v1.0 is required by circrna_finder
samtools0_dir = os.path.join(tools_dir, 'samtools-0.1.20')
samtools0_tar = '0.1.20.tar.gz'
samtools0_url = 'https://github.com/samtools/samtools/archive/' + samtools0_tar
samtools0_target = [os.path.join(tools_dir, samtools0_tar), 
                    os.path.join(samtools0_dir, 'samtools')]
samtools0 = env.Command(samtools0_target, 
                        [], 
                        ['wget -O $TARGET  ' + samtools0_url,
                        'tar -xf ${TARGETS[0]} -C ${TARGETS[0].dir}',
                        ' && '.join(['cd ' + samtools0_dir, 
                                     'make']) + 
                        '; cd ' + Dir('.').abspath])
env.Command(os.path.join(bin_dir, 'samtools_v0', "${SOURCE.file}"), samtools0[1], SymLink)

# GNU PARALLEL
parallel_dir = os.path.join(tools_dir, 'parallel')
parallel_tar = 'parallel-20200922-0.tar.bz2'
parallel_url = 'https://anaconda.org/conda-forge/parallel/20200922/download/linux-64/' + parallel_tar
parallel_target = [os.path.join(tools_dir, parallel_tar), 
                   os.path.join(parallel_dir, 'bin', 'parallel')]
parallel = env.Command(parallel_target,
                      [],
                      ['wget -O ${TARGETS[0]} ' + parallel_url,
                       'tar -xf ${TARGETS[0]} -C ' + parallel_dir])
for t in parallel[1:]:
    env.Command(os.path.join(bin_dir, "${SOURCE.file}"), t, SymLink)


