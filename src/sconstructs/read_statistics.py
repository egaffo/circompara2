'''
This is a SConscript script that performs  FASTQC

Software dependencies:
 * FASTQC

When called from a SConscript it imports the following variables:
 * env
 * read_statistics_readset
  
'''

import os

Import('*')

try:
    # these are the variables passed with 'exports' from a calling SConscript
    env = env_read_statistics.Clone()
    readset = env['READSET']#read_statistics_readset
except NameError, ne:
    vars = Variables('vars.py')
    vars.Add('READS', 'FASTQ read file to process, either in plain text or gzipped (.gz)', 
             'reads.fq')
    vars.Add('READSTAT_METHODS', 'Comma separated list of methods to use for read statistics. '\
             'Currently supported: fastqc', 'fastqc')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)
    # These are the variables given from the command line when the SConscript is called
    # standalone
    readset = env['READS']

SRC_DIR = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')

env.SetDefault(READSTAT_METHODS = 'fastqc')

readset_basename = os.path.splitext(os.path.basename(readset))[0]
readset_ext = os.path.splitext(os.path.basename(readset))[1]

if readset_ext == '.gz':
    pre_cmd = 'zcat $SOURCE | '
else:
    pre_cmd = 'cat $SOURCE | '

stats = []
## COMPUTE STATISTICS
fastqc = []
if('fastqc' in env['READSTAT_METHODS'].split(',')):
    ## perform FASTQC analyses (default)
    fastqc_dir   = 'fastqc_stats'
    env_fastqc = env.Clone()
    #fastqc_readset = File(readset)
    env_fastqc['READSET'] = File(readset)
    fastqc = env.SConscript(os.path.join(fastqc_dir, 'fastqc.py'), 
                                  src_dir = SRC_DIR,
                                  variant_dir = fastqc_dir, duplicate = 0,
                                  exports = 'env_fastqc')
    Clean('.', fastqc_dir)

Return('stats fastqc')
