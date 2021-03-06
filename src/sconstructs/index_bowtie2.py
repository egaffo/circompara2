'''
This SConscript impleemnts the building of a genome index for the Bowtie2 read aligner
Import:
    * env_index_bowtie2 environment including CPUS, GENOME, and EXTRA_PARAMS variables.
                       The GENOME environment variable must report absolute paths of 
                       the FASTA files
Returns:
    * ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
'''

import os

def SymLink(target, source, env):
    os.symlink(os.path.abspath(str(source[0])), os.path.abspath(str(target[0])))

Import('*')

try:
    env = env_index_bowtie2
except NameError as ne:
    print 'index_bowtie2.scons: command line execution'
    vars = Variables('vars.py')
    vars.Add('CPUS', 'Number of cpus to use for multi thread run', '1')
    vars.Add('GENOME', 'The list of input FASTA files composing the genome sequence. '\
            'Comma separated', '')
    vars.Add('EXTRA_PARAMS', 'Extra parameters for bowtie2-build', '')
    
    env = Environment(ENV=os.environ, SHELL = '/bin/bash', variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables", unknown.keys()
        Exit(1)
    
SCONSCRIPT_HOME = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')
    
CPUS        = env['CPUS']
GENOME      = env['GENOME'].split(',')
EXTRA_PARAMS= env['EXTRA_PARAMS']

source  = [File(f).abspath for f in GENOME]
target_basename = '_'.join([os.path.splitext(os.path.basename(f))[0] for f in GENOME])
index_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
target  = [target_basename + suffix for suffix in index_suffixes]

link_genome_fasta = env.Command(os.path.join(str(File(target[0]).dir), '${SOURCES[0].file}'), 
                                File(GENOME), 
                                SymLink)

command = '''bowtie2-build -f --seed 1 $(--threads ''' + CPUS +\
          ''' $) '''  + EXTRA_PARAMS + ' ' + ','.join(source) +\
          ''' ${TARGETS[0].dir}''' + os.path.sep + target_basename 
index = env.Command(target, [source, link_genome_fasta], command)

Return('index')
