'''
This SConscript impleemnts the building of a genome index for the SEGEMEHL read aligner.
Current implementation builds one single index file.

Import:
    * env_index_segemehl environment including CPUS, GENOME, and EXTRA_PARAMS variables.
                       The GENOME environment variable must report absolute paths of 
                       the FASTA files
Returns:
    * index.idx
'''

import os

Import('*')

try:
    env = env_index_segemehl
except NameError as ne:
    print 'index_segemehl.scons: command line execution'
    vars = Variables('vars.py')
    vars.Add('GENOME', 'The list of input FASTA files composing the genome sequence. '\
            'Comma separated', '')
    vars.Add('EXTRA_PARAMS', 'Extra parameters for segemehl', '')
    
    env = Environment(ENV=os.environ, SHELL = '/bin/bash', variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables", unknown.keys()
        Exit(1)
    
SCONSCRIPT_HOME = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')

GENOME      = env['GENOME'].split(',')
EXTRA_PARAMS= env['EXTRA_PARAMS']

source  = [File(f).abspath for f in GENOME]
target_basename = '_'.join([os.path.splitext(os.path.basename(f))[0] for f in GENOME])
index_suffixes = ['.idx']
target  = [target_basename + suffix for suffix in index_suffixes]
command = '''segemehl.x -d ''' + ' '.join(source) +\
          ''' -x ${TARGET}''' + EXTRA_PARAMS
index = env.Command(target, source, command)

Return('index')
