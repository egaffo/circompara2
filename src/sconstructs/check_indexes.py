'''
This SConscript set which genome indexes have to be built 
for the various read aligners and triggers the required builds,
plus the gene annotation in GenePred format when requested.
It finally returns the index/annotation files and an 
environment with the formatted index parameters.
'''

import os

Import('*')
try:
    env = env_check_indexes.Clone()
except NameError as ne:
    print 'set_indexes.scons: command line execution.'

    vars = Variables('vars.py')
    vars.Add('CIRCRNA_METHODS', '', '')
    vars.Add('CPUS', 'Number of cpus to use for multi thread run', '1')
    vars.Add('GENOME', 'The list of input FASTA files composing the genome sequence. '\
                       'Comma separated.', '')
    vars.Add('HISAT2_EXTRA_PARAMS', 'Extra parameters for htseq2-build', '')
    vars.Add('BWA_EXTRA_PARAMS', 'Extra parameters for bwa index', '')
    vars.Add('BOWTIE2_EXTRA_PARAMS', 'Extra parameters for bowtie2-build', '')
    vars.Add('SEGEMEHL_EXTRA_PARAMS', 'Extra parameters for segemehl.x -x', '')
    vars.Add('STAR_EXTRA_PARAMS', 'Extra parameters for STAR --genomeGenerate', '')
    vars.Add('BOWTIE_EXTRA_PARAMS', 'Extra parameters for bowtie-build', '')
    
    
    env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables", unknown.keys()
        Exit(1)

    env['SCONSCRIPT_HOME'] = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')

annotation_dir = 'indexes'

## BUILD READ ALIGNER PROGRAM GENOME INDEXES IF NOT PROVIDED BY THE USER
hisat2_index_suffixes   = ['.1.ht2', '.2.ht2', '.3.ht2', '.4.ht2',
                          '.5.ht2', '.6.ht2', '.7.ht2', '.8.ht2']
bowtie2_index_suffixes  = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', 
                           '.rev.1.bt2', '.rev.2.bt2']
bowtie_index_suffixes   = ['.1.ebwt', '.2.ebwt', 
                           '.3.ebwt', '.4.ebwt',	
                           '.rev.1.ebwt', '.rev.2.ebwt']
bwa_index_suffixes      = ['.amb', '.ann', '.bwt', '.pac', '.sa']

indexes ={}
genome_indexes_to_build = []

if env['GENOME_INDEX'] == '':
    genome_indexes_to_build.append('HISAT2')
else:
    indexes['HISAT2'] = [env['GENOME_INDEX'] + suffix for suffix in hisat2_index_suffixes]

if not 'circular' in env['BYPASS']:
    ## check for each of the specified methods whether the genome index 
    ## was given as input or it has to be built
    if any([f in env['CIRCRNA_METHODS'] for f in ['testrealign', 
                                                  'circexplorer2_segemehl']]):
        if env['SEGEMEHL_INDEX'] == '':
            # set the index to be built
            genome_indexes_to_build.append('SEGEMEHL')
        else:
            # set the given index to the actual index list
            indexes['SEGEMEHL'] = [env['SEGEMEHL_INDEX']]
    
    if any([f in env['CIRCRNA_METHODS'] for f in ['ciri', 
                                                  'circexplorer2_bwa']]):
        if env['BWA_INDEX'] == '':
            genome_indexes_to_build.append('BWA')
        else:
            indexes['BWA'] = [env['BWA_INDEX'] + suffix for suffix in bwa_index_suffixes]
    
    if 'findcirc' in env['CIRCRNA_METHODS']:
        if env['BOWTIE2_INDEX'] == '':
            genome_indexes_to_build.append('BOWTIE2')
        else:
            indexes['BOWTIE2'] = [env['BOWTIE2_INDEX'] + suffix for suffix in bowtie2_index_suffixes] 
    
    if any([f in env['CIRCRNA_METHODS'] for f in ['circexplorer2_star', 
                                                  'dcc',
                                                  'circrna_finder']]):
        if env['STAR_INDEX'] == '':
    	    genome_indexes_to_build.append('STAR')
        else:
            indexes['STAR'] = [os.path.join(env['STAR_INDEX'], 'SA')]
    
    if 'circexplorer2_tophat' in env['CIRCRNA_METHODS'] or \
    	'--bowtie1' in env['TOPHAT_PARAMS']:
    
        if env['BOWTIE_INDEX'] == '':
    	    genome_indexes_to_build.append('BOWTIE')
        else:
            indexes['BOWTIE'] = [env['BOWTIE_INDEX'] + suffix for suffix in \
    				bowtie_index_suffixes] 
    
if genome_indexes_to_build:

    env_build_indexes = env.Clone()
    env_build_indexes['INDEXES'] = ','.join(genome_indexes_to_build)
    env_build_indexes['CPUS']    = env['CPUS']
    env_build_indexes['GENOME']  = env['GENOME_FASTA']
    env_build_indexes['HISAT2_EXTRA_PARAMS']    = ''
    env_build_indexes['BWA_EXTRA_PARAMS']       = ''
    env_build_indexes['BOWTIE2_EXTRA_PARAMS']   = ''
    env_build_indexes['SEGEMEHL_EXTRA_PARAMS']  = ''
    env_build_indexes['STAR_EXTRA_PARAMS']      = ''
    env_build_indexes['BOWTIE_EXTRA_PARAMS']   = '' #'$( --threads $CPUS $)'

    indexes = SConscript(os.path.join(annotation_dir, 'build_indexes.py'),
                                src_dir = env['SCONSCRIPT_HOME'],
                                variant_dir = annotation_dir, duplicate = 0,
                                exports = '''env_build_indexes ''')

    if env['GENOME_INDEX'] == '':
        #remove the .1.ht2 suffix
        env.Replace(GENOME_INDEX = os.path.abspath(str(indexes['HISAT2'][0]))[0:-6])
    else:
        indexes['HISAT2'] = [env['GENOME_INDEX'] + suffix for suffix in \
							hisat2_index_suffixes]

    if not 'circular' in env['BYPASS']:
        if any([f in env['CIRCRNA_METHODS'] for f in ['testrealign', 'circexplorer2_segemehl']]):
            if env['SEGEMEHL_INDEX'] == '':
    	        env.Replace(SEGEMEHL_INDEX = os.path.abspath(str(indexes['SEGEMEHL'][0])))
            else:
                indexes['SEGEMEHL'] = [env['SEGEMEHL_INDEX']]
    
        if any([f in env['CIRCRNA_METHODS'] for f in ['ciri', 'circexplorer2_bwa']]):
            if env['BWA_INDEX'] == '':
                #remove the .amb suffix
            	env.Replace(BWA_INDEX = os.path.abspath(str(indexes['BWA'][0]))[0:-4])
            else:
                indexes['BWA'] = [env['BWA_INDEX'] + suffix for suffix in bwa_index_suffixes]
    
        if 'findcirc' in env['CIRCRNA_METHODS']:
            if env['BOWTIE2_INDEX'] == '':
                #remove .1.bt2 suffix
                env.Replace(BOWTIE2_INDEX = os.path.abspath(str(indexes['BOWTIE2'][0]))[0:-6])
            else:
                indexes['BOWTIE2'] = [env['BOWTIE2_INDEX'] + suffix for suffix \
    							 in bowtie2_index_suffixes]
                                 
        if any([f in env['CIRCRNA_METHODS'] for f in ['circexplorer2_star', 'dcc', 'circrna_finder']]):
            if env['STAR_INDEX'] == '':
                #index dir
                env.Replace(STAR_INDEX = os.path.dirname(os.path.abspath(str(indexes['STAR'][0]))))
            else:
                indexes['STAR'] = [os.path.join(env['STAR_INDEX'], 'SA')]
    
        if 'circexplorer2_tophat' in env['CIRCRNA_METHODS'] or \
    	   '--bowtie1' in env['TOPHAT_PARAMS']:
            if env['BOWTIE_INDEX'] == '':
                #remove .1.ebwt suffix
                env.Replace(BOWTIE_INDEX = os.path.abspath(str(indexes['BOWTIE'][0]))[0:-7])
            else:
                indexes['BOWTIE'] = [env['BOWTIE_INDEX'] + suffix for suffix in bowtie_index_suffixes]

if not 'circular' in env['BYPASS']:
    ## CREATE genePred FILE FOR CIRCexplorer2 IF ANNOTATION 
    ## WERE PROVIDED ONLY AS GTF
    if any([f in env['CIRCRNA_METHODS'] for f in ['circexplorer2_bwa', 'circexplorer2_star', 
    				                            'circexplorer2_segemehl', 'circexplorer2_tophat']]):
        if env['GENEPRED'] == '': 
            genePred_targets = [os.path.join(annotation_dir, t) for t in ['genePred.transcripts.info',
                                                                   '${SOURCES[0].filebase}.genePred',
                                                                   '${SOURCES[0].filebase}.genePred.wgn']]
            genePred = env.Command(genePred_targets, File(env['ANNOTATION']),
                                   ['gtfToGenePred -infoOut=${TARGETS[0]} ${SOURCES[0]} ${TARGETS[1]}',
                                    'cut -f2 ${TARGETS[0]} | grep -v geneId | '\
                                    'paste - ${TARGETS[1]} | sed "s_^\\t\([^\\t]*\)\\t\(.*\)_\\1\\t\\1\\t\\2_" '\
        			    '> ${TARGETS[2]}']
                                  )
            ## slow command
            #'cut -f1 ${TARGETS[1]} | grep -f - ${TARGETS[0]} | cut -f 9 | '\
            #'paste - ${TARGETS[1]} > ${TARGETS[2]}'
    
            env.Replace(GENEPRED = genePred[2])

results = {'ENV': env, 
           'INDEXES': indexes, 
           'GENEPRED': env['GENEPRED']}

Clean('.', annotation_dir)
Return('results')
