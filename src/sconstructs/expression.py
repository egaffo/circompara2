'''
This script functions as a dispatcher to the desired gene expression analysis methods.
It sets the necessary variables and launch the expression analyses.

Variables to export when calling from a SConscript:
 * expression_cpus
 * expression_annotation
 * mapping_file
 * sample_name
 * expression_cufflinks_params

'''

import os

Import('*')

try:
    env = env_expression.Clone()
except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('CPUS', 'Set number of CPUs', '4')
    vars.Add('ANNOTATION', 'The GFF/GTF file with gene annotation (e.g. from Ensembl)', 'exons.gtf') 
    vars.Add('ALIGNMENTS', 'The read alignment file in SAM/BAM format', 'sample.bam')
    vars.Add('SAMPLE', 'The sample name', 'sample')
    vars.Add('TOGGLE_TRANSCRIPTOME_RECONSTRUCTION', 'Set True to enable transcriptome '\
	         'reconstruction. Default only quantifies genes and transcripts from the given '\
	         'annotation GTF file', 
             'False')
    vars.Add('METHODS', 'Name or comma separated list of the method(s) to use. '\
                        'Currently suppported: stringtie',
             'stringtie')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

    CPUS         = env['CPUS']
    ANNOTATION   = env['ANNOTATION']
    ALIGNMENTS   = env['ALIGNMENTS']
    SAMPLE       = env['SAMPLE']
    #CUFFLINKS_PARAMS = env['CUFFLINKS_PARAMS']
    TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = False
    if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION'] == 'True':
    	TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = True
    
    env.Replace(METHODS = env['METHODS'].split(','))

expression_results = []
#expression = {'HTSEQ': None, 'CUFFLINKS': None, 'STRINGTIE': None}
expression = {'STRINGTIE': None}

SRC_DIR = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')

#mapping_file = ALIGNMENTS
#sample_name  = SAMPLE

mapping_file = env['ALIGNMENTS']
sample_name  = env['SAMPLE']

#htseq_counts_dir = 'htseq'
#cufflinks_dir = 'cufflinks'
stringtie_dir = 'stringtie'

#if 'htseq' in env['LINEAR_EXPRESSION_METHODS'] and not env['ANNOTATION'] == '':
#    ## RUN HTSEQ-COUNT
#    env_htseq_count = env.Clone()
#    if '--rna-strandness FR' in env['HISAT2_EXTRA_PARAMS']: 
#        env_htseq_count['STRANDED'] = 'yes'
#    else:
#        env_htseq_count['STRANDED'] = 'no'
#
#    ## GET READ COUNTS
#    htseq_count = SConscript(os.path.join(htseq_counts_dir, 'htseq_count.py'), 
#                             src_dir = env['SCONSCRIPT_HOME'], 
#                             variant_dir = htseq_counts_dir, duplicate = 0, 
#                             exports = '''env_htseq_count''')
#    expression_results.append(htseq_count)
#    expression['HTSEQ'] = htseq_count
#
#if 'cufflinks' in env['LINEAR_EXPRESSION_METHODS']:
#    ## RUN CUFFLINKS
#    cufflinks_annotation = env['ANNOTATION']
#    cufflinks_cpus = env['CPUS']
#    cufflinks_params = env['CUFFLINKS_PARAMS']
#    cufflinks_toggle_transcriptome_reconstruction = env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION']
#    cufflinks = SConscript(os.path.join(cufflinks_dir, 'cufflinks.py'),
#                           variant_dir = cufflinks_dir, src_dir = SRC_DIR, 
#                           duplicate = 0, exports = '''env mapping_file sample_name '''
#                           '''cufflinks_annotation cufflinks_cpus cufflinks_params '''
#			   '''cufflinks_toggle_transcriptome_reconstruction''')
#    expression_results.append(cufflinks)
#    expression['CUFFLINKS'] = cufflinks

if 'stringtie' in env['LINEAR_EXPRESSION_METHODS']:
    ## RUN STRINGTIE
    env_stringtie = env.Clone()
    stringtie = SConscript(os.path.join(stringtie_dir, 'stringtie.py'), 
                            variant_dir = stringtie_dir, src_dir = SRC_DIR,
                            duplicate = 0, exports = 'env_stringtie')

    expression_results.append(stringtie.values())
    expression['STRINGTIE'] = stringtie

    ## compute read raw counts for StringTie
    raw_counts_sources = [stringtie['TRANSCRIPTS_GTF'],
                          env['FASTQC_DATA']]
    raw_counts_cmd = os.path.join('''get_stringtie_rawcounts.R -g ${SOURCES[0]} '''\
                             '''-f ${','.join([str(s.abspath) for s in SOURCES[1:]])} '''\
                             '''-o ${TARGETS[0].dir}''', 
                             env['SAMPLE'] + '''_''')
    raw_counts_targets = [os.path.join(stringtie_dir, env['SAMPLE'] + "_" + f) \
                          for f in ['gene_expression_rawcounts.csv',
                                    'transcript_expression_rawcounts.csv']]
    raw_counts = env.Command(raw_counts_targets, 
                             raw_counts_sources, 
                             raw_counts_cmd)

    expression_results.append(raw_counts)
    expression['STRINGTIE']['RAW_COUNTS'] = raw_counts


#Clean('.', htseq_counts_dir)
#Clean('.', cufflinks_dir)
Clean('.', stringtie_dir)

Return('expression_results expression')

