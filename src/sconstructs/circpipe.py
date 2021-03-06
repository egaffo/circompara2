'''
This SConscript characterize the transcriptome from one RNA-seq sample.
It composes a pipeline that gives statistics on the raw data and
also processes the raw reads to map them on a reference genome, 
reconstruct transcripts, estimates gene/transcript expression, and 
also detects circRNAs by means of backsplice junctions.

Software dependencies are inherited from the CIRCOMPARA-SConscripts used:
 * read_statistics.py
 * preprocess.py
 * hisat2.py
 * expression.py
 * circrna_methods.py

When called from a SConscript it imports the following variables:
 * env_circpipe

'''

import os, itertools, re

## SET (DISPATCHER) SCRIPT NAMES. THESE ARE THE MAIN PIPELINE STEPS, EACH OF THEM MIGHT FIRE
## DIFFERENT ACTIONS/SCONSCRIPTS
read_statistics = 'read_statistics.py'
preprocess      = 'preprocess.py'
mapping         = 'hisat2.py'
expression      = 'expression.py'
circrna_methods = 'circrna_methods.py'

Import('*')

try:
    env = env_circpipe.Clone()
except NameError:
    varfile = ARGUMENTS.get('VARS', 'vars.py')
    vars = Variables(varfile)
    vars.Add('CPUS', 'Set number of CPUs', '4')
    vars.Add('PREPROCESSOR', 'The preprocessing method', 'trimmomatic')
    vars.Add('PREPROCESSOR_PARAMS', 
             '''Read preprocessor extra parameters. F.i. if Trimmomatic, an empty string '''\
             '''defaults to '''\
             '''MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:50 AVGQUAL:30 ''', 
             '')
    vars.Add('ADAPTER_FILE', 'FASTA file full path of the adapter sequence', '')
    vars.Add('ANNOTATION', 'Gene annotation (Ensembl GFF)', '')
    vars.Add('GENOME_FASTA', 'The FASTA file with the reference genome', 'genome.fa')
    vars.Add('GENOME_INDEX', '''The index of the reference genome for HISAT2''', '/path/to/index')
    vars.Add('READS', 'RNA-seq reads. Comma separated list if paired-end', 'reads.fa')

    vars.Add('SEGEMEHL_INDEX', '''The .idx index for segemehl''', 'genome.idx')
    vars.Add('BWA_INDEX', '''The index of the reference genome for BWA''','/path/to/index')
    vars.Add('BWA_PARAMS','Extra parameters for BWA','')
    #vars.Add('CIRI', 'The full path to the CIRI_vx.x.pl perl script', '')
    vars.Add('BOWTIE2_INDEX', '''The index of the reference genome for BOWTIE2''','/path/to/index')

    vars.Add('STAR_INDEX', 'The directory path where to find Star genome index', 
             '/path/to/Star/index/dir')
    vars.Add('GENEPRED', 'The genome annotation in GenePred format', 'genes.genePred')
    vars.Add('HISAT2_EXTRA_PARAMS', '''Extra parameters to add to the HISAT2 aligner fixed '''\
             '''parameters '--dta --dta-cufflinks --rg-id <SAMPLE> --no-discordant '''\
             '''--no-mixed'. For instance, '--rna-strandness FR' if stranded reads'''\
             ''' are used.''', '')
    vars.Add('CUFFLINKS_PARAMS', '''Cufflinks extra parameters. '''\
             '''F.i. '--library-type fr-firststrand' if dUTPs stranded library were used '''\
             '''for the sequencing''', '')
    vars.Add('CIRI_EXTRA_PARAMS', 'CIRI additional parameters', '')
    vars.Add('CIRCRNA_METHODS', 'Comma separated list of circRNA detection methods to use. '\
	     'Use all methods available as default', '')
    vars.Add('TOGGLE_TRANSCRIPTOME_RECONSTRUCTION', 'Set True to enable transcriptome '\
	     'reconstruction. Default only quantifies genes and transcripts from the given '\
	     'annotation GTF file', 'False')
    vars.Add('BYPASS', 'Skip analysis of linear/circular transcripts. This will also skip '\
	     'the analysis of linear-to-circular expression correlation.'\
         '{linear,circular}',
	     'False')

    env = Environment(variables = vars,
                      ENV = os.environ)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Run sample: unknown variables", unknown.keys()
        Exit(1)

    if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION'].tolower() == 'true':
    	env.Replace(TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = True)
    else:
    	env.Replace(TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = False)

    env.Replace(READS = env['READS'].split(','))
    env.Replace(CIRCRNA_METHODS = env['CIRCRNA_METHODS'].strip().split(','))

SRC_DIR = os.path.join(env['ENV']['CIRCOMPARA_HOME'], 'src', 'sconstructs')

results = []

## GET INPUT READ FILE FULL PATH
raw_reads = [File(f).abspath for f in env['READS']]

## COMPUTE STATISTICS ON RAW READS
read_statistics_dir = 'read_statistics'
raw_read_stats = []
fastqc_data = []
for readset in raw_reads:
    env_read_statistics = env.Clone()
    #read_statistics_readset = readset
    env_read_statistics['READSET'] = readset
    stat = SConscript(os.path.join(read_statistics_dir, read_statistics), 
                      src_dir=SRC_DIR, variant_dir=read_statistics_dir, duplicate=0,
                      exports='env_read_statistics')
    raw_read_stats.append(stat)
    ## explicit list of fastqc_data.txt files
    if('fastqc' in env['READSTAT_METHODS'].split(',')):
        fastqc_data.append(stat[1][4])

results.append(raw_read_stats)

build_dir = 'processings'
## EXECUTE SCRIPTS: 
## READ PREPROCESSING
env_preprocess = env.Clone()
preprocess  = env.SConscript(os.path.join(build_dir, preprocess), 
                             variant_dir = build_dir, src_dir = SRC_DIR, 
                             duplicate = 0, exports='env_preprocess')

#results.append([preprocess['READS'], preprocess['STATS']]) ##TODO:fix:adding stats makes stats dependency on genome index files
results.append([preprocess['READS']])
## set the clean read file paths, relative to the calling SConstruct (not the calling SConscript).
## Implementation note: here you should know what the preprocess SConscript returns to handle it 
## properly. For instance: Trimmomatic with paired-end reads gives four output read files and 
## you probably want to use the 0 and 2 indexed files of the returned list (i.e. the processed 
## reads with mate read). Instead, if no preprocessing was performed, the clean reads are the raw
## reads and the list is only two elements (for paired-end reads)

## get first-in-pair or the single-end preprocessed read file
clean_reads = [preprocess['READS'][0][0]]

## case: Trimmomatic paired-end
if env['PREPROCESSOR'] == 'trimmomatic':
    ## paired end case: add the second-in-pair clean read file
    if len(preprocess['READS'][0]) > 2:
        clean_reads.append(preprocess['READS'][0][2])
        fastqc_data = []
        if('fastqc' in env['READSTAT_METHODS'].split(',')):
            ## get clean reads' statistics only
            for cr in clean_reads:
                read_filename = re.sub("\.gz", "", os.path.basename(cr.abspath))
                fastqc_data.append(get_matching_nodes(preprocess['STATS'], 
                                   os.path.join('.*' + read_filename + '_fastqc',
                                   'fastqc_data\.txt'))[0])
                                   
## case: no preprocessing was performed
elif env['PREPROCESSOR'] == '':
    if len(preprocess['READS'][0]) > 1:
        clean_reads.append(preprocess['READS'][0][1])
else:
    clean_reads.append(preprocess['READS'][0][1])

mappings = {'BAM': None,
            'BAM_INDEX': None}

if not 'linmap' in env['BYPASS']:
    ## ALIGN TO GENOME
    env_hisat2 = env.Clone()
    env_hisat2['HISAT2_INDEX'] = env['GENOME_INDEX']
    env_hisat2['HISAT2_PARAMS'] = env['HISAT2_EXTRA_PARAMS'].split() if isinstance(env['HISAT2_EXTRA_PARAMS'], str) else env['HISAT2_EXTRA_PARAMS']
    env_hisat2.AppendUnique(HISAT2_PARAMS = ['--dta', '--dta-cufflinks', '--rg-id',
                                             env['SAMPLE'], '--no-discordant', 
                                             '--no-mixed'])
    env_hisat2['READS'] = clean_reads
    mappings    = env.SConscript(os.path.join(build_dir, mapping), 
                                 variant_dir = build_dir, src_dir = SRC_DIR, 
                                 duplicate = 0, exports = '''env_hisat2''')
    results.append([File(f) for f in mappings.items()])
    env.Replace(LINMAPS = mappings['BAM'])
    
elif env['LINMAPS']:    
    ## use pre-computed linear alignments, if given
    mappings = {'BAM': env['LINMAPS'],
                'BAM_INDEX': env['LINMAPS'].abspath + '.bai'}

env['ALIGNMENTS'] = mappings['BAM'] 

if env['ALIGNMENTS'] and (not 'linear' in env['BYPASS']):
    ## EXPRESSION ANALYSIS
    env_expression = env.Clone()
    env_expression['FASTQC_DATA'] = fastqc_data
    
    expression  = env.SConscript(os.path.join(build_dir, expression), 
                                 variant_dir = build_dir, src_dir = SRC_DIR, 
                                 duplicate = 0, exports = '''env_expression''')

    if 'stringtie' in env['LINEAR_EXPRESSION_METHODS']:
        results.append(expression[1]['STRINGTIE']['ALL_TARGETS'])
        results.append(expression[1]['STRINGTIE']['RAW_COUNTS'])
        # the ALL_TARGETS is a temporary hack. It will be removed when SConscript
        # results are passed as dictionaries/custom-objects and handled accordingly
        if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION']:
            env.Replace(ANNOTATION = expression[1]['STRINGTIE']['TRANSCRIPTS_GTF'].abspath)

    #if 'htseq' in env['LINEAR_EXPRESSION_METHODS']:
    #    results.append(expression[1]['HTSEQ'])

    #if 'cufflinks' in env['LINEAR_EXPRESSION_METHODS']:
    #    results.append(expression[1]['CUFFLINKS'])
    #    if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION']:
    #        env.Replace(ANNOTATION = expression[1]['CUFFLINKS'][0])

    ## if transcriptome was reconstructed, then update the genePred
    ## annoation file required by CIRCexplorer2
    if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION']:
        env_check_indexes = env.Clone()
        env_check_indexes['GENEPRED'] = ''
        genepred = env.SConscript(os.path.join(build_dir,
                                               'check_indexes.py'),
                                  variant_dir = build_dir, 
                                  src_dir = SRC_DIR,
                                  duplicate = 0,
                                  exports = '''env_check_indexes''')
        env.Replace(GENEPRED = genepred['GENEPRED'])

else:
	expression = [None, 
                  env.Command('mock_expression.txt', 
                              env['ALIGNMENTS'], 
                              'touch $TARGET')]

## CIRCRNA DETECTION
unmapped_dir = os.path.join(build_dir, 'unmapped_reads')

circrnas = [None, None]
if not 'linmap' in env['BYPASS']:
    ## Extract unmapped reads in FASTQ format from BAM alignments
    unmapped_reads_target = ['singleton.fastq.gz']
    #unmapped_reads_cmd = '''samtools fastq -f 12 -F 3328 -n -s >( gzip -c > ${TARGETS[0]} ) '''
    unmapped_reads_cmd = '''samtools fastq -f 12 -F 3328 -n -s ${str(TARGETS[0].abspath).replace('.gz', '')}  '''
    #SAMFLAG 3328 = 256+1024+2048
    #4: read unmapped
    #8: mate unmapped
    #256: not primary alignment
    #1024: read is PCR or optical duplicate
    #2048: supplementary alignment
    # samtools fastq:
    # -s FILE   write singleton reads to FILE [assume single-end]
    # -n        don't append /1 and /2 to the read name
    if len(env['READS']) > 1:
        paired_end_read_files = ['unmapped_1.fastq.gz', 'unmapped_2.fastq.gz']
        unmapped_reads_target.extend(paired_end_read_files)
        unmapped_reads_cmd = unmapped_reads_cmd +\
                             ''' -1 ${str(TARGETS[1].abspath).replace('.gz', '')} '''\
                             ''' -2 ${str(TARGETS[2].abspath).replace('.gz', '')} ${SOURCES[0]} '''\
                             '''&& gzip ${str(TARGETS[1].abspath).replace('.gz', '')} '''\
                             '''&& gzip ${str(TARGETS[2].abspath).replace('.gz', '')} '''
        unmapped_reads_cmd = unmapped_reads_cmd +\
                            ''' && gzip ${str(TARGETS[0].abspath).replace('.gz', '')}'''
        
        unmapped_reads = env.Command([os.path.join(unmapped_dir, f) for f in unmapped_reads_target], 
                                     [env['ALIGNMENTS'].abspath,
                                      mappings['BAM_INDEX']], 
                                     unmapped_reads_cmd)
    else:
        unmapped_reads = mappings['UNMAPPED_READS']
    
    results.append(unmapped_reads)

else:
    if len(env['READS']) > 1:
        unmapped_reads = ['mockfile.fq.gz'] + env['READS']
    else:
        unmapped_reads = env['READS']

if not 'circular' in env['BYPASS']:
    # run circRNA detection pipeline 
    env_sample_circrna_methods = env.Clone()

    if len(env['READS']) > 1:
        env_sample_circrna_methods.Replace(READS = [f.abspath for f in unmapped_reads[1:3]])
    else:
        env_sample_circrna_methods.Replace(READS = [unmapped_reads.abspath])
    
    ## set strandness parameter for circRNA detection methods
    ## check by HISAT parameters
    strandness_pattern = re.compile("--rna-strandness\s+[FR]{1,2}")
    if strandness_pattern.search(env['HISAT2_EXTRA_PARAMS']):
        ## Environment is propagated to the Sconscript
        env_sample_circrna_methods.AppendUnique(FINDCIRC_EXTRA_PARAMS = ['--strandpref', '--stranded'])
        # env_sample_circrna_methods.AppendUnique(DCC_EXTRA_PARAMS = ['-ss']) ## default DCC is stranded
        
        ## set automatically and accordingly the other methods' strandness parameters
        if 'RF' in env['HISAT2_EXTRA_PARAMS'].split():
            env_sample_circrna_methods.AppendUnique(TOPHAT_PARAMS = ['--library-type', 'fr-firststrand'])
        elif 'FR' in env['HISAT2_EXTRA_PARAMS'].split():
            env_sample_circrna_methods.AppendUnique(TOPHAT_PARAMS = ['--library-type', 'fr-secondstrand'])
    
    circrnas = env.SConscript(os.path.join(build_dir, circrna_methods), 
                            src_dir = SRC_DIR, 
                            variant_dir = build_dir, duplicate = 0,
                            exports = '''env_sample_circrna_methods''')
    results.append(circrnas[0]) ## the list result format

    if not 'linear' in env['BYPASS']:
        if env['TOGGLE_TRANSCRIPTOME_RECONSTRUCTION']:
            ## could be improved by setting only the reconstructed 
            ## transcriptome file as the dependency, f.i.
            ## expression[1]['STRINGTIE']['TRANSCRIPTS_GTF'] if stringtie
            Depends(circrnas[0], expression[0]) 

results_dict = {'RAW_READ_STATS'        : raw_read_stats,
                'PREPROCESSING'         : preprocess,
                'LINEAR_EXPRESSION'     : expression[1],
                'CIRCULAR_EXPRESSION'   : circrnas[1], ## the dict results format
                'LINEAR_ALIGNMENTS'     : env['ALIGNMENTS'],
                'LINEAR_ALIGNMENTS_IDX' : mappings['BAM_INDEX'],
                'LINEAR_MAPPING_RES'    : mappings,
                'FASTQC_DATA'           : fastqc_data}

## CLEAN COMMANDS
Clean('.', read_statistics_dir)
Clean('.', build_dir)
Clean('.', unmapped_dir)

Return('results results_dict')
