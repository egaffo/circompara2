'''
COMPUTE BACKSPLICES' LINEAR COUNTS FOR CLR (CIRCULAR TO LINEAR RATIO)
'''

import os, re
Import('*')

def writeListFile(target, source, env):
   
    with open(target[0].path, 'w') as listfile:
        listfile.write('\n'.join([s.abspath for s in source]))

    return None


try:
    env = env_circrna_linexp.Clone()
except NameError as ne:
    vars = Variables('vars.py')
    vars.Add('CIRCRNAS', 'A GTF file with circRNA predictions', '')
    vars.Add('RUNS_DICT', '', '')


    env = Environment(ENV=os.environ, SHELL = '/bin/bash',
                      variables=vars)
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print("Run sample: unknown variables", list(unknown.keys()))
        Exit(1)

env['CIRCOMPARA_HOME'] = env['ENV']['CIRCOMPARA_HOME']

env.SetDefault(LIN_COUNTER = 'dcc')

strandness_pattern = re.compile("--rna-strandness\s+[FR]{1,2}")

if 'dcc' in env['LIN_COUNTER']:
    ## convert circRNA GTF to BED
    circ_coords_cmd = '''zcat ${SOURCES[0]} | '''\
                      '''sed -r 's/([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
                                '''([^\\t]+)\\t([^\\t]+)\\t([^\\t]+)\\t'''\
                                '''([^\\t]+)\\t([^\\t]+)\\tgene_id "([^"]+)";'''\
                                '''/echo -e "\\1\\t\\$$((\\4-1))\\t\\5\\t'''\
                                          '''@\\9@\\t\\6\\t\\7"/e' | '''\
                      '''sed -r 's/@/"/g' > $TARGET && '''\
                      '''sed -i '1ichr\\tstart\\tend\\tname\\tscore'''\
                      '''\\tstrand' $TARGET'''

    circ_coords = env.Command('circrnas.bed', 
                              env['CIRCRNAS'], 
                              circ_coords_cmd)

    ## use DCC to compute circRNA host gene linear expression
    btmcovstrand = '-N' ## check if stranded read alignment
    if strandness_pattern.search(env['HISAT2_EXTRA_PARAMS']):
        btmcovstrand = ''
            
    ## make file listing circRNA coordinates for each sample:
    ## we actually repeat the same circRNA coordinates file
    ## as it is valid for all samples
    env['N_SAMPLES'] = len(list(env['RUNS_DICT'].keys()))
    circrna_beds = env.Command('circrna_beds', 
                              [circ_coords for i in range(0, env['N_SAMPLES'])], 
                              writeListFile)
    
    ## make file listing linear mapping files for the samples
    bams = env.Command('bam_files', 
                       [env['RUNS_DICT'][s]['LINEAR_ALIGNMENTS'] for 
                        s in sorted(env['RUNS_DICT'].keys())],
                       writeListFile)
    
    ## compute linear counts and set the right names in table header
    dcclinexp_sources = [bams, circrna_beds, circ_coords]
    dcclinexp_target = 'LinearCount'
    dcclinexp_cmd = '''DCC -G ''' + btmcovstrand + \
                 ''' -A $GENOME_FASTA -O $TARGET.dir '''\
                 '''-B @${SOURCES[0].abspath} -C ${SOURCES[2].abspath} '''\
                 '''@${SOURCES[1].abspath} && '''\
                 '''sed -i '1s/.*/chr\\tstart\\tend\\t''' +\
                 '\\t'.join([s for s in sorted(env['RUNS_DICT'].keys())]) +\
                 '''/' $TARGET'''
    dcclinexp = env.Command(dcclinexp_target, 
                         dcclinexp_sources, 
                         dcclinexp_cmd)

    linexp = dcclinexp
    
if 'ccp' in env['LIN_COUNTER']:
    
    genome_file_cmd = '''fasta_len.py $SOURCE > $TARGET'''
    genome_file = env.Command('chr_len.genome', 
                              env['GENOME_FASTA'],
                              genome_file_cmd)
    ## TODO: is it safer to instead use the BAM header to get the chromosome sorting?
    ## Should it be done on the fly for each BAM file? This would imply to (unnecessarily?) 
    ## re-sort the circRNA file for each sample
    ## genome_file_cmd = '''samtools view -H $SOURCE | grep "^@SQ" | '''\
    ##                   '''sed -r 's/.*SN:([^\t]+)\tLN:([0-9]+)/\1\t\2/' > $TARGET'''
    ## genome_file = env.Command('chr_len.genome', env['RUNS_DICT'][s]['LINEAR_ALIGNMENTS']], genome_file_cmd)
    
    btmcovstrand = '' ## check if stranded read alignment
    strandness_pattern = re.compile("--rna-strandness\s+[FR]{1,2}")
    if strandness_pattern.search(env['HISAT2_EXTRA_PARAMS']):
        btmcovstrand = '-s'
    
    ## TODO: improve the counting by considering only linear spliced reads
    ## and linear reads spanning the splice site, while not counting reads
    ## exactly matching splice site bound (for which we cannot decide 
    ## whether they belong to the circular or the linear transcript)
    #for s in sorted(runs_dict.keys()):
    #    ## keep only spliced reads
    #    target_1 = "lin_spliced_reads_in_bks.tab"
    #    sources_1 = [runs_dict[s]['LINEAR_ALIGNMENTS'], 
    #                 circexp['SNP_UNIQUE_CIRC']]
    #    cmd_1 = '''samtools view ${SOURCES[0]} | grep "N" | bedtools '''\
    #            '''coverage -counts -sorted ''' + btmcovstrand + \
    #            ''' -a ${SOURCES[1]} -b stdin > $TARGET '''
    #    lin_spliced_reads_in_bks = env.Command(target_1, sources_1, cmd_1)
    
    #    ## keep only linear reads spanning backsplice
    #    ## Hint: use bedtools slop to consider bases outside the 
    #    ## backsplice start/end sites
    #    target_2 = "lin_reads_spanning_bks.tab"
    #    sources_2 = sources_1.append(genome_file)
    #    cmd_2 = '''samtools view ${SOURCES[0]} | grep -v "N" | '''\
    #            '''bedtools coverage -counts -sorted ''' + btmcovstrand + \
    #            ''' -a <( bedtools slop -b 1 -i ${SOURCES[1]} ) '''\
    #            '''-b stdin > $TARGET '''
    #    lin_reads_spanning_bks = env.Command(target_2, sources_2, cmd_2)
    
    ## TODO: preliminary select read aligned to circRNAs
    # bedtools intesect sn_unique_circ.gtf f1.bam f2.bam
     
    ## count alignments (also) on circRNA outer positions.
    ## The -s in bedtools flank is not necessary since the flanking 
    ## regions are the same size in this case. Yet, it maybe useful 
    ## for chromosome terminal postions....   
    #btmcov_sources = [[env['CIRCRNAS'], genome_file], 
    #                  [env['RUNS_DICT'][s]['LINEAR_ALIGNMENTS'] for 
    #                             s in sorted(env['RUNS_DICT'].keys())]]
    #btmcov_target = 'bks_linear_counts.tab'
    #btmcov_cmd = '''bedtools flank -i ${SOURCES[0]} -g ${SOURCES[1]} -s -b 1 | '''\
    #             '''bedtools multicov ''' + btmcovstrand + \
    #             ''' -bed stdin -bams ${SOURCES[2:]} > ${TARGET} '''\
    #             '''&& sed -i '1ichr\\tsource\\tfeature\\tstart\\tend'''\
    #             '''\\tscore\\tstrand\\tframe\\tname\\t''' +\
    #             '\\t'.join([s for s in sorted(env['RUNS_DICT'].keys())]) +\
    #             '''' $TARGET'''
    #btmcov = env.Command(btmcov_target, 
    #                     btmcov_sources, 
    #                     btmcov_cmd)
    #linexp = btmcov

    sorted_sn_cmd = '''bedtools flank -i ${SOURCES[0]} -g ${SOURCES[1]} -s -b 1 | '''\
                    '''bedtools sort -faidx ${SOURCES[1]} -i stdin > $TARGET'''
                    #'''sort -k1,1 -k4,4n > $TARGET'''
    sorted_sn = env.Command('sorted_sn_circ.gtf',
                            [env['CIRCRNAS'], genome_file],
                            sorted_sn_cmd)
    res = {}
    for s in sorted(env['RUNS_DICT'].keys()):
        ## TODO: preliminary select read aligned to circRNAs
        # bedtools intesect env['CIRCRNAS'] f1.bam | bedtools coverage
        btmcov_sources = [sorted_sn, genome_file,
                          env['RUNS_DICT'][s]['LINEAR_ALIGNMENTS']]
        btmcov_target = s + '_bks_linear_counts.tab'
        #btmcov_cmd = '''bedtools coverage -counts -sorted ''' + btmcovstrand + \
        #             ''' -g ${SOURCES[1]} -a ${SOURCES[0]} -b ${SOURCES[2]} > ${TARGET} '''
        btmcov_cmd = '''parallel --pipepart -k -j $CPUS -a ${SOURCES[0]} --block -10 "''' + \
                     '''bedtools coverage -counts -sorted ''' + btmcovstrand + \
                     ''' -g ${SOURCES[1]} -a stdin -b ${SOURCES[2]}" > ${TARGET} '''
        res[s] = env.Command(btmcov_target, 
                             btmcov_sources, 
                             btmcov_cmd)

    linexp = env.Command('bks_linear_counts.tab',
                         [list(res.values())], 
                         "grep -H '.' ${SOURCES} > $TARGET")
   
Return('linexp')
