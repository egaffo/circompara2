'''
This SConscript filters segemehl results using custom script to obtain backsplices.
Formerly, it predicted circular RNAs by means of the
testrealign.x program of segemehl [1]
    
[1] Hoffmann, S. et al. 
    
    A multi-split mapping algorithm for circular RNA, splicing, 
    trans-splicing and fusion detection. 
    
    Genome Biology 15, R34 (2014).

'''

import os

Import('*')

try:
    env = env_testrealign.Clone()
except NameError:
    vars = Variables('vars.py')
    vars.Add('SAMPLE', 'The sample name', 'sample')
    vars.Add('BED', 'The BED file output by Segemehl', 
             '[sample.sngl.bed, sample.trns.txt]')
    
    env = Environment(ENV=os.environ,
                      variables=vars)
    
    Help(vars.GenerateHelpText(env))
    unknown = vars.UnknownVariables()
    if unknown:
        print "Unknown variables:", unknown.keys()
        Exit(1)

## EXTRACT CIRCULAR SPLICE EVENTS, COLLECT SPLICE JUNCTIONS 
## AND FILTER BY MAPPING QUALITY
collect_splice_cmd = "filter_segemehl.R -i ${SOURCES[0]} -t ${SOURCES[1]} $TESTREALIGN_PARAMS "

## let just count fragments and not single read (mates)
#if env['FIX_READ_HEADER']: #env['CCP_COUNT_MODE'] == 'reads'
#    #count separately read mates
#    collect_splice_cmd = collect_splice_cmd + "-m "

collect_splice_cmd = collect_splice_cmd + "-o ${TARGETS[0]} -r ${TARGETS[1]} -l ${TARGETS[2]}"

collect_splice_targets = ["splicesites.bed", 
                          "${SAMPLE}.circular.reads.bed.gz", 
                          "${SAMPLE}.old.segemehl.format.bed"]
collect_splice = env.Command(collect_splice_targets, 
                             env['BED'],
                             collect_splice_cmd)

## COLLECT BACKSPLICED READS
bks_reads_cmd = '''zcat ${SOURCES[0]} | cut -f 4 | sort | uniq -c | '''\
                '''sed -E "s/[^0-9]*([0-9]+)[ ]+([^ ]+)[ ]*/\\1\\t\\2/g" | '''\
                '''sort -k1,1nr > ${TARGETS[0]} '''

bks_reads = env.Command([env['SAMPLE'] + '.testrealign.bks.reads'], 
                        collect_splice_targets[1],
                        bks_reads_cmd)

circ_bed = collect_splice[0]
circ_reads = collect_splice[1]

results = {'CIRCRNAS':    circ_bed,
           #'CIRC_SN_BED': bed,
           'BKS_READS':   bks_reads[0],
           'CIRC_READS':  circ_reads,
           'OLD_BED': collect_splice[2]}

Return('results')

