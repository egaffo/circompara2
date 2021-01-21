META            = 'meta.csv'
GENOME_FASTA    = '../annotation/CFLAR_HIPK3.fa'
ANNOTATION      = '../annotation/CFLAR_HIPK3.gtf' 
CPUS            = '4'

## pre-computed index and annotation files 
#GENOME_INDEX = "../indexes/hisat2/CFLAR_HIPK3"
#SEGEMEHL_INDEX = "../indexes/segemehl/CFLAR_HIPK3.idx"
#BWA_INDEX = "../indexes/bwa/CFLAR_HIPK3"
#BOWTIE2_INDEX = "../indexes/bowtie2/CFLAR_HIPK3"
#BOWTIE_INDEX = "../indexes/bowtie/CFLAR_HIPK3"
#STAR_INDEX = "../indexes/star/CFLAR_HIPK3"
#GENEPRED = "../annotation/CFLAR_HIPK3.genePred.wgn"

#PREPROCESSOR    = 'trimmomatic'
#PREPROCESSOR_PARAMS = 'MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:35 AVGQUAL:30'
#FIX_READ_HEADER = 'True'
#HISAT2_EXTRA_PARAMS = '--rna-strandness RF' # stranded libraries
#CE2_PARAMS = ['--no-fix'] #suggested not to set '--no-fix' in real datasets

#MIN_READS = 2 #default
#CIRCRNA_METHODS = 'dcc,ciri,findcirc,testrealign,'\
#                  'circexplorer2_star,circexplorer2_segemehl,'\
#                  'circexplorer2_bwa,circexplorer2_tophat,circrna_finder'
#
#TOGGLE_TRANSCRIPTOME_RECONSTRUCTION = 'True'

## aligners' custom parameters
## parameters from CIRI
#BWA_PARAMS = ['-T', '19', '-c', '1']
## parameters from CIRCexplorer2
#SEGEMEHL_PARAMS = ['-M','1'] #'-D', '0', '-Z', '20'
#TOPHAT_PARAMS = ['--max-multihits', '1']#'--zpacker','pigz' 
## parameters used in DCC manual example
#STAR_PARAMS = ['--outFilterMultimapNmax', '1', 
#               '--outSJfilterOverhangMin', '15', '15', '15', '15',
#               '--alignSJoverhangMin', '15',
#               '--alignSJDBoverhangMin', '15',
#               '--seedSearchStartLmax', '30',
#               '--outFilterScoreMin', '1',
#               '--outFilterMatchNmin', '1',
#               '--outFilterMismatchNmax', '2',
#               '--chimSegmentMin', '15',
#               '--chimScoreMin', '15',
#               '--chimScoreSeparation', '10',
#               '--chimJunctionOverhangMin', '15'] 

#LIN_COUNTER = 'ccp' #'dcc'
#
#DCC_EXTRA_PARAMS = ['-fg', '-M', '-Nr', 1, 1, '-F', '-ss']
#TESTREALIGN_PARAMS = ['-q', 'median_1'] ## suggested 'median_40'
#FINDCIRC_EXTRA_PARAMS = ['--best-qual', '0'] #suggested '40' 
#SAM_SORT_MM = '1G'
##BYPASS = 'linear'
#CCP_COUNTS = 'True'
#
#CIRC_MAPPING = "{'SE':['STAR','TOPHAT','BOWTIE2'],'PE':['BWA','SEGEMEHL']}"
