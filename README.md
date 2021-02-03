---
Title:      CirComPara2  
Subtitle:   CircRNA detection from RNA-seq data using multiple methods  
Project:    CirComPara2  
Author:     Enrico Gaffo  
Affiliation:    Compgen - University of Padova  
Web:        http://compgen.bio.unipd.it  
Date:       January 20, 2021  
---

# CirComPara2

CirComPara2 is a computational pipeline to detect, quantify, and correlate expression of linear and circular RNAs from RNA-seq data that combines multiple circRNA-detection methods.

<!--TODO: more exhaustive description -->

## Quick install

Execute the following commands to download and install (locally) in your system the scripts and tools required to run circompara2. 
If something goes wrong with the installation process try to manually install the software as described below.

Download and extract [the latest release of CirComPara][circompara_pack_link], or clone the GIT repository, enter circompara2 directory and run the automatic installer script:  

```bash
git clone http://github.com/egaffo/circompara2
cd circompara2
./install_circompara
```

### Test your installation

```bash
cd test_circompara/analysis
../../circompara2
```

If you plan to use single-end reads, test with:  

```bash
cd test_circompara/analysis_se
../../circompara2
```

### Add  circompara2 to your environment

Once completed the installation, if you do not want to type the whole path to the circompara2 executable each time, you can update your `PATH` environment variable. From the terminal type the following command (replace the `/path/to/circompara2/install/dir` string with circompara2's actual path)   

```bash
export PATH=/path/to/circompara2/install/dir:$PATH
```

Another way is to link circompara2's main script in your local `bin` directory  

```bash
cd /home/user/bin
ln -s /path/to/circompara2/install/dir/circompara2
```

## circompara2 Docker image

A [Docker image of CirComPara](http://hub.docker.com/r/egaffo/circompara-docker/) is available from DockerHub.

To pull the image:

```bash
docker pull egaffo/circompara2-docker
```
    
You'll find the instructions on how to use the docker image at https://hub.docker.com/r/egaffo/circompara2-docker.

# How to use

## Set your analysis project

This section shows how to set your project directory and run the analysis.
To run an analysis usually you want to specify your data (the sequenced reads in FASTQ format) and a reference genome in FASTA format.

### Compose META file

You have to specify read files and sample names in a metadata table file. The file format is a comma separated text file with the following header line:  

    file,sample
  
Then, each row corresponds to a read file. If you have paired-end sequenced samples write one line per file with the same sample name. 

An example of the metadata table:

file|sample
----|------
/home/user/reads_S1_1.fq|S1
/home/user/reads_S1_2.fq|S1
/home/user/reads_S2_1.fq|S2
/home/user/reads_S2_1.fq|S2


and metadata file content:

    file,sample
    /home/user/reads_S1_1.fq,S1
    /home/user/reads_S1_2.fq,S1
    /home/user/reads_S2_1.fq,S2
    /home/user/reads_S2_1.fq,S2

In the meta file you can also specify the adapter sequences to preprocess the reads, just add an `adapter` column with the adpter file.

file|sample|adapter
----|------|-------
/home/user/reads_S1_1.fq|S1|/home/user/circompara/adapter.fa
/home/user/reads_S1_2.fq|S1|/home/user/circompara/adapter.fa


### Specify the reference genome file

A required parameter is the reference genome. You can either pass the reference genome from the command line

```bash
./circompara2 "GENOME_FASTA='/home/user/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa'"
```

or by setting the `GENOME_FASTA` parameter in the `vars.py` file; e.g.:

```bash
GENOME_FASTA = '/home/user/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
```


### Specify options in vars.py

Although parameters can be set from command line (sorrounded by quotes), you can set them in the `vars.py` file, which must be placed into the directory where circompara2 is called.  
Below there is the full list of the parameters:

```
META: The metadata table file where you specify the project samples, etc.
    default: meta.csv

ANNOTATION: Gene annotation file (like Ensembl GTF/GFF)
    default: 

GENOME_FASTA: The FASTA file with the reference genome
    default: 

CIRCRNA_METHODS: Comma separated list of circRNA detection methods to use. Repeated values will be collapsed into unique values. Currently supported: ciri, dcc, circrna_finder, find_circ, circexplorer2_star, circexplorer2_bwa, circexplorer2_tophat, circexplorer2_segemehl, testrealign (a.k.a. Segemehl). Set an empty string to use all methods available (including deprecated methods). 
    default: ciri,find_circ,circexplorer2_star,circexplorer2_bwa,circexplorer2_segemehl,circexplorer2_tophat,dcc

CPUS: Set number of CPUs
    default: 1

GENEPRED: The genome annotation in GenePred format
    default: 

GENOME_INDEX: The index of the reference genome for HISAT2
    default: 

SEGEMEHL_INDEX: The .idx index for segemehl
    default: 

BWA_INDEX: The index of the reference genome for BWA
    default: 

BOWTIE2_INDEX: The index of the reference genome for BOWTIE2
    default: 

STAR_INDEX: The directory path where to find Star genome index
    default: 

BOWTIE_INDEX: The index of the reference genome for BOWTIE when using CIRCexplorer2_tophat
    default: 

HISAT2_EXTRA_PARAMS: Extra parameters to add to the HISAT2 aligner fixed parameters '--dta --dta-cufflinks --rg-id <SAMPLE> --no-discordant --no-mixed --no-overlap'. For instance, '--rna-strandness FR' if stranded reads are used.
    default: --seed 123

BWA_PARAMS: Extra parameters for BWA
    default: -T 19

SEGEMEHL_PARAMS: SEGEMEHL extra parameters
    default: -D 0

TOPHAT_PARAMS: Extra parameters to pass to TopHat
    default: 

STAR_PARAMS: Extra parameters to pass to STAR
    default: --runRNGseed 123 --outSJfilterOverhangMin 15 15 15 15 --alignSJoverhangMin 15 --alignSJDBoverhangMin 15 --seedSearchStartLmax 30 --outFilterScoreMin 1 --outFilterMatchNmin 1 --outFilterMismatchNmax 2 --chimSegmentMin 15 --chimScoreMin 15 --chimScoreSeparation 10 --chimJunctionOverhangMin 15

BOWTIE2_PARAMS: Extra parameters to pass to Bowtie2 in addition to -p $CPUS --reorder --score-min=C,-15,0 -q
    default: --seed 123

STRINGTIE_PARAMS: Stringtie extra parameters. F.i. '--rf' assumes a stranded library fr-firststrand, to be used if dUTPs stranded library were sequenced  
    default:  

CIRI_EXTRA_PARAMS: CIRI additional parameters
    default: 

DCC_EXTRA_PARAMS: DCC additional parameters
    default: -fg -M -F -Nr 1 1 -N

CE2_PARAMS: Parameters to pass to CIRCexplorer2 annotate
    default:

TESTREALIGN_PARAMS: Segemehl/testrealign filtering parameters-q indicates the minimum median quality of backsplices ends (like the Haarz parameter)
    default: -q median_1

FINDCIRC_EXTRA_PARAMS: Parameters for find_circ.py. Additional parameters: --best-qual INT is used to filter find_circ results according to best_qual_left and best_qual_right fields >= INT. Default: INT = 40. --filter-tags TAG is used to filter lines of find_circ.py output (sites.bed). Repeat it if multiple consecutive filter tags has to be applied.
    default: --best-qual 40 --filter-tags UNAMBIGUOUS_BP --filter-tags ANCHOR_UNIQUE

CFINDER_EXTRA_PARAMS: Parameters for CircRNA_finder 
    default:

PREPROCESSOR: The read preprocessing tool to use. Currently, only "trimmomatic" is supported.Leave empty for no read preprocessing.
    default: 

PREPROCESSOR_PARAMS: Read preprocessor extra parameters. F.i. if Trimmomatic, an empty string defaults to MAXINFO:40:0.5 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:50 AVGQUAL:30 
    default: 

LINEAR_EXPRESSION_METHODS: The method to be used for the linear expression estimates/transcriptome reconstruction. To run more methods use a comma separated list. However, only the first method in the list will be used in downstream processing. Currently supported methods: stringtie,cufflinks,htseq.  
    default: stringtie  

TOGGLE_TRANSCRIPTOME_RECONSTRUCTION: Set True to enable transcriptome reconstruction. Default only quantifies genes and transcripts from the given annotation GTF file
    default: False

READSTAT_METHODS: Comma separated list of methods to use for read statistics. Currently supported: fastqc
    default: fastqc

MIN_METHODS: Number of methods that commmonly detect a circRNA to define the circRNA as reliable. If this value exceeds the number of methods specified, it will be set to the number of methods.
    default: 2

MIN_READS: Number of reads to consider a circRNA as expressed
    default: 2

BYPASS: Skip analysis of linear/circular transcripts. This will also skip the analysis of linear-to-circular expression correlation. The circular analysis includes the pre-filtering of linearly mapping reads. If you want to analyze reads already filtered for linear mappings you should set "linear,linmap". Choose among linear and or linmap, circular. NB: you still have to set the --rna-strandness parameter into the HISAT_EXTRA_PARAMS if you have stranded alignments/reads.
    default:

LINMAPS: You can specify here the path to pre-computed files of linearly aligned reads. This will skip read pre-processing and linear alignments (use jointly to BYPASS linmap to get also circular-to-linear analysis). Mind that the alignments must be in BAM format and the .bai mapping file index file must be in the same directory. NB: you still have to set the --rna-strandness parameter into the HISAT_EXTRA_PARAMS if you have stranded alignments/reads. You need to set a Python dict-like string parameter with sampleName and the corresponding BAM file. E.g: {"SAMPLE1": "sample1/hisat2.bam", "SAMPLE2": "sample2/hisat2.bam"}
    default: None

CIRC_MAPPING: By default (SE), linearly unmapped reads arealigned as single-end reads to search for circRNA backsplices. Set PE to align as paired-end reads by each circRNA method aligner. You can also specify each aligner's mode, or just which aligner has to use the PE mode,  with the syntax for Python dictionaries {'SE':['ALN1','ALN2'],'PE':['ALN3','ALN4','ALNn']} or simply {'PE':['ALN1','ALN2']} if you want just ALN1 and ALN2 tu align as PE. Supported aligners are BWA,SEGEMEHL,STAR and TOPHAT. BOWTIE2 is also supported but it is run only in single-end mode as it serves only Findcirc.
    default: {'SE':['STAR','TOPHAT','BOWTIE2'],'PE':['BWA','SEGEMEHL']}

LIN_COUNTER: The method to estimate circRNA-host gene linear expression. Available are using the DCC [dcc], or the CirComPara [ccp] method.
    default: ccp

FIX_READ_HEADER: Trim FASTQ headers to the read ids. Recommended when processing SRA datasets.
    default: True

UNSTRANDED_CIRCS: Force unstranded circRNAs even if stranded library was used
    default: False

SAM_SORT_MM: Value for samtools sort -m option
    default: 768M

QRE_FIND: (Experimental) Set True to toggle analysis of QKI response elements sequences
    default: False

CCP_COUNTS: Set the strategy to estimate circRNA expression.
    default: UN
```  

## Run the analysis

To trigger the analyses you simply have to call the `./circompara2` script in the analysis directory. Remember that if you used the `vars.py` option file, this has to be in the analysis directory.  


```bash
cd /home/user/circrna_analysis
/home/user/circompara2/circompara2
```

### Additional options from the Scons engine:

* *Basic execution*: run the analysis as a linear pipeline, i.e. no parallel task execution, and stop on errors

```bash
/path/to/circompara2/dir/circompara2
```

* *Show parameters*: to show the parameters set before actually run the analysis, use `-h`:
```bash
/path/to/circompara2/dir/circompara2 -h
```

* *Dryrun*: to see which commands will be executed without actually execute them, use the `-n` option. NB: many commands will be listed, so you should redirect to a file or pipe to a reader like `less`
```bash
/path/to/circompara2/dir/circompara2 -n | less -SR
```

* *Multitasks*: the `-j` option specifies how many **tasks** can be run in parallel. N.B: "j x CPUS <= available cores", i.e: the j option value times the CPUS parameter value should not be greater than the number of CPU cores available, unless you want to overload your machine. 
```bash
/path/to/circompara2/dir/circompara2 -j4
```

* *Ignore errors*: keep executing the tasks even when some of them fails. Caveat: this can break downstream analyses
```bash
/path/to/circompara2/dir/circompara2 -i
```

* *Combine options*: to set multiple options you must sorround them with quotes
```bash
/path/to/circompara2/dir/circompara2 "-j4 -i"
```

## Output files

* Statistics on the read filtering steps and alignments can be found into `read_statistics` directory. A report is saved in the `processing_and_mapped_read_counts.csv` file.  
* Results regarding circRNAs (expression matrices, etc.) will be saved into the `circular_expression/circrna_analyze` directory: the reliable circRNA expression matrix is in the `reliable_circexp.csv` file; gene annotation associated to the circRNAs is stored in the files under the `circular_expression/circrna_collection/circrna_gene_annotation` directory.  
* CircRNA parent gene linear expression is saved in the `circular_expression/circrna_analyze/ccp_bks_linexp.csv` file.  
* Gene expression values for each gene and sample will be saved in the `linear_expression/linear_quantexp/geneexp/` directory: `gene_expression_FPKM_table.csv` file reports FPKMs and `gene_expression_analysis.html` file reports summary analysis.  
* Linear transcript sequences are saved as a multi-FASTA file into the `linear_expression/transcript_sequences` directory.

# Advanced features

## Make genome indexes for multiple instances of circompara2: the `make_indexes` utility

Building the genome indexes for each mapper can take lot of computing time. However, the same indexes can be used in different circompara2 runs, saving time and disk space. In circompara2's package the `src/utils/bash/make_indexes` script can be used to automatically build the genome index (and gene annotation formats) for each of the supported read aligner, and save them into a directory. In addition, it gives the parameter values to be set to use the index files to be shared.  
Example commands using the test data follows:  
```bash
cd test_circompara
mkdir genome_indexes
cd genome_indexes
../../src/utils/bash/make_indexes "-j2 GENOME=../annotation/CFLAR_HIPK3.fa ANNOTATION=../annotation/CFLAR_HIPK3.gtf" 
```

The above commands will eventually generate a `annotation_vars.py` file that can be appended to the `vars.py` file of your project so that circompara2 will skip the building of genome indexes. Note that `make_indexes` can use the same options provided by Scons showed above: `-j 2` option will allow the script to build two indexes in parallel.  

```bash
cd test_circompara
## clear circompara2 files in the test directory
cd analysis
../../circompara2 -c
cd ..
## overwrite the vars.py file omitting the genome and annotation parameters
grep -v "GENOME\|ANNOTATION" vars.py > analysis/vars.py
## append the parameters for the genome, the annotation and the genome indexes
## generated by the make_indexes utility
cat genome_indexes/annotation_vars.py >> analysis/vars.py
## run the test analysis
cd analysis
../../circompara2
```

## Stranded libraries

Some tools in circompara2 require special parameters to handle properly stranded reads. circompara2 allows to specify such parameters
Example: include the following parameters if you used the Illumina TruSeq Stranded Total RNA Library Prep Kit with Ribo-Zero Human/Mouse/Rat

    HISAT2_EXTRA_PARAMS = "--rna-strandness FR "

<!-- Experimental
## Fusion genes and fusion circular RNAs
If you want to analyze fusion genes and enable detection of fusion circular RNAs (f-circRNAs) you have to include a `translocation` column to the metadata file. This field specifies the genomic coordinates of the gene pair involved in the fusion. You do not have to specify fusion breakpoints as the transcript structure is inferred by the transcript reconstruction algorithm (currently [Cufflinks][]).   
Gene coordinates must be defined as follow:

    chr:start-end:strand&chr:start-end:strand

with strand being either + or -; note that the two genes' coordinates are separated by a `&` character.
More gene pairs can be specified in one row, you just have to separate the pairs by a `#` character
Below, a metadata file example where the MLL-AF4 fusion gene is specified, as well as its "complementary" translocation:

    file,sample,condition,translocation
    /home/user/reads_S1_1.fq,S1,WT,11:34488-128832:+&4:21033000-21243056:+#4:21033000-21243056:+&11:34488-128832:+
    /home/user/reads_S1_2.fq,S1,WT,
    /home/user/reads_S2_1.fq,S2,MU,
    /home/user/reads_S2_1.fq,S2,MU,

As you can note from the exampple, you do not have to specify the gene pair in each line. Moreover, different samples can be set with different gene pairs.
Fusion gene analysis will then been performed on all samples for all the fusion gene set. 




## Advanced parameters: the vars.py file

Type

    circompara_circompara2 -h
    
to show an help of all parameters

## Advanced features output

**TODO** Output for f-circRNAs ...

### Fusion genes
When enabled, the fusion gene analysis will generate "synthetic" chromosomes ... **TODO**    

 -->

TODO: update software versions!

Software|Website|Version
--------|-------|-------:
Scons|http://www.scons.org|3.1.2
Trimmomatic|http://www.usadellab.org/cms/?page=trimmomatic|0.39
FASTQC|http://www.bioinformatics.babraham.ac.uk/projects/fastqc/|0.11.9
HISAT2|http://ccb.jhu.edu/software/hisat2/index.shtml|2.1.0
STAR|http://github.com/alexdobin/STAR|2.6.1e
BWA|http://bio-bwa.sourceforge.net/|0.7.15-r1140
Bowtie2|http://bowtie-bio.sourceforge.net/bowtie2/index.shtml|2.4.1
Bowtie|http://bowtie-bio.sourceforge.net/index.shtml|1.1.2
TopHat|http://ccb.jhu.edu/software/tophat/index.shtml|2.1.0
Segemehl|http://www.bioinf.uni-leipzig.de/Software/segemehl/|0.3.4
CIRI|http://ciri.sourceforge.io/|2.0.6
CIRCexplorer2|http://github.com/YangLab/CIRCexplorer|2.3.8
find_circ|http://github.com/marvin-jens/find_circ|1.2
BEDtools|http://bedtools.readthedocs.io|2.29.2
Samtools|http://www.htslib.org/|1.10

### Errors with R packages

If you get error messages from R packages of your already installed circompara2, maybe some update occurred in your R system. Try to re-install all circompara2 R package dependencies by using the `reinstall_R_pkgs` command.  

<!-- ## Details on circompara2 architecture and implementation

circompara2 is part of the [circompara][circompara_link] project, which is a collection of scripts consisting of modules that can be assembled to compose computational pipelines for RNA-seq data analysis.   
Each circompara module consists of a script written within the [Scons][scons_link] building tool (say a 'Sconscript') that can implement its function also by wrapping other software. Each module can be used either standalone or included in a pipeline of commands, as it was done for circompara2.

The core engine is the Scons build tool, which manage the various steps of the analysis.

**TODO**
-->

# How to cite
If you used CirComPara2 for your analysis, please add the following citation to your references:  

citation [circompara_article]


[scons_link]: http://scons.org/
[circompara_git_link]: http://github.com/egaffo/circompara2 "circompara Git repository"
[circompara_pack_link]: http://github.com/egaffo/circompara2/releases/latest "circompara package"
[test_data_link]: http://github.com/egaffo/circompara2
[circompara_link]: http://github.com/egaffo/circompara2
[circompara_article]: url

