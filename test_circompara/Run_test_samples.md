To run this example you will need a Linux platform with the Docker engine installed.

## 1. Get test input data

Use the SRA toolkit <a href="https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump">fastq-dump</a> script to download Drosophila RNA-seq data in FASTQ format
```bash
mkdir SRR1197368
cd SRR1197368
fastq-dump --split-files --gzip SRR1197368
```

The above command will generate two compressed FASTQ files in your directory: SRR1197368_1.fastq.gz and SRR1197368_2.fastq.gz.  

## 2. Download genome sequence and gene annotation 
```bash
wget http://ftp.ensembl.org/pub/current_fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz
gunzip Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/current_gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.104.gtf.gz
gunzip Drosophila_melanogaster.BDGP6.32.104.gtf.gz
```

## 3. Write the vars.py file to look like
```python
META            = 'meta.csv'
GENOME_FASTA    = '/data/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa'
ANNOTATION      = '/data/Drosophila_melanogaster.BDGP6.32.104.gtf' 
CPUS            = '4'
```

Change the number of CPUS according to your workstation/server specs.

## 4. Write the meta.csv file, which will look like 
```
sample,file
SRR1197368,/data/SRR1197368_1.fastq.gz
SRR1197368,/data/SRR1197368_2.fastq.gz
```

Your directory now contains the following files:

- SRR1197368_1.fastq.gz
- SRR1197368_2.fastq.gz
- Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa
- Drosophila_melanogaster.BDGP6.32.104.gtf
- meta.csv
- vars.py

## 5. Run CirComPara2 through Docker, using 4 parallel tasks. Mind that 16 cores will be used, in total (max 4 cores per each parallel task).  
```bash
docker run -u `id -u` --rm -it -v $(pwd):/data egaffo/circompara2:v0.1.2 '-j4'
```

## 6. Check results

TODO
CircRNA expression

Backsplice linearly spliced read counts

Gene expression and input files for tximport

Summary of read processing statistics

## Running time
