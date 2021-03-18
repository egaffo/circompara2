#!/usr/bin/parallel --shebang-wrap bash

## usage:
# find dir_with_the_BAM_files/ | grep "_hisat2.bam$" > bam_list.txt
# get_unmapped_reads_from_bam.sh < bam_list.txt
# pigz *.fq

SAMPLE=`basename $1 | sed "s@\.bam@@"`
samtools fastq -f 12 -F 3328 -n -s $SAMPLE\_singleton.fq -1 $SAMPLE\_unmapped_1.fq -2 $SAMPLE\_unmapped_2.fq $1

