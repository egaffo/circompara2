#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

options(scipen=999)

option_list <- list(
    make_option(c("-i", "--input"), action = "store", type = "character",
                default = "sample.sngl.bed",
                help = "Segemehl result table file"),
    make_option(c("-q", "--minqual"), action = "store", type = "character",
                default = "median_10",
                help = paste("Filter function and minimum mapping quality of backsplice-reads. ",
                       "Function supported: ",
                       "- median: for each backsplice compute the read mapping median quality and discard backsplices with median quality <= N",
                       "- mean: as for median, but compute the mean instead",
                       "- any: first, discard any alignment with mapping quality <= N. Then, report the median as backsplice quality",
                       sep = "\n")),
    make_option(c("-o", "--count_output"), action = "store", type = "character",
                default = "splicesites.bed",
                help = "The filtered circrnas in BED format. The 7th column is just a copy of the score field and reports the read count."),
    make_option(c("-r", "--reads_output"), action = "store", type = "character",
                default = "sample.circular.reads.bed.gz",
                help = "The output file"),
    make_option(c("-m", "--keep_mates"), action = "store_true", #type = "logical",
                default = F,
                help = "Whether to include /1 or /2 in read names to distinguish readmates"),
    make_option(c("-t", "--trns_file"), action = "store", type = "character",
                default = "unmapped_1.fastq.trns.txt",
                help = "The *trns.txt file output by Segemehl"),
    make_option(c("-l", "--legacy_output"), action = "store", type = "character",
                default = "old_splicesites.bed",
                help = "The filename for Segemehl splicesites.bed in Segemehl < v0.3 format")
)

parser <- OptionParser(usage = paste0("%prog -i sample.sngl.bed -q median_10 -o splicesites.bed",
                                      " -r sample.circular.reads.bed.gz -t unmapped_1.fastq.trns.txt"),
                       option_list = option_list,
                       description = "Filter segemehl results by alignments' quality and compute backsplices' read count")
arguments <- parse_args(parser, positional_arguments = F)
input <- arguments$input
count_output <- arguments$count_output
reads_output <- arguments$reads_output
min.qual <- arguments$minqual
trns.file <- arguments$trns_file
old_output <- arguments$legacy_output

split.qual.par <- strsplit(min.qual, "_", fixed = T)[[1]]
if(length(split.qual.par) > 1){
    qual.func <- split.qual.par[1]
    if(!qual.func %in% c("median", "mean", "any")){
        message("Error in evaluating quality function in '", min.qual, "': setting default 'median'")
        qual.func <- "median"
    }
    minqual <- suppressWarnings(as.numeric(split.qual.par[2]))
    if(is.na(minqual)){
        message("Error in evaluating quality value in '", min.qual, "': setting default '10'")
        minqual <- 10
    }
}else{
    message("Error in evaluating quality function '", min.qual, "': setting default 'median_10'")
    qual.func <- "median"
    minqual <- 10
}

# Start and end position indicate the genomic range of the predicted intron.
# The name has the format (read-group;type;read-name;mate-status), the bed
# score is the alignment score of the respective alignment. The type is either
# 'R' (in case of a regular, collinear split), 'C' (circular split) or 'B' (backsplice)

sege_circ <- fread(cmd = paste0('grep ";B\\|C;" ', input),
                   header = F, skip = 1)

# Add backspliced reads from the transplice file (backsplices > 20000 bp length)
trns.fields <- c("chr","pos","strand","start.in.read","align.length","algin.edist","score")

trns <- fread(trns.file, header = F, sep = "\t")
trns[, paste0(trns.fields,"_L") := tstrsplit(V1, ",", type.convert = T)]
trns[, paste0(trns.fields,"_R") := tstrsplit(V2, ",", type.convert = T)]

trns.samechr.reads <-
    trns[chr_L == chr_R &
             strand_L == strand_R &
             align.length_L >= 20 & align.length_R >= 20]

trns.bks.reads <-
    rbindlist(list(trns.samechr.reads[strand_L == "+" & pos_L > pos_R, ## select backsplices on forward strand
                                      .(chr = chr_L,
                                        start = as.integer(pos_R) - 1, ## scale by 1 to comply with BED format
                                        end = as.integer(pos_L) + as.integer(align.length_L) - 1,
                                        read.name = V3,
                                        score = ifelse(score_L > score_R,
                                                       as.integer(score_L),
                                                       as.integer(score_R)),
                                        strand = strand_L)],
                   trns.samechr.reads[strand_L == "-" & pos_L < pos_R,  ## select backsplices on reverse strand
                                      .(chr = chr_L,
                                        start = as.integer(pos_L) - 1, ## scale by 1 to comply with BED format
                                        end = as.integer(pos_R) + align.length_R - 1,
                                        read.name = V3,
                                        score = ifelse(score_L > score_R,
                                                       as.integer(score_L),
                                                       as.integer(score_R)),
                                        strand = strand_L)]),
              use.names = T)[, .(V1 = chr, V2 = start, V3 = end,
                                            V5 = score, V6 = strand, read.name)]

splicesites.bed <- data.table()
reads_output.dt <- data.table()

if(qual.func == "any"){
    sege_circ <- sege_circ[V5 >= minqual]
    trns.bks.reads <- trns.bks.reads[V5 >= minqual]
}

if(nrow(sege_circ) > 0){

    if(arguments$keep_mates){
        sege_circ[, c("read.group", "type", "read.name",
                      "mate.status"):=(tstrsplit(V4, ";"))][, read.name :=
                                                                paste0(read.name, "/",
                                                                       mate.status)][, `:=`(V4 = NULL,
                                                                                            read.group = NULL,
                                                                                            type = NULL,
                                                                                            mate.status = NULL)]
    }else{
        sege_circ[, c("read.group", "type", "read.name",
                      "mate.status"):=(tstrsplit(V4, ";"))][, `:=`(V4 = NULL,
                                                                   read.group = NULL,
                                                                   type = NULL,
                                                                   mate.status = NULL)]
    }

    sege_circ <- rbindlist(list(sege_circ,
                                trns.bks.reads),
                           use.names = T)

    ## remove duplicated lines/alignments
    sege_circ <- sege_circ[, .(multi.mapping = .N,
                               map.qual = max(V5)),
                           by = .(chr = V1, left = V2,
                                  right = V3, read.name,
                                  strand = V6)]

    if(qual.func == "any" || qual.func == "median"){
        splicesites.bed <-
            sege_circ[, .(n = .N,
                          map.qual = median(map.qual)),
                      by = .(chr, left, right,
                             strand)][, .(chr, left, right, n, map.qual, strand,
                                          score = n)][order(chr, left,
                                                            right)]
    }

    if(qual.func == "mean"){
        splicesites.bed <-
            sege_circ[, .(n = .N,
                          map.qual = mean(map.qual)),
                      by = .(chr, left, right,
                             strand)][, .(chr, left, right, n, map.qual, strand,
                                          score = n)][order(chr, left,
                                                            right)]
    }

    if(qual.func == "median" || qual.func == "mean"){
        ## filter by quality
        splicesites.bed <- splicesites.bed[map.qual >= minqual]
    }

    ## report selected backsplices' read names
    # reads_output.dt <- sege_circ[, .(chr, left, right, read.name, map.qual, strand)]
    reads_output.dt <- sege_circ[splicesites.bed[, .(chr, left, right, strand)],
                                 on = c("chr", "left", "right",
                                        "strand")][, .(chr, left, right, read.name,
                                                       map.qual, strand)]
}

## write backsplice counts
write.table(x = splicesites.bed,
            file = count_output,
            row.names = F, quote = F, sep = "\t", col.names = T)

## write gzipped file for circular reads
reads_output.gz <- gzfile(reads_output, "w")
write.table(x = reads_output.dt,
            file = reads_output.gz,
            row.names = F, quote = F, sep = "\t", col.names = F)
close(reads_output.gz)

## write old segemehl format
write.table(x = splicesites.bed[, .(V4 = paste("splits", sum(n), sum(n), sum(n), "C", "P", sep=":")),
                                by = .(chr, left, right)][order(chr, left, right),
                                                          .(chr, left = left+1, right, V4, 0, "+")],
            file = old_output,
            row.names = F, quote = F, sep = "\t", col.names = F)

