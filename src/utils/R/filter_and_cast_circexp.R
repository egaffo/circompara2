#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
  make_option(c("-i", "--input"), action = "store", type = "character",
              help = "CircRNA expression table file in 'long' format."),
  make_option(c("-r", "--min_reads"), action = "store", type = "integer",
              default = 2,
              help = "The minimum detection read threshold for a circRNA  (in at least one sample)"),
  make_option(c("-m", "--min_methods"), action = "store", type = "integer",
              default = 2,
              help = "Keep circRNAs commonly detected by >= m circRNA detection methods (in at least one sample)"),
  make_option(c("-o", "--output"), action = "store", type = "character",
              default = "./",
              help = "Output file name")
)

parser <- OptionParser(usage = paste0("%prog -i bks.counts.union.csv -r 2 -m 2 ",
                                      "-o circexp_count_matrix.csv"),
                       option_list = option_list)
arguments <- parse_args(parser, positional_arguments = F)

input_file <- arguments$input
min_reads <- arguments$min_reads
min_methods <- arguments$min_methods
output_file <- arguments$output

tab <- fread(input_file, 
             header = T, 
             showProgress = F)[read.count >= min_reads & 
                                 n_methods >= min_methods]

if("strand" %in% colnames(tab)){
  tab[, circ_id := paste0(chr, ":", start, "-", end, ":", strand)]
}else{
  tab[, circ_id := paste0(chr, ":", start, "-", end)]
}

## TODO: print out some statistics of the filtering

fwrite(x = dcast(tab, circ_id ~ sample_id, value.var = "read.count", fill = 0),
       file = output_file, 
       sep = "\t", 
       row.names = F, 
       col.names = T)
