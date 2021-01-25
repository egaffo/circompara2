#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

option_list <- list(
  make_option(c("-o", "--results.dir"), action = "store", type = "character",
              default = "./",
              help = "The output directory path"),
  make_option(c("-t", "--transcripts.gtf.files"), action = "store", type = "character",
              # default = "",
              help = "The transcripts.gtf files as output by Stringtie"),
  make_option(c("-r", "--gene_raw_counts_list"), action = "store", type = "character",
              # default = "",
              help = ""),
  make_option(c("-g", "--gene.xpr.file"), action = "store", type = "character",
              # default = "",
              help = ""),
  make_option(c("-x", "--trx_raw_counts_list"), action = "store", type = "character",
              # default = "",
              help = "")
)

parser <- OptionParser(usage = "%prog ",
                       option_list = option_list,
                       description = "")

arguments <- parse_args(parser, positional_arguments = F)

results.dir             <- arguments$results.dir
gene.xpr.file           <- readLines(arguments$gene.xpr.file)
transcripts.gtf.files   <- readLines(arguments$transcripts.gtf.files)
gene_raw_counts_list    <- readLines(arguments$gene_raw_counts_list)
trx_raw_counts_list     <- readLines(arguments$trx_raw_counts_list)

## assume input is from StringTie

gene.xpr.list <- lapply(gene.xpr.file, fread)
names(gene.xpr.list) <- sub("_gene_abund.tab", "", 
                            sapply(gene.xpr.file, basename, 
                                   USE.NAMES = F))
gene.xpr <- rbindlist(gene.xpr.list, use.names = T, idcol = "sample")

## TPM
gene.xpr.tpm <- dcast(data = gene.xpr,
                      formula = `Gene ID` + `Gene Name` + Reference + Strand + Start + End ~ sample,
                      value.var = "TPM", 
                      fill = 0)

fwrite(x = gene.xpr.tpm, 
       file = file.path(results.dir, "gene_expression_TPM_table.csv"), 
       sep = "\t",
       col.names = T,
       row.names = F)

## N reads (a.k.a. coverage)
gene.xpr.nreads <- dcast(data = gene.xpr, 
                         formula = `Gene ID` + `Gene Name` + Reference + Strand + Start + End ~ sample, 
                         value.var = "Coverage", 
                         fill = 0)
fwrite(x = gene.xpr.nreads, 
       file = file.path(results.dir, "gene_expression_Nreads_table.csv"), 
       sep = "\t", 
       col.names = T,
       row.names = F)

genes.read_group_tracking <- gene.xpr[, .(sample_id = sample, tracking_id = `Gene ID`, FPKM)]

## merge gene_raw_counts_list
gene_raw_counts <- lapply(gene_raw_counts_list, fread)
names(gene_raw_counts) <- sub("_gene_expression_rawcounts.csv", "",
                              sapply(gene_raw_counts_list, basename,
                                     USE.NAMES = F))
gene_raw_counts <- rbindlist(gene_raw_counts, use.names = T, idcol = "sample_id")

fwrite(x = dcast(gene_raw_counts,
                 formula = gene_id ~ sample_id,
                 value.var = "raw.reads", 
                 fill = 0),
       file = file.path(results.dir, "gene_expression_rawcounts_table.csv"), 
       sep = "\t",
       col.names = T,
       row.names = F)

## merge trx_raw_counts_list
trx_raw_counts <- lapply(trx_raw_counts_list, fread)
names(trx_raw_counts) <- sub("_transcript_expression_rawcounts.csv", "",
                             sapply(trx_raw_counts_list, basename,
                                    USE.NAMES = F))
trx_raw_counts <- rbindlist(trx_raw_counts, use.names = T, idcol = "sample_id")

fwrite(x = dcast(trx_raw_counts[raw.reads > 0, ],
                 formula = gene_id + transcript_id ~ sample_id,
                 value.var = "raw.reads", 
                 fill = 0),
       file = file.path(results.dir, "transcript_expression_rawcounts_table.csv"), 
       sep = "\t", 
       col.names = T,
       row.names = F)

# expressed_genes <- genes.read_group_tracking[, FPKM := round(FPKM, digits = 8)][FPKM > 0]

gene.xpr.fpkm <- 
  dcast(data = gene.xpr,
        formula = `Gene ID` + `Gene Name` + Reference + Strand + Start + End ~ sample,
        value.var = "FPKM", 
        fill = 0)

colnames(gene.xpr.fpkm)[colnames(gene.xpr.fpkm) == "Gene ID"] <- "gene"

gene.xpr.fpkm <- 
  gene.xpr.fpkm[, `:=`("Gene Name" = NULL, 
                       "Reference"= NULL, 
                       "Strand" = NULL, 
                       "Start" = NULL, 
                       "End" = NULL)][]

expressed_genes_table.file <- file.path(results.dir, "gene_expression_FPKM_table.csv")
fwrite(x = gene.xpr.fpkm, 
       file = expressed_genes_table.file, 
       row.names = F, 
       col.names = T, 
       sep = "\t")
