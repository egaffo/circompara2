#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

ncpus <- as.integer(Sys.getenv('CPUS'))
if (is.na(ncpus)) ncpus <- 1
setDTthreads(threads = ncpus)

option_list <- list(
    make_option(c("-i", "--input"), action = "store", type = "character",
                default = "circular.reads.bed.gz.txt",
                help = "File with the list of circular.reads.bed.gz files to merge, each sample and method (default: circular.reads.bed.gz.txt)"),
    make_option(c("-q", "--min_methods"), action = "store", type = "integer",
                default = 2,
                help = "Minimum number of methods (default: 2)"),
    make_option(c("-o", "--output_prefix"), action = "store", type = "character",
                default = "bks.counts.",
                help = "A prefix for output file names (default: bks.counts.intersect.csv, bks.counts.union.csv, bks.counts.union.intersected.csv)"),
    # make_option(c("-f", "--format"), action = "store", type = "character",
    #             default = "BED",
    #             help = "The file format of output files {GTF,BED}"),
    make_option(c("-s", "--stranded"), action = "store_true", default = F,
                help = "Set this option if circRNA strand has to be considered"),
    make_option(c("-c", "--circrnas_gtf"), action = "store", type = "character", default = NA,
                help = "A circrnas.gtf file or a text file listing circrnas.gtf file paths to merge"),
    make_option(c("-t", "--xprtypes"), action = "store", type = "character", default = "UN_IN_IU_MD",
                help = paste("An underscore separated list of the strategy(ies) that will be used to combine and report",
                             "the expression estimates. The options available are:",
                             "UN = combine and count all the unique backsplice junction read fragments (BJRs) from any circRNA detection method;",
                             "IN = count only the BJRs commonly identified by at least 'min_methods' methods",
                             "IU = at least 1 BJRs is required to be commonly identified by at least 'min_methods' methods, then count all other BJRs",
                             "MD = compute the median of the read count reported by each method that detected the circRNA",
                             "ME = (not yet implemented) compute the mean of the read count reported by each method that detected the circRNA"))
)

parser <- OptionParser(usage = "%prog -i circular.reads.bed.gz.txt -q 2 -o bks.counts.",
                       option_list = option_list,
                       description = "Compute backsplices' read counts")

arguments <- parse_args(parser, positional_arguments = F)

input <- arguments$input
output_prefix <- arguments$output_prefix
min_methods <- arguments$min_methods
# file.ext <- tolower(arguments$format)
file.ext <- "csv"
stranded <- arguments$stranded
circrnas.gtf.files <- arguments$circrnas_gtf
strategies <- unlist(strsplit(arguments$xprtypes, split = "_"))

circular.reads.bed.gz.txt <- readLines(input)

## filter out empty files
keep <- sapply(circular.reads.bed.gz.txt, file.exists, simplify = T)
if(any(!keep))print(paste0(circular.reads.bed.gz.txt[!keep], " skipped because not exisitng"))
circular.reads.bed.gz.txt <- circular.reads.bed.gz.txt[keep]

keep <- sapply(circular.reads.bed.gz.txt, function(x){file.size(x)>0}, simplify = T)
if(any(!keep))print(paste0(circular.reads.bed.gz.txt[!keep], " skipped because 0 size"))
circular.reads.bed.gz.txt <- circular.reads.bed.gz.txt[keep]

keep <- sapply(circular.reads.bed.gz.txt,
               function(x){length(grep("^$", x = readLines(x), invert = T))>0},
               simplify = T)
if(any(!keep))print(paste0(circular.reads.bed.gz.txt[!keep], " skipped because empty"))
circular.reads.bed.gz.txt <- circular.reads.bed.gz.txt[keep]


bks.read.method <-
    unique(rbindlist(sapply(circular.reads.bed.gz.txt, fread, header = F,
                            col.names = c("chr", "start", "end", "read_id", "score", "strand"),
                            simplify = F),
                     idcol = "sample_id")[, `:=`(sample_id = sub(".circular.reads.bed.gz", "",
                                                                 basename(sample_id)),
                                                 circ_method = sub(".+circRNAs/([^/]+)/.*", "\\1",
                                                                   sample_id))][])
## fix DCC coordinates
bks.read.method[circ_method == "dcc", start := start -1L]

## fix possible duplicate names of reads (this will count sequenced fragments, not reads)
## given by "\1", "\2" mate name extension, and/or "revcomp_of_" prefix when stranded reads
bks.read.method[, read_id := sub("\\\\[12]$", "", read_id)]
bks.read.method[, read_id := sub("^revcomp_of_", "", read_id)]

bks.read.method <- unique(bks.read.method)

## save collected reads
filename <- paste0(output_prefix, "collected_reads.", file.ext)
write.csv(x = bks.read.method[, .(sample_id, chr, start, end, read_id, strand, circ_method)],
          file = filename,
          row.names = F)

##
if(stranded){
    bks.reads <-
        bks.read.method[, .(n_methods = length(unique(circ_method)),
                            circ_methods = list(unique(circ_method))),
                        by = .(sample_id, chr, start, end, read_id,
                               strand)]
}else{
    bks.reads <-
        bks.read.method[, .(n_methods = length(unique(circ_method)),
                            circ_methods = list(unique(circ_method))),
                        by = .(sample_id, chr, start, end,
                               read_id)]
}

if("IN" %in% strategies){
    ## STRATEGY A)
    ## for each backsplice, count only reads detected by >= min_methods
    ## This should improve reliability of read counts
    ## n_methods_partials: comma-separated list of read.count@n_methods.
    ##                     Example
    ##                     33@5,56@2,8@3,1@2 means that for 98 reads in total:
    ##                     - 33 reads were commonly detected by 5 methods,
    ##                     - 56 reads were commonly detected by 2 methods,
    ##                     - 8 reads were commonly detected by 3 methods,
    ##                     - 1 reads were commonly detected by 2 methods (in a different combination with respect to the other 56 reads)
    ## methods_partials: comma-separated list of read.count@method_names.
    ##                   As for n_methods_partials, giving the method names combinations
    if(stranded){
        # bks.read.counts.intersect <-
        #     bks.reads[n_methods >= min_methods, .N,
        #               by = .(sample_id, chr, start, end, strand, methods,
        #                      n_methods)][, .(read.count = sum(N),
        #                                      n_methods_partials = paste0(paste0(N, "@", n_methods),
        #                                                                  collapse = ","),
        #                                      methods_partials = paste0(paste0(N, "@", methods),
        #                                                                collapse = ",")),
        #                                  by = .(sample_id, chr, start, end, strand)]

        bks.read.counts.intersect <-
            bks.reads[n_methods >= min_methods,
                      .(read.count = .N,
                        circ_methods = list(unique(unlist(circ_methods)))),
                      by = .(sample_id, chr, start, end,
                             strand)][, .(n_methods = length(unlist(circ_methods)),
                                          circ_methods = paste0(sort(unlist(circ_methods)),
                                                                collapse = "|")),
                                      by = .(sample_id, chr, start, end, strand,
                                             read.count)]

    }else{

        # bks.read.counts.intersect <-
        #     bks.reads[n_methods >= min_methods, .N,
        #               by = .(sample_id, chr, start, end, methods,
        #                      n_methods)][, .(read.count = sum(N),
        #                                      n_methods_partials = paste0(paste0(N, "@", n_methods),
        #                                                                  collapse = ","),
        #                                      methods_partials = paste0(paste0(N, "@", methods),
        #                                                                collapse = ",")),
        #                                  by = .(sample_id, chr, start, end)]

        bks.read.counts.intersect <-
            bks.reads[n_methods >= min_methods,
                      .(read.count = .N,
                        circ_methods = list(unique(unlist(circ_methods)))),
                      by = .(sample_id, chr, start,
                             end)][, .(n_methods = length(unlist(circ_methods)),
                                       circ_methods = paste0(sort(unlist(circ_methods)),
                                                             collapse = "|")),
                                   by = .(sample_id, chr, start, end,
                                          read.count)]

    }


    filename <- paste0(output_prefix, "intersect.", file.ext)
    fwrite(x = bks.read.counts.intersect,
           file = filename,
           row.names = F, col.names = T, sep = "\t")
}

if("UN" %in% strategies){
    ## STRATEGY B)
    ## count all reads for each backsplice and keep it if >= min_methods detected
    ## any read for it, i.e: keep the backsplice and its read count if different
    ## methods detected it also if by different read ids
    ## F.i, if min_methods = 2:
    ## - keep bks Y with count "100 reads" in total: 50 reads are detected only
    ## by DCC and 50 only by CIRI, no reads detected by both methods
    ## - discard bks Z with 100 reads detected only by DCC
    ## This should improve sensibility of detection
    if(stranded){
        bks.read.counts.union <-
            bks.reads[, .(read.count = .N,
                          circ_methods = list(unique(unlist(circ_methods)))),
                      by = .(sample_id, chr, start, end,
                             strand)][, .(n_methods = length(unlist(circ_methods)),
                                          circ_methods = paste0(sort(unlist(circ_methods)),
                                                                collapse = "|")),
                                      by = .(sample_id, chr, start, end,
                                             strand, read.count)]
        #[n_methods >= min_methods] ## do not filter already
    }else{
        bks.read.counts.union <-
            bks.reads[, .(read.count = .N,
                          circ_methods = list(unique(unlist(circ_methods)))),
                      by = .(sample_id, chr, start,
                             end)][, .(n_methods = length(unlist(circ_methods)),
                                       circ_methods = paste0(sort(unlist(circ_methods)),
                                                             collapse = "|")),
                                   by = .(sample_id, chr, start, end,
                                          read.count)]
        #[n_methods >= min_methods] ## do not filter already
    }

    filename <- paste0(output_prefix, "union.", file.ext)
    fwrite(x = bks.read.counts.union,
           file = filename,
           row.names = F, col.names = T, sep = "\t")
}

if("IU" %in% strategies){
    ## STRATEGY C)
    ## for each backsplice, if one read was detected by >= min_methods then count
    ## also reads detected by < min_methods
    ## F.i, if min_methods = 2:
    ## - keep bks Y with count "100 reads" in total: 50 reads are detected only
    ## by DCC and 49 only by CIRI, plus 1 read which is detected by both methods;
    ## - discard bks Z with count "100 reads" in total: 50 reads are detected only
    ## by DCC and 50 only by CIRI, no reads detected by both methods
    if(stranded){
        bks <- unique(bks.reads[n_methods >= min_methods, .(sample_id, chr, start, end, strand)])
        # bks.read.counts.union.intersected <-
        #     bks.reads[bks,
        #               on = c("sample_id", "chr", "start", "end",
        #                      "strand")][, .(read.count = .N),
        #                                 by = .(sample_id, chr, start, end, strand)]

        bks.read.counts.union.intersected <-
            bks.reads[bks,
                      .(read.count = .N,
                        circ_methods = list(unique(unlist(circ_methods)))),
                      by = .(sample_id, chr, start, end,
                             strand),
                      on = c("sample_id", "chr", "start",
                             "end"),
                      nomatch = NULL][, .(n_methods = length(unlist(circ_methods)),
                                          circ_methods = paste0(sort(unlist(circ_methods)),
                                                                collapse = "|")),
                                      by = .(sample_id, chr, start, end, strand,
                                             read.count)]
    }else{
        bks <- unique(bks.reads[n_methods >= min_methods, .(sample_id, chr, start, end)])
        # bks.read.counts.union.intersected <-
        #     bks.reads[bks, on = c("sample_id", "chr", "start",
        #                           "end")][, .(read.count = .N),
        #                                   by = .(sample_id, chr, start, end)]

        bks.read.counts.union.intersected <-
            bks.reads[bks,
                      .(read.count = .N,
                        circ_methods = list(unique(unlist(circ_methods)))),
                      by = .(sample_id, chr, start,
                             end),
                      on = c("sample_id", "chr", "start",
                             "end"),
                      nomatch = NULL][, .(n_methods = length(unlist(circ_methods)),
                                          circ_methods = paste0(sort(unlist(circ_methods)),
                                                                collapse = "|")),
                                      by = .(sample_id, chr, start, end,
                                             read.count)]

    }

    filename <- paste0(output_prefix, "union.intersected.", file.ext)
    fwrite(x = bks.read.counts.union.intersected,
           file = filename,
           row.names = F, col.names = T, sep = "\t")

}


# ## save in a matrix-like format
# bks.read.counts <- bks.read.counts.union
# bks.read.counts.wide <-
#     dcast(data = bks.read.counts,
#           formula = chr + start + end + strand ~ sample_id,
#           value.var = "read.count",
#           fill = 0)
#
# filename <- paste0(output_prefix, "")
# write.csv(x = ,
#           file = filename,
#           row.names = F)

if("MD" %in% strategies){
    if(!is.na(circrnas.gtf.files)){
        if(ncol(fread(circrnas.gtf.files, showProgress = F, nrows = 1)) == 1){
            ## case: text file listing circrnas.gtf files
            circrnas.gtf.files <- readLines(circrnas.gtf.files)
        }

        ## read circRNA results: filter low expressed (less than min_reads reads) circRNAs
        colClasses <- c("factor", "factor", "character", "integer",
                        "integer", "integer", "factor", "character", "character")
        circrnas.gtf <- rbindlist(lapply(circrnas.gtf.files, fread, data.table = T,
                                         colClasses = colClasses, showProgress = F),
                                  use.names = T)

        if(stranded){
            circrnas.gtf[, `:=` (sample_id = sub('.*sample_id "([^"]*)".*', "\\1", V9),
                                 circ_id = paste0(V1, ":", V4-1, "-", V5, ":", V7))][, V9 := NULL]

            ## sumup expression by circ_id (within the same detection method and sample)
            circrnas.gtf <-
                circrnas.gtf[, .(V6 = sum(V6)),
                             by = .(V1, V2, V4, V5, V7, sample_id, circ_id)]
        }else{
            circrnas.gtf[, `:=` (sample_id = sub('.*sample_id "([^"]*)".*', "\\1", V9),
                                 circ_id = paste0(V1, ":", V4-1, "-", V5))][, V9 := NULL]

            ## sumup expression by circ_id (within the same detection method and sample)
            circrnas.gtf <-
                circrnas.gtf[, .(V6 = sum(V6)),
                             by = .(V1, V2, V4, V5, sample_id, circ_id)]
        }
        circrna.median.reads <-
            circrnas.gtf[, .(read.count = as.numeric(median(V6)),
                             n_methods = length(unique(V2)),
                             circ_methods = paste0(sort(unique(V2)),
                                                   collapse = "|")),
                         by = .(sample_id, circ_id)]

        filename <- paste0(output_prefix, "median.", file.ext)
        fwrite(x = circrna.median.reads,
               file = filename,
               row.names = F, col.names = T, sep = "\t")
    }
}
