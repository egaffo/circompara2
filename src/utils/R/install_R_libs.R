#!/usr/bin/env Rscript

packs <- c("optparse", "data.table", #"BH", "svglite",  
	   "plyr",
           #"scales", "reshape2", "DT",
           #"ggplot2", "ggthemes", "viridis",
           #"RSvgDevice", "rmarkdown", "knitr", "VennDiagram",
	   #"tidyr", "pheatmap", 
           #"stringi",
           "bedr")
ncpus <- as.integer(Sys.getenv('CPUS'))
if(is.na(ncpus)){
    ncpus <- 1
}
install.packages(packs, repos="https://cloud.r-project.org/", dependencies = T, Ncpus = ncpus)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#BiocManager::install(c("DESeq2", "ReportingTools", "ballgown"))
#BiocManager::install(c("ReportingTools", "ballgown"))
BiocManager::install(c("ballgown"), Ncpus = ncpus)

