library(data.table)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggrepel)
script.dir <- "./"
functions <- file.path(script.dir,"functions.R")
source(functions)

args <- commandArgs(TRUE)

corr_files=args[1] ##A tab-separated file contains unique sample ID and location of correlation.csv 'sample\tpath'
comparisons=args[2] ##A tab-separated file contains unique sample ID and case-control status 'sample\tstatus'
output_dir=args[3] ##Output directory

#import sample file containing unique study sample IDs and the path to the correlation.csv output for each sample.
corr_files <- fread(files, header=T)

#import annotation file containing unique sample ID and the case or control status for each sample.
comparisons <- fread(comparisons, header=T)

#import correlation.csv outputs for all samples
dat <- lapply(corr_files$path, function(x) fread(x, header=T))

#extract cell type column, annotate with unique sample ID, and generate rank order.
for (i in 1:length(dat)){
        dat[[i]] <- dat[[i]][,1]
        colnames(dat[[i]]) <- c("cell_type")
        dat[[i]]$sample <- files[i,]$sample
        dat[[i]]$rank <- seq(1, nrow(dat[[i]]), by=1)
}
              
#combine output
dat <- do.call(rbind, dat)

#annotate combined output
dat.anno <- left_join(comparison, dat)
dat.anno$status <- as.factor(dat.anno$status)
write.table(dat.anno, paste0(output_dir, "/combined_data.txt"), col.names=T, row.names=F, quote=F, sep="\t")
              
#compare ranks between cases and controls for each cell type and visualize as volcano plot
plotDifferentialCellTypes(dat.anno, output_dir)        
