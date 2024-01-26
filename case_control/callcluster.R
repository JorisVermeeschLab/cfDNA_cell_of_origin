library(data.table)
library(Rtsne)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(igraph)
library(FNN)
library(cowplot)
script.dir <- "./"
functions <- file.path(script.dir,"functions.R")
source(functions)

args <- commandArgs(TRUE)

output_dir <- args[1]
perplexity <- args[3]
initk <- args[4]
walksteps <- args[5]
myseed <- args[6]

#pivot wider combined data file
combined <- fread(paste0(output_dir, "/combined_data.txt"), header=T)
combined$group <- paste0(combined$sample, ":", combined$status)
combined.wider <- combined %>% select(-sample, -status) %>% pivot_wider(names_from = cell_type, values_from = rank) %>% as.data.frame()
rownames(combined.wider) <- combined.wider$group
combined.wider <- combined.wider %>% select(-group) %>% data.matrix()

set.seed(myseed)
plotWalkTrapClusterTsne(combined.wider, initk, walksteps, perplexity, output_dir)