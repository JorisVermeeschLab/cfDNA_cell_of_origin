library(data.table)
library(reshape2)
library(Rtsne)
library(parallel)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(igraph)
library(FNN)
library(Matrix)
library(cowplot)

outputdir <- args[1]
outprefix <- args[2]
perplexity <- args[3]
initk <- args[4]
walkstep <- args[5]
myseed <- args[6]

#pivot wider combined data file
combined.data <- fread(paste0(outputdir, "/combined_data.txt"), header=T)
combined.data <- combined.data %>% select(cell_type, sample, rank, class)
combined.data$group <- paste0(combined.data$sample, ":", combined.data$class)
