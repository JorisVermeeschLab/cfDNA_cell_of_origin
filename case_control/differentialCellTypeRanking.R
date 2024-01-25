library(data.table)
library(rstatix)
library(tidyverse)
library(ggrepel)
library(rstatix)
library(ggplot2)

args <- commandArgs(TRUE)

files=args[1]
comparisons=[2]

#import sample file containing unique study sample IDs as the first column ("sample") and the path to the correlation.csv output for each sample in the second column ("path").
files <- fread(files, header=T)

#import annotation file containing unique sample ID in the first column ("sample") and case-control status in the second column ("status").
comparisons <- fread(comparisons, header=T)

#import correlation.csv outputs for all samples
dat <- lapply(files$path, function(x) fread(x, header=T))

#extract cell type column, annotate with unique sample ID, and generate rank order.
for (i in 1:length(dat)){
        dat[[i]] <- dat[[i]][,1]
        colnames(dat[[i]]) <- c("cell_type")
        dat[[i]]$sample <- files[i,]$sample
        dat[[i]]$rank <- seq(1, nrow(dat[[i]]), by=1)
}
              
#combine output
dat <- do.call(rbind, dat)
dat.anno <- left_join(comparison, dat)
dat.anno$status <- as.factor(dat.anno$status)

#compare ranks between cases and controls for each cell type
alltruegroup <- c()
cell <- unique(dat.anno$cell_type)
           
for(i in 1:length(cell)) {
  df <- subset(dat.anno, cell_type %in% cell[i])
  pheno.counts <- as.data.table(table(df$status))
  labels <- paste0(pheno.counts[1,1],".n.", pheno.counts[1,2], ".vrs.", pheno.counts[2,1], ".n.", pheno.counts[2,2])
  print(labels)
  print(paste0("reference group = ", pheno.counts[2,1]$V1))
  df.wilcox <- wilcox_test(df, rank ~ status, ref.group = pheno.counts[2,1]$V1, p.adjust.method="fdr")
  df.mean <- df %>% group_by(status) %>% summarize(Mean = mean(rank, na.rm=TRUE)) %>% as.data.frame()
  df.mean.FC <- df.mean$Mean[2]/df.mean$Mean[1] %>% as.data.frame()
  print(paste0("foldchange control: ", df.mean[2,1]))
  colnames(df.mean.FC) <- c("foldchange")
  #df.mean.FC$cell_type_tissue <- cell[i]
  df.mean.FC$cell_type <- cell[i]
  df.all <- cbind(df.wilcox, df.mean.FC)
  alltruegroup <- rbind(alltruegroup, df.all)
}
  alltruegroup$p.adj <- p.adjust(alltruegroup$p, method="fdr")
  alltruegroup <- alltruegroup[order(alltruegroup[,10]),]
  #alltruegroup <- left_join(alltruegroup, comp)
  alltruegroup$diffexpressed <- "Not Sig"
  alltruegroup$diffexpressed[alltruegroup$foldchange > 1 & alltruegroup$p.adj < 0.05] <- "Up"
