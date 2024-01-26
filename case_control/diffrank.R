library(data.table)
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggrepel)
script.dir <- "./"
functions <- file.path(script.dir,"functions.R")
source(functions)

args <- commandArgs(TRUE)

corr_files=args[1] 
comparisons=[2]
output_dir=[3]

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
              
#compare ranks between cases and controls for each cell type
allcomp <- c()
cell <- unique(dat.anno$cell_type)
           
for(i in 1:length(cell)) {
  df <- subset(dat.anno, cell_type %in% cell[i])
  pheno.counts <- as.data.table(table(df$status))
  labels <- paste0(pheno.counts[1,1],".n.", pheno.counts[1,2], ".vrs.", pheno.counts[2,1], ".n.", pheno.counts[2,2])
  print(labels)
  print(paste0("reference group = ", pheno.counts[2,1]$V1))
  df.wilcox <- wilcox_test(df, rank ~ status, ref.group = pheno.counts[2,1]$V1)
  df.mean <- df %>% group_by(status) %>% summarize(Mean = mean(rank, na.rm=TRUE)) %>% as.data.frame()
  df.mean.FC <- df.mean$Mean[2]/df.mean$Mean[1] %>% as.data.frame()
  colnames(df.mean.FC) <- c("foldchange")
  df.mean.FC$cell_type <- cell[i]
  df.all <- cbind(df.wilcox, df.mean.FC)
  allcomp <- rbind(allcomp, df.all)
}

#multiple testing correction
allcomp$p.adj <- p.adjust(allcomp$p, method="fdr")

#order by p.adj
allcomp <- allcomp[order(allcomp[,10]),]

#export output
write.table(allcomp, paste0(output_dir, "/", labels, ".txt"), sep="\t", col.names=T, row.names=F, quote=F)

#create color and labels for volcano plot
allcomp$diffranked <- "Not Sig"
allcomp$diffranked[allcomp$foldchange > 1 & allcomp$p.adj < 0.05] <- "Up"
allcomp$diffranked[allcomp$foldchange < 1 & allcomp$p.adj < 0.05] <- "Down"

allcomp$delabel <- NA
allcomp[allcomp$diffranked == "Down" | allcomp$diffranked == "Up",]$delabel <- allcomp[allcomp$diffranked == "Down" | allcomp$diffranked == "Up",]$cell_type
       
#set colors for volcano plot
mycolors <- c("#541352FF", "#10a53dFF", "grey80")
names(mycolors) <- c("Down", "Up", "Not Sig") 

#volcano plot
ggplot(data=allcomp, aes(x=log2(foldchange), y=-log10(p), color=diffranked, label=delabel)) +
    geom_point() +
    theme_minimal() +
    geom_label_repel(size=2)+
    scale_color_manual(values = mycolors)+
    theme_classic() + ggtitle(labels)
ggsave(paste0(output_dir, "/", labels, ".pdf"), height=7, width=7)
              
