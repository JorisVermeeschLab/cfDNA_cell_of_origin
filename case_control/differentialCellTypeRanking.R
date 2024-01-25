library(data.table)
library(rstatix)
library(tidyverse)
library(ggrepel)
library(rstatix)
library(ggplot2)

args <- commandArgs(TRUE)

files=args[1]
comparisons=[2]
output_data=[3]
output_plot=[4]

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

#annotate combined output
dat.anno <- left_join(comparison, dat)
dat.anno$status <- as.factor(dat.anno$status)

#compare ranks between cases and controls for each cell type
allcomp <- c()
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
write.table(alltruegroup, paste0(output_data, "/", labels, ".txt"), sep="\t", col.names=T, row.names=F, quote=F)

#create data label for volcano plot
allcomp$diffexpressed <- "Not Sig"
allcomp$diffexpressed[allcomp$foldchange > 1 & allcomp$p.adj < 0.05] <- "Up"
allcomp$diffexpressed[allcomp$foldchange < 1 & alltruegroup$p.adj < 0.05] <- "Down"

#set colors for volcano plot
mycolors <- c("#541352FF", "#10a53dFF", "grey80")
names(mycolors) <- c("Down", "Up", "Not Sig") 

#volcano plot
ggplot(data=allcomp, aes(x=log2(foldchange), y=-log10(p), color=diffexpressed)) +
    geom_point() +
    theme_minimal() +
    geom_label_repel(data=filter(allcomp, p.adj<0.05), aes(label=cell_type),size=0.05, label.padding = 0.05)+
    scale_color_manual(values = mycolors)+
    theme_classic() + ggtitle(labels)
ggsave(paste0(output_plot, "/", labels, ".pdf"), height=7, width=7)
              
