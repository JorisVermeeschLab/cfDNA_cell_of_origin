plotDifferentialCellTypes <- function(dat, output_dir) {
  
  allcomp <- c()
  cell <- unique(dat$cell_type)

  #run wilcox and compute mean foldchange for each cell type
  for(i in 1:length(cell)) {
    df <- subset(dat, cell_type %in% cell[i])
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
    geom_label_repel(size=2) +
    scale_color_manual(values = mycolors) +
    theme_classic() + 
    ggtitle(labels)
  ggsave(paste0(output_dir, "/", labels, ".pdf"), height=7, width=7)
}
