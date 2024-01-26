plotDifferentialCellTypes <- function(dat, output_dir) {

  allcomp <- c()
  cell <- unique(dat$cell_type)

  #run two-sided wilcox test and compute mean foldchange for each cell type
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

plotWalkTrapClusterTsne <- function(dat, initk, walksteps, perplexity, output_dir){

  #run pca
  dtmatpca <- prcomp(dat)
  
  #get number of PCs
  eigenvalues <- dtmatpca$sdev^2
  num_PCs <- sum(eigenvalues > 1)
  cat("Number of PCs with eigenvalues > 1:", num_PCs)

  #get distance matrix and walktrap community clusters
  knn.norm <- get.knn(dtmatpca$x[,1:num_PCs],k=initk,algo="kd_tree")
  knn.norm <- data.frame(from = rep(1:nrow(knn.norm$nn.index), initk), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
  nw.norm <- graph_from_data_frame(knn.norm, directed = FALSE)
  nw.norm <- simplify(nw.norm)
  wt.norm <- cluster_walktrap(nw.norm,steps = walksteps)

  #run tsne
  tsne <- Rtsne(dtmat, check_duplicates = FALSE, pca = TRUE, initial_dims = num_PCs, perplexity=perplexity, dims=3)
  embedding <- as.data.frame(tsne$Y)
  embedding$phenotype <- rownames(dtmat)
  embedding <- embedding %>% separate(phenotype, sep=":", into=c("sample", "phenotype"))
  embedding$walktrap <- as.factor(wt.norm$membership)
  summary(embedding$walktrap)

  #plot walktrap clusters as barchart
  myclustercount <- cbind(as.data.frame(embedding$phenotype),as.data.frame(wt.norm$membership))
  colnames(myclustercount) <- c("Annotation","Cluster")
  myclustercount1 <- myclustercount %>% group_by(Annotation) %>% count(Cluster, name="Count")
  myclustercount2 <- myclustercount1 %>% group_by(Cluster) %>% add_tally(Count, name="Total")
  myclustercount3 <- myclustercount2 %>% mutate(Percentage=Count/Total)
  myclustercount3 <- as.data.frame(myclustercount3)
  myclustercount3$Cluster <- as.factor(myclustercount3$Cluster)
  myclustercount3tot <- myclustercount3[!duplicated(myclustercount3[,c("Cluster","Total")]),]
  
  bartot <- ggplot(myclustercount3tot) + geom_bar(aes(y = Total,x = Cluster), fill="lightblue",stat="identity", width= 0.8)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(),legend.position = "none",axis.title.x = element_blank())
  barclust <- ggplot(myclustercount3) + geom_bar(aes(y = Percentage, x = Cluster, fill = Annotation), stat="identity",width = 0.8)+theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank())+theme(legend.position="none") + scale_fill_viridis(option = "D", discrete = TRUE)
  barchart <- plot_grid(bartot, barclust, ncol = 1, align = "v",rel_heights=c(1,3))

  #plot tsne with walktrap clusters annotated
  p1 <- ggplot(embedding, aes(x=V1, y=V2, color=phenotype, shape=walktrap)) +
    geom_point(size=2.7) +
    guides(colour = guide_legend(override.aes = list(size=6))) +
    xlab("dim1") + ylab("dim2") +
    ggtitle(" ") +
    theme_light(base_size=13) + scale_color_viridis(option = "D", discrete = TRUE) +
    scale_shape_manual(values = 0:length(unique(wt.norm$membership)))

  #combine barchart and tsne plot
  pdf(paste0(output_dir, "/walktrap_tsne.pdf"), height=7, width=7)
  print(plot_grid(barchart, p1, ncol = 1, rel_heights=c(1,1)))
  dev.off()
}

