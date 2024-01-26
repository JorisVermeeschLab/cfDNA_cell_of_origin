# Case-control analysis

The case-control analysis uses the correlation.csv outputs generated per sample by the gene expression analysis presented here https://github.com/shendurelab/cfDNA. The gene expression analysis should be performed using a matrix of average gene expression values summarised by cell type per tissue from the Tabula Sapiens consortium. This gene expression matrix has been provided in data for future users. The correlation.csv outputs contain cell types ordered by strength of correlation between their gene expression values and gene-level fast-fourrier transformed window protection scores in the cfDNA sample.

**usage:** Rscript rankdiff.R 
