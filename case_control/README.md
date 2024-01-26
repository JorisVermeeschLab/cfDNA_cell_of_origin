# Case-control analysis

The case-control analysis uses the correlation.csv outputs generated per sample by the gene expression analysis presented here https://github.com/shendurelab/cfDNA. The gene expression analysis should be performed using a data matrix of average normalized gene expression values summarised by cell type per tissue from the Tabula Sapiens consortium (https://datasets.cellxgene.cziscience.com/bfd80f12-725c-4482-ad7f-1ed2b4909b0d.rds). The correlation.csv outputs contain cell types ordered by strength of correlation between their gene expression values and gene-level fast-fourrier transformed window protection scores in the cfDNA sample.

**usage:** Rscript callrankdiff.R <corr_files> <case_control_comparisons> <output_dir>

**usage:** Rscript callcluster.R <my_seed> <perplexity_tnse> <initk_cluster> <walksteps_cluster> <output_dir>
