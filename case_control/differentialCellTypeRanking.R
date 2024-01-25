library(data.table)
library(rstatix)
library(tidyverse)
library(ggrepel)
library(rstatix)
library(ggplot2)

args <- commandArgs(TRUE)

files=args[1]
comparisons=[2]

#Read in a sample file containing all unique study sample IDs as the first column. The second column should be a path to the correlation.csv output from the https://github.com/shendurelab/cfDNA gene expression analysis for that sample. The gene expression analysis should be performed using average gene expression values summarised by cell type per tissue from the Tabula Sapiens consortium.     
