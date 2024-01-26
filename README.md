# Cell type signatures in cell free DNA fragmentation profiles reveal disease biology

The scripts provide details for the analyses performed in the manuscript entitled 'Cell type signatures in cell free DNA fragmentation profiles reveal disease biology'.

> This is not the release of a software package. All scripts and binaries are provided as is, without any warranty. We are only providing this information and code in addition to the description of methods to make it easier to reproduce the analyses.

## Versions
```
R version 4.1.3
```

## Bioinformatic processing
Sequencing data was demultiplexed, quality checked, and adapters were trimmed using fastp (version 0.12.4). Raw reads were aligned to the human reference genome GRCh38 using the Burrows-Wheeler aligner (version 0.7.17).

## Data availability
Whole genome sequencing data can be downloaded from the European Genome-Phenome Archive under EGAC00001003120. Nucleosome calls for study samples can be found in data.

## Analysis
For case-control related analyses, please see the folder case_control
usage: Rscript callCluster.R <test_demo_sample_List> <test_demo_label_List> <gipseqcount_option_input> all 1 noscale <testprefix_out> yes <PCs_number> <k_neareast_neighbors_number> <walkstep_number> <seed_number>

usage: Rscript callPCAmodel.R <test_demo_sample_List> <test_demo_label_List> <gipseqcount_option_input> all 1 noscale <PCs_number> <performance_output> <model_dir> <model_name> <test_demo_predict_sample_List> <test_demo_predict_label_List> <prediction_output>

For prediction related analyses, please see the folder classification

usage: Rscript callPCAmodel.R <test_demo_sample_List> <test_demo_label_List> <gipseqcount_option_input> all 1 noscale <PCs_number> <performance_output> <model_dir> <model_name> <test_demo_predict_sample_List> <test_demo_predict_label_List> <prediction_output>
