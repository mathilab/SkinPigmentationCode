# SkinPigmentationCode

## Description
This repository contains code used to create the main figures for the paper
"The evolution of skin pigmentation associated variation in West Eurasia" as
well as key results. Data files that are used in some of the pipelines
are available at the link in the manuscript under "Code availability." Note
that some of the pipelines were meant to be run on a computing cluster, hence
commands for job submissions. These commands may require specific configuration
for one's own cluster to run properly.

Required software: PLINK, Python 3, R, BCFtools


## dataPrep
Steps for data processing from the downloaded Reich Lab dataset
of all published ancient individuals (v37.2), which constitutes the
capture-shotgun data referred to in the paper. The 1240K VCF of the smaller
shotgun dataset was generated from the capture-shotgun dataset but also
included additional samples not present in the release, which we pulled down
using apulldown.py and then merged. However, for manually curated SNPs in the
shotgun dataset samples, we pulled down these SNPs for all individuals since
not all SNPs are present on the 1240K capture array.

Required software: convertf, PLINK, gdc (https://github.com/mathii/gdc),
gdc3 (https://github.com/mathii/gdc3), BCFtools, R, Python 3, bedtools

Required data: 1000 Genomes, UK Biobank GWAS summary statistics

## fig1
Pipelines to create subfigures of Fig. 1. For UK Biobank GWAS skin pigment
SNPs, run ukbTimeSeries.sh from command line with arguments listed in
ukbTimeSeries_arg.txt. For manually curated SNPs, run manualTimeSeries.sh
from command line with arguments listed in manualTimeSeries_arg.txt. Also plots
results from fig2 to create subfigures of Fig. 2.

## fig2
Pipeline to generate statistics for individual SNP logistic regression. Run
runLogRegAllSNPs with arguments as shown in runLogRegAllSNPs_args.txt.

For the actual figure plots, this is generated as part of
the ukbTimeSeries.sh and manualTimeSeries.sh pipelines using the statistics
generated here in the fig1 folder.

## fig3
Pipeline to create subfigures of Fig. 3. Run plotPBSJointDiffTrees.sh with
arguments as shown in args.txt. 1000 Genomes GBR, YRI, and CHB alternate allele
frequencies must be calculated beforehand using getAlFqVCF.py and also combined
with combineFqTables.R with GBR being group 1.

## fig4
These scripts are for plotting the panels of Fig. 4. To plot (A) use
plot_figure_4A.R, (B) use plotQx.R, and (C) getUKBioPVE.R. For (A), we
ran manualTimeSeries.sh but with lists of SNPs where we manually removed top
effect size SNPs to create separate regression models. For (C), we ran
the Qx test from Berg & Coop 2014 (https://github.com/jjberg2/PolygenicAdaptationCode)
beforehand to get statistics from different sets of SNPs.

## resampleScore
Contained here is the pipeline for obtaining empirical p-values from a genomic
null distribution of frequency-matched SNPs for the genetic score time series
analysis in fig1. This pipeline was configured for running on our local
computing cluster. Run resampleScore.sh to begin the pipeline.
