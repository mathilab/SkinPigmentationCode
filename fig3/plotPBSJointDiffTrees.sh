#!/bin/sh
## plotPBSJointDiffTrees.sh
#` Plot both distributions of PBS

source $1

cd $wd

echo "Calculating alt allele frequencies for ancient group..."
python3 $scripts/getAlFqVCF.py -v $vcf -o anc_fq.csv -d F

## Filter sites by alleles sampled in ancient
Rscript --vanilla $scripts/filterByFq.R anc_fq.csv $filter_snps_th filt_anc_fq.csv

## Combine results into one table for PBS calculation
echo "Combining tables..."
Rscript --vanilla $scripts/combineFqTables.R filt_anc_fq.csv $eas_fq $afr_fq anc_comb_fq.tsv

## Filter comparison EUR-EAS-AFR fq table by sites filtered in ancient
# awk 'FNR>1 {print $1}' anc_comb_fq.tsv > filter_sites.txt
# Rscript --vanilla $scripts/filterCombinedFqTable.R filter_sites.txt $eur_comb_table filt_eur_table.csv

## Calculate Fst and PBS over windows of [X] SNPs ##
echo "Calculating present-day PBS genome-wide..."
Rscript --vanilla $scripts/calcPBSJointly.R anc_comb_fq.tsv $eur_comb_table $window pbs_table.tsv

## Get percentiles of trait SNPs based on PBS of windows
echo "Calculating PBS for trait SNPs..."
Rscript --vanilla $scripts/getPBSPercentJoint.R all_fq.tsv $window $trait_snp pbs_table.tsv trait_pbs.tsv

## Make BED file for all PBS windows
awk 'FNR>1 {printf "chr%s\t%s\t%s\n", $1,$2,$3}' pbs_table.tsv > pbs_all.bed

## Get nearest gene for top hits in the graph
bedtools sort -i pbs_all.bed > pbs_all_sort.bed
bedtools closest -t first -a pbs_all_sort.bed -b $refgene > closest.tsv

## Plot joint distribution of present-day EUR and ancient PBS
#` Do this in the Rscript plotJoint2PBStrees.R
