#!/bin/sh
## Pipeline to generate time series from UK Biobank GWAS

source $1 # arguments for pipeline in format of ukbTimeSeries_arg.txt

## Filter UK Biobank data
cd $dir

## Filter UKB summary statistics by SNPs and/or p-value
wd=p_${p}_1240k_${restrict_1240k}
mkdir $wd
uksub=uk_subset_p_${p}_1240k_${restrict_1240k}.tsv
echo "Rscript --vanilla filterUkTable.R"
Rscript --vanilla $scripts/filterUkTable.R $uk $capture $p ${restrict_1240k} \
$wd/$uksub

## Clump or LD prune SNPs
#` *** USE PLINK 1.9 ***
cd $wd
#` Make list of significant hits from UKB GWAS
uksub_snps=snps_${uksub}
awk 'FNR>1 {print $1}' $uksub | awk -F ":" '{printf "%s\t%s\n", $1,$2}' >\
  $uksub_snps
#` Subset GBR VCF by significant hits from UKB
vcf_sig=${gbr_vcf%%.*}_sig.vcf.gz
echo "bcftools view -v snps -R $uksub_snps -O z -o $vcf_sig $dir/$gbr_vcf"
bcftools view -v snps -R $uksub_snps -O z -o $vcf_sig $dir/$gbr_vcf
#' Convert VCF to bed file
plink --vcf $vcf_sig --out ${vcf_sig%%.*}
#` Run plink clumping or LD pruning
c_dir=clump_${clump}_r2_$r2
mkdir $c_dir
if $clump;
  then
    # reformat UKB summary statistics for PLINK
    Rscript --vanilla $scripts/formatUKBtoClump.R $uksub $vcf_sig uk_clump.tsv
    # run clumping
    echo "plink --clump uk_clump.tsv"
    plink --clump uk_clump.tsv --clump-p1 $p --clump-r2 $r2 \
    --bfile ${vcf_sig%%.*} --out $c_dir/gbr
    # extract rsIDs
    awk '{print $3}' $c_dir/gbr.clumped > $c_dir/1kgp_gbr.clump.in
    ## Pick most significant SNPs from clumps with p-val of 0
    Rscript --vanilla $scripts/pickTopFromClump.R uk_clump.tsv $c_dir/gbr.clumped $uksub
    prune=final_clump.txt
  else
    # run LD pruning
    echo "plink --indep-pairwise 50 5 $r2"
    plink --indep-pairwise 50 5 $r2 --bfile ${vcf_sig%%.*} --out \
    $c_dir/1kgp_gbr
    prune=1kgp_gbr.prune.in
  fi

cd $c_dir

## Subset UKB summary statistics by prune.in file
ukprune=uk_pruned_r2_${r2}_1240k_${restrict_1240k}_p_${p}.tsv
echo "Rscript --vanilla subsetGwasLD.R"
Rscript --vanilla $scripts/subsetGwasLD.R ../$uksub ../$prune $ukprune 1

## Create master file from UKB SNPs
#` Also do physical pruning on the file
master=ukb_master_${out}.csv
echo "Rscript --vanilla $scripts/makeMaster.R"
Rscript --vanilla $scripts/makeMaster.R $ukprune $master ../$prune $refgene $prune_d

## Filter VCF files
now=modern_human_sub.vcf
anc=ancient_human_sub.vcf
snps=snps_master.tsv

awk -F, 'FNR>1 {printf "%s\t%s\n", $1,$2}' $master > $snps
bcftools view -v snps --max-alleles 2 -R $snps -o $now $dir/$now_vcf
bcftools view -v snps --max-alleles 2 -R $snps -o $anc $dir/$anc_vcf

## Merge VCF files (optional) and make VCF list
vcf=vcflist_${out}.txt
if $merge;
then
  ## Merge ancient and modern into one VCF
  bgzip -f $now
  bgzip -f $anc
  bcftools index modern_human_sub.vcf.gz
  bcftools index ancient_human_sub.vcf.gz
  bcftools merge -O v ancient_human_sub.vcf.gz modern_human_sub.vcf.gz > combined_${out}.vcf
  ## Make VCF list
  #` Average 2 chromosomes for modern diploid samples
  printf "%s\t1\n" combined_${out}.vcf > $vcf
else
  ## Make VCF list
  #` Diploid sample broken into two chromosomes
  vcf=vcflist_ukb.txt
  printf "%s\t0\n" $now > $vcf
  printf "%s\t1" $anc >> $vcf
fi

## Calculate PRS for ancient and modern samples
#` Score as a proportion of dark alleles
score=prs_modern_ancient_${out}.csv
echo "Calculating PRS with no weights"
python3 $scripts/getPRSHaploid.py -m $master -v $vcf -b 0 -e 6 -o $score

#` Score using effect sizes
score_b=prs_beta_modern_ancient_${out}.csv
echo "Calculating PRS with weights by beta"
python3 $scripts/getPRSHaploid.py -m $master -v $vcf -b 1 -e 6 -o $score_b

## Obtain mean PRS for 1KGP populations
avg_prs=avg_prs_1kgp.tsv
Rscript --vanilla $scripts/findAvgPRS1KGP.R $score $dir/$id1kgp $avg_prs $scripts
avg_prsb=avg_prs_beta_1kgp.tsv
Rscript --vanilla $scripts/findAvgPRS1KGP.R $score_b $dir/$id1kgp $avg_prsb $scripts

## Run regression
#` With beta coefficients
mkdir beta
cd beta
rscript --vanilla $scripts/runPRSregression.R ../$score_b $dir/$anno $dir/$pc T ../$avg_prsb $use_anc_calls $out $scripts $dir/$loc

#` Score as a proportion
cd ..
mkdir fraction
cd fraction
rscript --vanilla $scripts/runPRSregression.R ../$score $dir/$anno $dir/$pc F ../$avg_prs $use_anc_calls $out $scripts $dir/$loc

## Demography/selection plots
cd ..

## Plot ancestry vs time p-value figure and get ancestry selected SNPs
mkdir anc-time
Rscript --vanilla $scripts/plotAncTimePval.R $selscan ukb_master_${out}.csv ./anc-time/${out}

#` Reformat SNP list for extracting from selection scan
awk '{printf "%s_%s\n", $1,$2}' snps_master.tsv > snps_master_${out}.tsv
#' Extract model coefficient for each allele (estimates done before on HPC)
echo "head -n 1 $selscan > sel_snps_${out}.tsv && grep -w -f snps_master_.tsv"
head -n 1 $selscan > sel_snps_${out}.tsv && grep -w -f snps_master_${out}.tsv $selscan >> sel_snps_${out}.tsv
