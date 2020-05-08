#!/bin/sh
## Pipeline to generate time series from manually curated SNPs

source $1 # arguments for pipeline in format of manualTimeSeries_arg.txt

cd $dir

## Filter VCF files
now=modern_human_sub.vcf
anc=ancient_human_sub.vcf
snps=snps_master.tsv

awk -F, 'FNR>1 {printf "%s\t%s\n", $1,$2}' $master > $snps
if $filter;
then
  bcftools view -v snps -R $snps -o $now $now_vcf
  bcftools view -v snps -R $snps -o $anc $anc_vcf
else
   cp $now_vcf ./${now}
   cp $anc_vcf ./${anc}
fi

## Make VCF list
vcf=vcflist.txt
if $merge;
then
  ## Merge ancient and modern into one VCF
  bgzip -f $now
  bgzip -f $anc
  bcftools index modern_human_sub.vcf.gz
  bcftools index ancient_human_sub.vcf.gz
  bcftools merge -O v ancient_human_sub.vcf.gz modern_human_sub.vcf.gz > combined.vcf
  ## Make VCF list
  #` Average 2 chromosomes for modern diploid samples
  printf "%s\t1\n" combined.vcf > $vcf
else
  ## Make VCF list
  #` Diploid sample broken into two chromosomes
  vcf=vcflist_ukb.txt
  printf "%s\t0\n" $now > $vcf
  printf "%s\t1" $anc >> $vcf
fi

## Calculate PRS for ancient and modern samples
#` Score as a proportion of dark alleles
score=prs_modern_ancient.csv
python3 $scripts/getPRSHaploid.py -m $master -v $vcf -b 0 -e 6 -o $score

#` Score using effect sizes
if $weight;
then
score_b=prs_beta_modern_ancient.csv
python3 $scripts/getPRSHaploid.py -m $master -v $vcf -b 1 -e 6 -o $score_b
fi

## Obtain mean PRS for 1KGP populations
avg_prs=avg_prs_1kgp.tsv
Rscript --vanilla $scripts/findAvgPRS1KGP.R $score $id1kgp $avg_prs $scripts

if $weight;
then
avg_prsb=avg_prs_beta_1kgp.tsv
Rscript --vanilla $scripts/findAvgPRS1KGP.R $score_b $id1kgp $avg_prsb $scripts
fi

## Run regression for modern and ancient populations
#` With beta coefficients
if $weight;
then
mkdir beta
cd beta
rscript --vanilla $scripts/runPRSregression.R ../$score_b $anno $pc \
T ../$avg_prsb $use_anc_calls $out $scripts $loc
fi

#` Score as a proportion
cd $dir
mkdir fraction
cd fraction
rscript --vanilla $scripts/runPRSregression.R ../$score $anno $pc F ../$avg_prs $use_anc_calls $out ${scripts} $loc

## Demography/selection plots
cd ..
#` Reformat SNP list for extracting from selection scan
awk '{printf "%s_%s\n", $1,$2}' snps_master.tsv > snps_master_.tsv
#' Extract model coefficient for each allele (estimates done before on HPC)
echo "head -n 1 $selscan > sel_snps.tsv && grep -w -f snps_master_.tsv"
head -n 1 $selscan > sel_snps.tsv && grep -w -f snps_master_.tsv \
  $selscan >> sel_snps.tsv
#` Plot the p values of ancestry against time
mkdir anc-time && cd anc-time
echo "Rscript --vanilla plotAncestryTime.R"
Rscript --vanilla $scripts/plotAncTimePval.R $selscan $master anc_time
