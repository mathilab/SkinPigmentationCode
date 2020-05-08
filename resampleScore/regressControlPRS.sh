## regressControlPRS.sh
#` Part of the resampleScore.sh pipeline
#` Run PRS ancient DNA time regressions

source $1 # arguments for pipeline in format of regressControlPRS.sh
module load R/3.2.1
module load python/3.6.3

cd $wd/workspace
mkdir $2 && cd $2

## make randomized control master snp annotation file
Rscript --vanilla $scripts/makeRandomMaster.R ${wd}/sets/${2}.txt $master $all_snp_fq ${wd}/label.txt ${wd}/masters/${2}.csv $use_b

## Filter VCF files
now=modern_human_sub.vcf
anc=ancient_human_sub.vcf
snps=snps_master.tsv

awk -F, 'FNR>1 {printf "%s\t%s\n", $1,$2}' ${wd}/masters/${2}.csv > $snps
bcftools view -v snps -m2 -M2 -R $snps -o $now $now_vcf
bcftools view -v snps -m2 -M2 -R $snps -o $anc $anc_vcf

## Make VCF list
vcf=vcflist.txt
touch $vcf
printf "%s\t0\n" $now >> $vcf
printf "%s\t1" $anc >> $vcf

## Calculate PRS for ancient and modern samples
#` Score as a proportion of dark alleles
score=prs_modern_ancient.csv
python3 $scripts/getPRSHaploid.py -m ${wd}/masters/${2}.csv -v $vcf -b $use_b -e 6 -o $score

## Run regression for modern and ancient populations
#` Polygenic Score
Rscript --vanilla $scripts/runPRSControl.R $score $anno $pc $use_b $loc ${wd}/out/ $2 $scripts $use_anc_calls

## Delete workspace and sets
# rm -rf $wd/workspace/$2 ${wd}/sets/${2}.txt
