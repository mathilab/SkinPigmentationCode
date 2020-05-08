## runNullReg.sh
#` Part 2 of resampleScore.sh pipeline

source $1

cd ${wd}
## Make table of genome-wide control SNPs for each trait SNP
paste -d'\t' ${wd}/control/*.txt > controls.tsv

## Make label of trait SNPs for assignment
touch label.txt
echo "snps" > label.txt
for f in ${wd}/control/*.txt
  do
    filename=$(basename -- "$f")
    name=${filename%_*}
    awk -v snp="${name}" '{print $0 "\t" snp}' label.txt > label.tmp && mv label.tmp label.txt
  done

## Remove control folder
# rm -rf ${wd}/control/

## Make output folders for different regression models
mkdir out/anc/
mkdir out/anc/pc && mkdir out/anc/plain
if $use_anc_calls;
then
  mkdir out/anc/ef && mkdir out/anc/hg && mkdir out/anc/yam
fi
# mkdir out/anc-1kg && mkdir out/anc-1kg/pc && mkdir out/anc-1kg/plain

## Run PRS score regressions with time for control sets of SNPs
i=0
mkdir masters
mkdir sets
mkdir workspace
while read -r set
  do
    ((i+=1))
    echo $set > sets/${i}.txt
    bsub -e ${wd}/error/$i.txt -o ${wd}/stdout/$i.txt -J "set_${i}" $scripts/regressControlPRS.sh $run $i
  done < controls.tsv

## Merge summary statistics into a table
bsub -M 15000 -R -w "done(set_*) && done(null)" -J "merge" $scripts/mergeSummaryStats.sh $run
