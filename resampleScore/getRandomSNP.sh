## getRandomSNP.sh
#` Part 1 of resampleScore.sh pipeline

source $1

cd $wd
## Get list of SNPs around the frequency of a variant
awk ' FNR>1 {print $1} ' ${master_fq} > snps.tsv

mkdir control && cd control
mkdir pick-err
while read -r snp # get list of {n} randomized control SNPs per trait SNP
  do
    fq_line=$(grep -w $snp -m1 ${master_fq})
    low=$(awk '{print $4}' <<< $fq_line)
    hi=$(awk '{print $5}' <<< $fq_line)
    bsub -e pick-err/${snp}.txt -J "snp_${snp}" ${scripts}/pickRandomSNPs.sh $snp $snp_fq $low $hi $n $scripts
  done < ${wd}/snps.tsv

## Continue to part 2
bsub -w "done(snp*) && done(getRandom)" -J "null" ${scripts}/runNullReg.sh $run
