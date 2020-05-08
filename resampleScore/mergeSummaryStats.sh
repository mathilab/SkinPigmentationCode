## mergeSummaryStats.sh
#` Part 3 of resampleScore.sh pipeline

source $1

cd ${wd}/out

## Merge just ancient
cat anc/plain/* > ${wd}/results/anc_plain.tsv
cat anc/pc/* > ${wd}/results/anc_pc.tsv

## Merge stratified
if $use_anc_calls;
then
  cat anc/ef/* > ${wd}/results/ef_plain.tsv
  cat anc/hg/* > ${wd}/results/hg_plain.tsv
  cat anc/yam/* > ${wd}/results/yam_plain.tsv
fi

## Merge ancient + 1KGP
# cat anc-1kg/plain/* > ${wd}/results/anc-1kg_plain.tsv
# cat anc-1kg/pc/* > ${wd}/results/anc-1kg_pc.tsv

## Clean up unneeded files
cd $wd
# rm -rf masters sets workspace
