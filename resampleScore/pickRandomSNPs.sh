## pickRandomSNPs.sh
#` Part of resampleScore.sh pipeline

## Run Python script to get indices of first/last occurence of low/high SNP
module load python/3.6.3
o=$(python ${6}/binarySearch.py -s $2 -l $3 -h $4)
echo $o
low_i=$(echo $o | awk ' {print $1} ')
hi_i=$(echo $o | awk ' {print $2} ')

## Make list of SNPs
awk -v low_i=$low_i -v hi_i=$hi_i 'NR >= low_i && NR <= hi_i {print $0}' $2 | awk ' {printf "%s\n", $1} ' > ${1}.txt

## Pick {n} random SNPs from table
python ${6}/pickRandomEntry.py -n $5 -l ${1}.txt -o ${1}_random.txt

rm ${1}.txt
