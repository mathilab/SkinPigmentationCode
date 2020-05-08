#!/bin/sh
# Add header with sample names to files after split

D=$1
WD=$2

cd $WD
nfiles=$(ls | wc -l)
files=($(ls))
head -n 1 ${files[0]} > $D/header.txt

# skip first file as it has header
for i in $(eval echo "{1..$nfiles}")
do
  f=${files[i]}
  awk '{print $0}' $D/header.txt | cat - $f > /tmp/out && mv /tmp/out $f
done
