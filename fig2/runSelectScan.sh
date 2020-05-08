#!/bin/sh
# Generate beta and p values across genome from aDNA samples
# for all capture array SNPs along with alternate allele frequencies
# RUN ON CLUSTER

WD=$1
ANNO_PATH=$2
RSCRIPT=$3
PCA=$4

cd $WD
mkdir run
mkdir warn
mkdir log
mkdir error

module load R/3.5.1
for f in ${WD}/split-csv/*
do
  bsub -e error/error_${f##*/}.txt -o log/log_${f##*/}.txt "Rscript --vanilla $RSCRIPT $f $ANNO_PATH run/${f##*/}.tsv warn/${f##*/}.tsv $PCA"
done
