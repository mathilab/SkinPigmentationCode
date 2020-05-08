## runLogRegAllSNPs.sh
#` Selection scan using aDNA time series
#` Pipeline for getting logistic regression model coefficents genome-wide

source $1

cd $d
## Create alternate allele matrix for all SNPs on capture array
mkdir alt-matrix
module load python/3.4.2
echo "Making alternate allele matrix..."
python $scripts/getAltAlMat.py -v $vcf -o alt-matrix/all.csv -a $apulldown

## Split these SNPs into 5k line files on cluster
mkdir split-csv
split -l 5000 alt-matrix/all.csv split-csv/capture_all_

## Add header row with sample names to splitted CSV files
$scripts/addHeader.sh $d $d/split-csv

## Run selection scan on cluster
$scripts/runSelectScan.sh $d $anno $scripts/runLogRegAlt.R $pca

## Afterwards, run combineSelScanOut.sh and mergeScanFq.R
