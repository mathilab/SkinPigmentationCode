## resampleScore.sh
#` Start of resampling PRS pipeline
#` Get distribution of score statistics based on logistic regression for
#` frequency matched alleles using ancient samples
#` *** RUN ON HPC CLUSTER ***

source $1

## Check for the following scripts
cd $scripts
#` Binary search and randomized selection
if test ! -f pickRandomEntry.py;
then
  echo "pickRandomEntry.py missing. Exiting..."
  exit 1
fi

if test ! -f binarySearch.py;
then
  echo "binarySearch.py missing. Exiting..."
  exit 1
fi

#` PRS calculation
if test ! -f getPRSHaploid.py ;
then
  echo "getPRSHaploid.py missing. Exiting..."
  exit 1
fi

cd $wd

## Make folders
mkdir error
mkdir stdout
mkdir warn
mkdir out
mkdir results

## Pick frequency matched SNPs genome-wide (DEPRACATED)
#` add dark allele to table and search frequency and upper/lower bound
# echo "Rscript --vanilla ${scripts}/makeDarkFqTable.R $master $all_snp_fq $err master_fq.tsv"
# bsub -M 20000 -J "makeDark" -e makedark.err Rscript --vanilla ${scripts}/makeDarkFqTable.R $master $darkfq_ref $err ${wd}/master_fq.tsv

## Continue to part 1
bsub -e err.txt -J "getRandom" ${scripts}/getRandomSNP.sh $run
