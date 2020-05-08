## makeDarkFqTable.R
#` Outputs table with dark, derived alleles and frequencies of SNPs from master
#` Run this as input to resampleScore.sh pipeline for UKB/manual SNPs

args = commandArgs(trailingOnly = TRUE)

master_path <- args[1]
fq_path <- args[2]
err <- as.numeric(args[3])
out_path <- args[4]

master.df <- read.csv(file=master_path, header=T, as.is=T)
fq.df <- read.csv(file=fq_path, header=T, as.is=T)

master.df$ID <- paste(master.df$CHROM, master.df$POS, sep='_')

# add derived fq to trait SNP table 
df <- merge(master.df, fq.df, by='ID')
i <- which(df$Alt.x == df$Der)
j <- which(df$Alt.x != df$Der)

df$derfq <- NA
df$derfq[i] <- as.numeric(df$Fq)[i]
df$derfq[j] <- 1 - as.numeric(df$Fq)[j]

# make lower and upper search fq variables
df$high <- df$derfq + err 
df$low <- df$derfq - err 

# export table
final.df <- df[ , c('ID','Dark','Der','low','high')]
write.table(final.df, file=out_path, quote=F, row.names=F, sep='\t')