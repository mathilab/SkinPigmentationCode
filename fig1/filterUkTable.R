## Filter UK Biobank summary statistics with SNPs on array

library(data.table)

args = commandArgs(trailingOnly = TRUE)
uk_path <- args[1]
filter_path <- args[2]
p <- as.numeric(args[3])
filter <- args[4]
out <- args[5]
 
# uk.df <- read.table(uk_path, header=T, as.is=T)
uk.df <- fread(uk_path, header=T, stringsAsFactors=F)

# filter by p-value
uk.df <- uk.df[which(uk.df$pval<p), ]

if (filter == T | filter == "T") {
  # filter by snp array
  filter.df <- read.table(filter_path, header=F, as.is=T)
  uk.snp <- strsplit(uk.df$variant, ":")
  uk_snp <- rep(0, length(uk.snp))
  uk_snp <- lapply(uk.snp, function(x) {
    paste(x[1],x[2],sep=':')
  })
  final.df <- uk.df[which(uk_snp %in% filter.df[ , 1]), ]
} else {
  final.df <- uk.df
}

# remove low quality imputation
final.df <- final.df[final.df$low_confidence_variant==F, ]

# save to table
write.table(final.df, file=out, row.names=F, quote=F)