## Subset GWAS results by Plink LD pruned files

args = commandArgs(trailingOnly = TRUE)
gwas_path <- args[1]
prune_path <- args[2]
out_path <- args[3]
is_uk <- args[4] # if gwas file is UK Biobank

library(biomaRt)

snps.in <- read.delim(prune_path, header=F, as.is=T)
gwas.table <- read.table(gwas_path, as.is=T, header=T)

# convert rsID to chr:pos
hsapiens.grch37 <- useEnsembl(biomart = "snp", dataset="hsapiens_snp", GRCh=37)
snp.loc <- getBM(attributes = c("chr_name", "chrom_start"), 
                 filters = "snp_filter", 
                 values = snps.in$V1,
                 mart = hsapiens.grch37)
colnames(snp.loc) <- c('Chr','bp')
snp.loc$id <- paste(snp.loc$Chr, snp.loc$bp, sep=':')

# subset gwas table by pruned in SNPs
if(is_uk %in% c(1,'yes','y')) {
  uk.snp <- strsplit(gwas.table[ , 1], ":")
  uk_snp <- rep(0, length(uk.snp))
  uk_snp <- lapply(uk.snp, function(x){paste(x[1],x[2],sep=':')})
  final.df <- gwas.table[which(uk_snp %in% snp.loc$id), ]
} else {
  final.df <- gwas.table[which(gwas.table[ , 1] %in% snp.loc$id), ]
}

# remove indels from table
uk.df <- final.df
variant <- strsplit(uk.df[,1], split = ":")
uk.df$a1 <- lapply(variant, function(x){
  x[3]})
uk.df$a2 <- lapply(variant, function(x){
  x[4]})
uk.df$a1 <- unlist(uk.df$a1)
uk.df$a2 <- unlist(uk.df$a2)

remove <- nchar(uk.df$a1)>1 | nchar(uk.df$a2)>1
final.df <- final.df[!(remove), ]

write.table(final.df, file=out_path, sep="\t", row.names=F, quote=F)
