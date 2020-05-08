## Convert UK table to format for PLINK clumping
#` Requires VCF with ID column

args = commandArgs(trailingOnly=TRUE) 
uk_path <- args[1]
vcf_path <- args[2]
out_path <- args[3]

uk.df <- read.table(uk_path, stringsAsFactors=F, header=T)
vcf.df <- read.table(vcf_path, as.is=T, header=F)

# add chr and pos to uk df
variant <- strsplit(uk.df[,1], split = ":")
uk.df$id <- lapply(variant, function(x){
  paste(x[1],x[2],sep=':')})
uk.df$id <- unlist(uk.df$id)
chr_pos <- strsplit(uk.df$id, ":")
uk.df$chr <- lapply(chr_pos, function(x){x[1]})
uk.df$pos <- lapply(chr_pos, function(x){x[2]})
uk.df$chr <- unlist(uk.df$chr)
uk.df$pos <- unlist(uk.df$pos)
uk.df$chr <- as.integer(uk.df$chr)
uk.df$pos <- as.integer(uk.df$pos)

# link genomic coordinates to rsID  
colnames(vcf.df)[1:3] <- c('chr','pos','SNP')
uk.vcf.df <- merge(uk.df, vcf.df[ , 1:3], by = c('chr', 'pos'))

# write table of UKB summary statistics for clumping
final.df <- data.frame(uk.vcf.df$chr, 
                       uk.vcf.df$pos, 
                       uk.vcf.df$SNP, 
                       uk.vcf.df$pval)
colnames(final.df)[1:4] <- c('chr','pos','SNP','P')

write.table(final.df, file=out_path, row.names=F, quote=F)