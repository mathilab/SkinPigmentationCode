## makeMaster.R
## Make master file for UKB SNPs

args = commandArgs(trailingOnly = TRUE)
uk_path <- args[1]
out_path <- args[2]
rsid_path <- args[3]
gene_path <- args[4]
prune_d <- as.numeric(args[5])

# set.seed(42)
library(biomaRt)

# read data
uk.df = read.table(uk_path, as.is=T, header=T)
snps.in <- read.delim(rsid_path, header=T, as.is=T)

# get chr:pos, ref, alt allele
variant <- strsplit(uk.df[ , 1], ":")
uk.df$ref <- lapply(variant, function(x){x[3]})
uk.df$alt <- lapply(variant, function(x){x[4]})
uk.df$id <- lapply(variant, function(x){
                paste(x[1],x[2],sep=':')})
uk.df$ref <- unlist(uk.df$ref)
uk.df$alt <- unlist(uk.df$alt)
uk.df$id <- unlist(uk.df$id)

# get rsID, ancestral allele
hsapiens.grch37 <- useEnsembl(biomart = "snp", dataset="hsapiens_snp", GRCh=37)
bm_att <- c("chr_name", "chrom_start", "refsnp_id", "allele_1")
snp.bm <- getBM(attributes = bm_att, 
                 filters = "snp_filter", 
                 values = snps.in[ , 1],
                 mart = hsapiens.grch37)
colnames(snp.bm) <- c('chr','pos','rsid','anc')
snp.bm$id <- paste(snp.bm$chr, snp.bm$pos, sep=':')

# get dark/light allele
uk.df$dark <- NA
uk.df$light <- NA

uk.df$dark[which(uk.df$beta>0)] <- uk.df$alt[which(uk.df$beta>0)]  
uk.df$dark[which(uk.df$beta<0)] <- uk.df$ref[which(uk.df$beta<0)] 
uk.df$light[which(uk.df$beta>0)] <- uk.df$ref[which(uk.df$beta>0)]  
uk.df$light[which(uk.df$beta<0)] <- uk.df$alt[which(uk.df$beta<0)]

# create output data frame
uk.df <- merge(uk.df, snp.bm, by='id')
uk.df$der <- NA
uk.df$major <- NA
final.df <- data.frame(uk.df$chr,uk.df$pos,uk.df$rsid,uk.df$dark,uk.df$light,
                       uk.df$ref,uk.df$alt,uk.df$anc,uk.df$der, 
                       uk.df$minor_allele,uk.df$major,uk.df$beta,uk.df$se,
                       uk.df$pval, uk.df$tstat,stringsAsFactors=F)

# fill in miscellaneous info
final.df$gene <- NA 
final.df$pop <- 'EUR'
final.df$source <- 'UKBiobank'

# rename columns
col_names <- c('CHROM', 'POS', 'RSID', 'Dark', 'Light','Ref', 'Alt', 'Anc', 
               'Der', 'Minor', 'Major', 'Beta', 'SE',	'p', 'tstat', 'NearestGene', 'Pop','Source')
colnames(final.df) <- col_names

# sort data frame
final.df <- final.df[with(final.df, order(CHROM, POS)), ]

# find nearest gene
bed.df <- data.frame(final.df$CHROM, final.df$POS, final.df$POS)
bed.df$final.df.CHROM <- paste('chr', bed.df$final.df.CHROM, sep='')
write.table(bed.df, file='ukb_snps.bed', row.names=F, col.names=F, quote=F, sep='\t')
system("bedtools sort -i ukb_snps.bed > ukb_snps_sorted.bed")

bed_command <- paste('bedtools closest -t first -a ukb_snps_sorted.bed -b', gene_path, '> nearest.tsv', sep=' ')
system(bed_command)

near.df <- read.table('nearest.tsv', as.is=T)
final.df$NearestGene <- near.df$V7

# fill in derived allele
x <- which(final.df$Anc != final.df$Ref)
final.df$Der[x] <- final.df$Ref[x]
y <- which(final.df$Anc == final.df$Ref)
final.df$Der[y] <- final.df$Alt[y]

# fill in major allele 
x <- which(final.df$Minor != final.df$Ref)
final.df$Major[x] <- final.df$Ref[x]
y <- which(final.df$Minor == final.df$Ref)
final.df$Major[y] <- final.df$Alt[y]

# take absolute value of tstat
final.df$abs_t <- abs(final.df$tstat)
                      
# Perform physical pruning 
ids <- c() # which snps to keep 
for (i in 1:22) {
  # subset by chr
  chr.df <- final.df[final.df$CHROM==i, c('POS','abs_t','RSID')]
  # skip to next iteration if chr not present
  if (nrow(chr.df)==0) {
    next
  }
  # order by t-stat and pos ascending
  chr.df <- chr.df[order(chr.df$abs_t, chr.df$POS, decreasing=T), ]
  # iteratively remove nearby SNPs from top to bottom 
  for (j in 1:nrow(chr.df)) {
    # stop pruning chromosome at end
    if (j==nrow(chr.df)) break
    cur_pos <- chr.df$POS[j]
    # remove nearby SNPs
    keep <- chr.df$POS<(cur_pos-prune_d) | chr.df$POS>(cur_pos+prune_d)
    keep[j] <- TRUE
    chr.df <- chr.df[keep, ]
  }
  # record to big list
  ids <- c(ids, chr.df$RSID)
}
  
  
# subset df by physical pruning list 
final.df <- final.df[final.df$RSID %in% ids, -c(19,15)]

write.table(final.df, out_path, sep=',', quote=F, row.names=F)
