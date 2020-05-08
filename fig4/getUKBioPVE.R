## Calculate proportion of variance explained for each SNP from UK Biobank

library(data.table)
library(scales)
library(ggplot2)
library(ggrepel)

# Original UKB GWAS summary statistics TSV file
uk_path <- '1717.gwas.imputed_v3.both_sexes.tsv'

# CSV of UK Biobank skin color SNPs in master file layout
master_path <- ''

uk.df <- fread(uk_path, stringsAsFactors=F)
master.df <- read.csv(master_path, as.is=T)

# Parse through UK ID
variant <- strsplit(uk.df$variant, ":")
uk.df$ref <- lapply(variant, function(x){x[3]})
uk.df$alt <- lapply(variant, function(x){x[4]})
uk.df$id <- lapply(variant, function(x){
  paste(x[1],x[2],sep=':')})
uk.df$ref <- unlist(uk.df$ref)
uk.df$alt <- unlist(uk.df$alt)
uk.df$id <- unlist(uk.df$id)

master.df$id <- paste(master.df$CHROM, master.df$POS, sep=':')
# Add MAF to master table
pve.df <- merge(master.df, uk.df[,c('id','minor_AF','se','minor_allele')], by='id')

# Calculate PVE
pve.df$ve <- 2 * pve.df$Beta^2 * pve.df$minor_AF * (1-pve.df$minor_AF)
pve.df$pve <- pve.df$ve / (pve.df$ve + 2 * 356530 * pve.df$se^2 * pve.df$minor_AF * (1-pve.df$minor_AF))


 
## Calculate PVE between populations
# TSV file with alternate allele frequencies for YRI 
yri_fq_path <- ''

yri.df <- read.table(yri_fq_path, as.is=T)

colnames(yri.df)[c(1,5)] <- c('id','yri_alt_fq')

pve.df$id <- paste(pve.df$CHROM, pve.df$POS, sep='_')

# add YRI alt allele frequencies
pve.df <- merge(pve.df, yri.df[,c(1,5)], by='id')

# get YRI fq based on UKB minor allele
pve.df$yri_fq <- NA
i <- pve.df$Alt==pve.df$minor_allele
j <- pve.df$Alt!=pve.df$minor_allele

pve.df$yri_fq[i] <- pve.df$yri_alt_fq[i]
pve.df$yri_fq[j] <- 1 - pve.df$yri_alt_fq[j]
 
# absolute variation b/w populations
pve.df$ve_yri_uk <- pve.df$Beta^2 * ((0.5 * (pve.df$minor_AF + pve.df$yri_fq)) * (1 - (0.5 * (pve.df$minor_AF + pve.df$yri_fq))) - 0.5 * pve.df$minor_AF * (1-pve.df$minor_AF) - 0.5 * pve.df$yri_fq * (1-pve.df$yri_fq)) 




ggplot(data=pve.df, aes(x=ve, y=ve_yri_uk)) + 
  geom_point(color="black", alpha=0.65, pch=16, aes(size=abs(Beta))) +
  theme_classic() + scale_size_area(max_size = 5) +
  xlab('VE within UK') + ylab('VE between UK and YRI') +
  theme(legend.position="none") +
  geom_label_repel(aes(label=NearestGene), size=3, alpha=0.7, nudge_y = 0.002, nudge_x = 0.001, force = 50,
                   data=pve.df[pve.df$ve_yri_uk>0.001 | abs(pve.df$Beta)>0.2, ],
                   min.segment.length = 0.015,
                   segment.alpha = 0.7)

out_4b <- ''
ggsave(out_4b,width=3,height=3)  

# write.table(pve.df, file='pve.csv', row.names=F, col.names=T, quote=F, sep=',')
