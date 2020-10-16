## admixFqCombinedTables.R
#` Make combined tables for each ancient group polarized by alternate allele
#` frequency from ADMIXTURE

admix_fq <- '' # Table with SNP ID and P matrix frequencies
group2_path <- 'chb_fq.csv' 
outgroup_path <- 'yri_fq.csv'
admix_prop <- 'admix_input.3.Q' # Q matrix from ADMIXTURE
alt_path <- '' # PLINK file with alternate alleles
a1_path <- '' # Note: recode admix_input alleles into ACTG

admix.df <- read.table(admix_fq, header=F, as.is=T)
group2.df <- read.csv(group2_path, header=T, as.is=T)
outgroup.df <- read.csv(outgroup_path, header=T, as.is=T)
prop.df <- read.table(admix_prop, as.is=T)
alt.df <- read.table(alt_path, as.is=T)
a1.df <- read.table(a1_path, as.is=T)

## Get sample sizes for ancient groups based on total ADMIXTURE proportions
#` Note: this is not necessarily integer
yam_n <- sum(prop.df$V1)
ana_n <- sum(prop.df$V3)
hg_n <- sum(prop.df$V2)

colnames(admix.df) <- c('rsid','ID','YAM','HG','ANA')
colnames(group2.df)[4:6] <- c('group2_alt','group2_n','group2_fq')
colnames(outgroup.df)[4:6] <- c('out_alt','out_n','out_fq')

outgroup_2.df <- merge(outgroup.df, group2.df, by = c('ID','chr','pos'))

# Polarize ancient frequencies from A2 allele to alternate allele
alt.df$ID <- paste(alt.df$V1, alt.df$V4, sep='_')
colnames(alt.df)[6] <- 'Alt'
a1.df$ID <- paste(a1.df$V1, a1.df$V4, sep='_')
colnames(a1.df)[6] <- 'A2'

admix.df <- merge(admix.df, alt.df[ , c('ID','Alt')], by='ID')
admix.df <- merge(admix.df, a1.df[ , c('ID','A2')], by='ID')

alt_isnt_a2 <- admix.df$A2 != admix.df$Alt

admix.df$YAM[alt_isnt_a2] <- 1 - admix.df$YAM[alt_isnt_a2]
admix.df$ANA[alt_isnt_a2] <- 1 - admix.df$ANA[alt_isnt_a2]
admix.df$HG[alt_isnt_a2] <- 1 - admix.df$HG[alt_isnt_a2]

# Export Steppe combined table
sp.df <- data.frame(ID=admix.df$ID, group1_alt=rep(NA, nrow(admix.df)), group1_n=yam_n, group1_fq=admix.df$YAM)
sp.df$group1_alt <- sp.df$group1_n * sp.df$group1_fq

all.df <- merge(outgroup_2.df, sp.df, by = 'ID')
all.df <- all.df[ , c(1,2,3,4,5,6,10,11,12,7,8,9)]

out <- 'yam_admix_fq.tsv'
write.table(all.df, file=out, quote=F, row.names=F, col.names=T, sep='\t')

# Export hunter-gatherer combined table
hg.df <- data.frame(ID=admix.df$ID, group1_alt=rep(NA, nrow(admix.df)), group1_n=hg_n, group1_fq=admix.df$HG)
hg.df$group1_alt <- hg.df$group1_n * hg.df$group1_fq

all.df <- merge(outgroup_2.df, hg.df, by = 'ID')
all.df <- all.df[ , c(1,2,3,4,5,6,10,11,12,7,8,9)]

out <- 'hg_admix_fq.tsv'
write.table(all.df, file=out, quote=F, row.names=F, col.names=T, sep='\t')

# Export Early Farmer combined table
ana.df <- data.frame(ID=admix.df$ID, group1_alt=rep(NA, nrow(admix.df)), group1_n=ana_n, group1_fq=admix.df$ANA)
ana.df$group1_alt <- ana.df$group1_n * ana.df$group1_fq

all.df <- merge(outgroup_2.df, ana.df, by = 'ID')
all.df <- all.df[ , c(1,2,3,4,5,6,10,11,12,7,8,9)]

out <- 'ana_admix_fq.tsv'
write.table(all.df, file=out, quote=F, row.names=F, col.names=T, sep='\t')