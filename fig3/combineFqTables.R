## combingroup2qTables.R
#` Combine the tables for groups' alt allele frequencies for PBS analysis

args = commandArgs(trailingOnly = TRUE)
group1_path <- args[1]
group2_path <- args[2]
outgroup_path <- args[3]
out_path <- args[4]

group1.df <- read.csv(group1_path, header=T, as.is=T)
group2.df <- read.csv(group2_path, header=T, as.is=T)
outgroup.df <- read.csv(outgroup_path, header=T, as.is=T)

colnames(group1.df)[4:6] <- c('group1_alt','group1_n','group1_fq')
colnames(group2.df)[4:6] <- c('group2_alt','group2_n','group2_fq')
colnames(outgroup.df)[4:6] <- c('out_alt','out_n','out_fq')

outgroup1.df <- merge(outgroup.df, group1.df, by = c('ID','chr','pos'))
all.df <- merge(outgroup1.df, group2.df, by = c('ID','chr','pos'))

write.table(all.df, file=out_path, quote=F, row.names=F, col.names=T, sep='\t')
