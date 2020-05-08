# Merge tables from selection scan and allele frequencies for modern European

args = commandArgs(trailingOnly = TRUE)
# arg1 is output of selections can
# arg2 contains frequencies
# arg3 is path for output file

scan.df = read.table(args[1], as.is=T, header=T)
fq.df = read.table(args[2], as.is=T)

colnames(fq.df) <- c('SNP', 'Fq')
final.df <- merge(scan.df, fq.df, by = 'SNP')
write.table(final.df, args[3], row.names=F, sep="\t", quote=F)