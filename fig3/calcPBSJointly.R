## calcPBS.R
#` Takes table from combineFqTables.R and calculates PBS for two different trees
#` Calculate PBS for windows of{X} SNPs

library(data.table)

calcFst <- function(p1,p2,n1,n2) {
  num <- (p1-p2)^2 - ( (p1*(1-p1))/(n1-1) - (p2*(1-p2))/(n2-1) )
  den <- (p1*(1-p2)) + (p2*(1-p1))
  return(data.frame(num, den))
}

calcAvgFst <- function(p1,p2,n1,n2) {
  nd.df <- calcFst(p1,p2,n1,n2)
  avg_num <- mean(nd.df[, 1], na.rm=T) 
  avg_den <- mean(nd.df[, 2], na.rm=T)
  if (avg_den == 0) {
    fst <- 0} else {
      fst <- avg_num / avg_den
    }
  return(fst)
}

args = commandArgs(trailingOnly = TRUE)

df_path1 <- args[1]
df_path2 <- args[2]
window <- as.integer(args[3])
out_path <- args[4]
  
df1 <- fread(df_path1, stringsAsFactors = F)
df2 <- fread(df_path2, stringsAsFactors = F)
half_w <- window / 2

# Rename columns 
colnames(df1)[c(7,8,9)] <- c('anc_alt','anc_n','anc_fq')
colnames(df2)[c(7,8,9)] <- c('eur_alt','eur_n','eur_fq')
  
# Merge both data frames 
df <- merge(df1, df2[,c(1,2,3,7,8,9)], by=c('ID','chr','pos'))

# Remove duplicate rows
df <- unique(df)

# Write table with all allele frequencies
write.table(df, 'all_fq.tsv', row.names=F, col.names=T, quote=F, sep='\t')

## Calculate Fst for all windows between all groups
# initialize table of window fst
fst.df <- data.frame(matrix(nrow = 0, ncol = 8))

for (i in 1:22) {
  chr.df <- df[which(df$chr == i), ]
  index <- seq(1+half_w, length(chr.df$ID)-half_w, window+1)
  temp_fst.df <- data.frame(matrix(nrow = length(index), ncol = 8))
  # order pos from small to large
  chr.df <- chr.df[order(chr.df$pos), ]
  # chrom
  temp_fst.df[ , 1] <- i 
  # start pos
  beg <- lapply(index, function(x) {
    unlist(chr.df[(x-half_w), 3])
  })
  temp_fst.df[ , 2] <- unlist(beg)
  # end pos
  end <- lapply(index, function(x) {
    unlist(chr.df[(x+half_w), 3])
  })
  temp_fst.df[ , 3] <- unlist(end)
  # average Fst group2 and ancient 
  fst_group2_anc <- lapply(index, function(x) {
      win <- (x - half_w):(x + half_w)
      calcAvgFst(chr.df[win, 'group2_fq'], chr.df[win, 'anc_fq'], chr.df[win, 'group2_n'], chr.df[win, 'anc_n'])
    })
  temp_fst.df[ , 4] <- unlist(fst_group2_anc)
  # average Fst group2 and Outgroup
  fst_group2_out <- lapply(index, function(x) {
    win <- (x - half_w):(x + half_w)
    calcAvgFst(chr.df[win, 'group2_fq'], chr.df[win, 'out_fq'], chr.df[win, 'group2_n'], chr.df[win, 'out_n'])
  })
  temp_fst.df[ , 8] <- unlist(fst_group2_out)
  # average Fst ancient and Outgroup
  fst_anc_out <- lapply(index, function(x) {
    win <- (x - half_w):(x + half_w)
    calcAvgFst(chr.df[win, 'anc_fq'], chr.df[win, 'out_fq'], chr.df[win, 'anc_n'], chr.df[win, 'out_n'])
  })
  temp_fst.df[ , 5] <- unlist(fst_anc_out)
  # average Fst group2 and European
  fst_group2_eur <- lapply(index, function(x) {
    win <- (x - half_w):(x + half_w)
    calcAvgFst(chr.df[win, 'group2_fq'], chr.df[win, 'eur_fq'], chr.df[win, 'group2_n'], chr.df[win, 'eur_n'])
  })
  temp_fst.df[ , 6] <- unlist(fst_group2_eur)
  # average Fst European and Outgroup
  fst_eur_out <- lapply(index, function(x) {
    win <- (x - half_w):(x + half_w)
    calcAvgFst(chr.df[win, 'eur_fq'], chr.df[win, 'out_fq'], chr.df[win, 'eur_n'], chr.df[win, 'out_n'])
  })
  temp_fst.df[ , 7] <- unlist(fst_eur_out)
  # bind to main table
  fst.df <- rbind(fst.df, temp_fst.df)
}                                                                                       
colnames(fst.df) <- c('chr','start','end','fst_group2_anc','fst_out_anc','fst_group2_eur','fst_out_eur','fst_group2_out')

# Calculate divergence time from Fst
fst.df$d_group2_anc <- -log(1 - fst.df$fst_group2_anc)
fst.df$d_group2_out <- -log(1 - fst.df$fst_group2_out)
fst.df$d_anc_out <- -log(1 - fst.df$fst_out_anc)
fst.df$d_group2_eur <- -log(1 - fst.df$fst_group2_eur)
fst.df$d_eur_out <- -log(1 - fst.df$fst_out_eur)

# Calculate PBS for group1 and group2
fst.df$pbs_anc <- (fst.df$d_group2_anc + fst.df$d_anc_out - fst.df$d_group2_out) / 2
fst.df$pbs_eur <- (fst.df$d_group2_eur + fst.df$d_eur_out - fst.df$d_group2_out) / 2

write.table(fst.df[,c('chr','start','end','fst_group2_anc','fst_out_anc','fst_group2_eur','fst_out_eur','fst_group2_out','pbs_anc','pbs_eur')], file=out_path, quote=F, row.names=F, col.names=T) 
