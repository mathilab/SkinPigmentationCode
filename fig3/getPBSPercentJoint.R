## getPBSPercent.R
## Get PBS for a window centered around specific trait SNP

args = commandArgs(trailingOnly = TRUE)
df_path <- args[1]
window <- as.integer(args[2])
trait_snps <- args[3]
pbs_path <- args[4]
out_path <- args[5]

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

# Read data
pbs.df <- fread(pbs_path, stringsAsFactors = F, header=T)
df <- fread(df_path, stringsAsFactors = F)
snps <- scan(trait_snps, what = character())
half_w <- window / 2

# Extract chr and pos from trait SNPs
snps.df <- data.frame(snps)
snps.df$chr <- lapply(snps, function(x) {
  strsplit(x, split="_")[[1]][[1]]
} )
snps.df$pos <- lapply(snps, function(x) {
  strsplit(x, split="_")[[1]][[2]]
} )

## Calculate Fst for all windows between all groups
# initialize table of window fst
trait.df <- data.frame(matrix(nrow = length(snps), ncol = 9))

for (i in 1:length(snps)) {
  chr.df <- df[df$chr == snps.df$chr[i]]
  # order pos from small to large
  chr.df <- chr.df[order(chr.df$pos), ]
  try({
  index <- which(chr.df$ID == snps[i])
  # print SNP if no index 
  if (length(index)==0) {
    print(snps[i])
  }
  # trait SNP ID
  trait.df[i, 9] <- snps[i]
  # chrom
  trait.df[i, 1] <- snps.df$chr[i]
  # start pos
  trait.df[i, 2] <- chr.df[(index-half_w), 3]
  # end pos
  trait.df[i, 3] <- chr.df[(index+half_w), 3]
  # average Fst group2 and group1
  win <- (index - half_w):(index + half_w)
  # remove negative values in window
  if (win[1] < 0) {
  print(snps[i])
  print('Window around SNP at chromosome edge!')
  }
  win <- win[win > 0]
  # average Fst group2 and ancient
  trait.df[i, 4] <- calcAvgFst(chr.df[win, 'group2_fq'], chr.df[win, 'anc_fq'], chr.df[win, 'group2_n'], chr.df[win, 'anc_n'])
  # average Fst group2 and Outgroup
  trait.df[i, 8] <- calcAvgFst(chr.df[win, 'group2_fq'], chr.df[win, 'out_fq'], chr.df[win, 'group2_n'], chr.df[win, 'out_n'])
  # average Fst ancient and Outgroup
  trait.df[i, 5] <- calcAvgFst(chr.df[win, 'anc_fq'], chr.df[win, 'out_fq'], chr.df[win, 'anc_n'], chr.df[win, 'out_n'])
  # average Fst group2 and European
  trait.df[i, 6] <- calcAvgFst(chr.df[win, 'group2_fq'], chr.df[win, 'eur_fq'], chr.df[win, 'group2_n'], chr.df[win, 'eur_n'])
  # average Fst ancient and European
  trait.df[i, 7] <- calcAvgFst(chr.df[win, 'eur_fq'], chr.df[win, 'out_fq'], chr.df[win, 'eur_n'], chr.df[win, 'out_n'])
  })
}

colnames(trait.df) <- c('chr','start','end','fst_group2_anc','fst_anc_out','fst_group2_eur','fst_eur_out','fst_group2_out','ID')

trait.df$d_group2_anc <- -log(1 - trait.df$fst_group2_anc)
trait.df$d_group2_out <- -log(1 - trait.df$fst_group2_out)
trait.df$d_anc_out <- -log(1 - trait.df$fst_anc_out)
trait.df$d_group2_eur <- -log(1 - trait.df$fst_group2_eur)
trait.df$d_eur_out <- -log(1 - trait.df$fst_eur_out)

# Calculate PBS for group1 and group2
trait.df$pbs_anc <- (trait.df$d_group2_anc + trait.df$d_anc_out - trait.df$d_group2_out) / 2
trait.df$pbs_eur <- (trait.df$d_group2_eur + trait.df$d_eur_out - trait.df$d_group2_out) / 2

# Calculate percentiles
perc_anc <- ecdf(pbs.df$pbs_anc)
perc_eur <- ecdf(pbs.df$pbs_eur)

trait.df$anc_perc <- 1 - perc_anc(trait.df$pbs_anc)
trait.df$eur_perc <- 1 - perc_eur(trait.df$pbs_eur)

write.table(trait.df[,c('chr','start','end','fst_group2_anc','fst_anc_out','fst_group2_eur','fst_eur_out','fst_group2_out','pbs_anc','pbs_eur','anc_perc','eur_perc','ID')], file=out_path, row.names=F, quote=F)
