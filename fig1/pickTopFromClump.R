## pickTopFromClump.R
#` Takes PLINK clumping output and looks for clumps with index SNPs with p-value
#` equal to 0 and clumped SNPs with same p-value. Selects new index SNP based on
#` largest test statistic. 

args = commandArgs(trailingOnly = TRUE)

ref_path <- args[1]
clump_path <- args[2]
uk_path <- args[3]

ref.df <- read.table(ref_path, as.is=T, header=T)
clump.df <- read.table(clump_path, as.is=T, header=T)
uk.df <- read.table(uk_path, as.is=T, header=T)

# Format UK summary stats table
variant <- strsplit(uk.df[ , 1], ":")
uk.df$id <- lapply(variant, function(x){
  paste(x[1],x[2],sep='_')})
uk.df$id <- unlist(uk.df$id)

clump_in <- clump.df$SNP[clump.df$P>0]
clump_in <- c(clump_in, clump.df$SNP[clump.df$P==0 & clump.df$SP2=='NONE'])

# Subset p-val = 0 index SNPs
p0.df <- clump.df[clump.df$P==0 & clump.df$SP2!='NONE', ]

# Extract rsIDs
for (i in 1:nrow(p0.df)) {
  sp2 <- p0.df[i,12]
  sp2 <- as.list(strsplit(sp2, split=',')[[1]])
  sp2_rsid <- lapply(sp2, function(x) {
    strsplit(unlist(x[1]), split="\\(")[[1]][1]
    })
  sp2_rsid <- unlist(sp2_rsid)
  # Add index SNP
  sp2_rsid <- c(sp2_rsid, p0.df$SNP[i])
  # Get test stat and p-values of clump members
  sp2.df <- ref.df[ref.df$SNP %in% sp2_rsid, ]
  # Keep only p-val = 0
  sp2.df <- sp2.df[sp2.df$P==0, ]
  sp2.df$id <- paste(sp2.df$chr, sp2.df$pos, sep='_')
  sp2.df <- merge(sp2.df, uk.df[,c('id','tstat')], by='id', all.x=T)
  top <- which(abs(sp2.df$tstat)==max(abs(sp2.df$tstat)))
  clump_in <- c(clump_in, sp2.df$SNP[top])
  }
clump_in <- c('SNP',clump_in)

write.table(clump_in, file='final_clump.txt', row.names=F, col.names=F,
            quote=F)
