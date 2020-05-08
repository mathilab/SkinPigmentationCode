## Run logistic regression on alternate frequency matrix for all SNPs 
#' Outputs beta coefficient and p-value for each SNP
#' Models ancestry with first 10 PCs
#` Update: 2.11.2020 Calculates partial r2 for ancestry and time terms

library(rsq)

args = commandArgs(trailingOnly = T) 
alt_path <- args[1] # matrix of alternate allele
anno_path <- args[2] # annotations with date and ancestry
export_path <- args[3]
bad_path <- args[4] # path to where warning/error SNPs go
pc_path <- args[5] # corresponding pca.evec file

# load data
alt.df = read.csv(alt_path, as.is = T, header = T, na.strings = '.', check.names=F)
anno.df = read.csv(anno_path, as.is = T, header = T)
pca.df <- read.table(pc_path, as.is=T)

# merge annotation, assignment and alternate allele datasets
snps <- alt.df[ , 1]
alt.df <- alt.df[ , 2:ncol(alt.df)]
alt.df <- data.frame(t(alt.df))
colnames(alt.df) <- snps
alt.df$Instance.ID <- rownames(alt.df)
# choose column for merging based on match proportion
match <- length(which(anno.df$Instance.ID %in% alt.df$Instance.ID))
if (match < nrow(anno.df)-1) {
  group.df <- merge(anno.df, alt.df, by.x = "ID", by.y = "Instance.ID")
} else {
  group.df <- merge(anno.df, alt.df, by = "Instance.ID")
}
group.df <- merge(group.df, pca.df, by.x = "Instance.ID", by.y = "V1")

# make data frame for output
final.df <- data.frame(matrix(ncol = 8, nrow = length(snps)))
colnames(final.df) <- c('SNP', 'beta_t', 'beta_t_anc', 'p_t', 'p_anc', 'n',
                        't_r2','anc_r2')
final.df$SNP <- snps

# run GLM on all SNPs
bad_snps = {}
for (x in snps) {
  bmodel_anc <- NA
  bmodel_time <- NA
  bmodel <- NA
  
  # set up input for glm
  alt <- group.df[ , which(colnames(group.df)==x)]
  ref <- abs(alt - 1)
  values <- cbind(alt, ref)
  colnames(values) <- c('alt', 'ref')
  
  # run glm alternate ~ time + 10 PCs
  tryCatch(
  bmodel <- glm(values ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, 
                    family = "binomial", data = group.df),
    warning = function(w) {
      bad_snps <<- append(bad_snps, x)
      },
    error = function(err){
      bad_snps <<- append(bad_snps, x)
    })
  
  # run glm alternate ~ time 
  tryCatch(
    bmodel_time <- glm(values ~ Date, family = "binomial", data = group.df),
    warning = function(w) {
      bad_snps <<- append(bad_snps, x)
    },
    error = function(err){
      bad_snps <<- append(bad_snps, x)
    })

  # run glm alternate ~ 10 PCs 
  tryCatch(
    bmodel_anc <- glm(values ~ V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, family = "binomial", data = group.df),
    warning = function(w) {
      bad_snps <<- append(bad_snps, x)
    },
    error = function(err){
      bad_snps <<- append(bad_snps, x)
    })
  
  if (is.na(bmodel)[1]==T | is.na(bmodel_anc)[1]==T | is.na(bmodel_time)[1]==T) {
    next
  }
  
  # get p-values
  test_anc <- anova(bmodel_time, bmodel, test='Chisq')
  test_time <- anova(bmodel_anc, bmodel, test='Chisq')
  # fill in table
  final.df$beta_t[which(final.df$SNP==x)] <- bmodel_time$coefficients[[2]]
  final.df$beta_t_anc[which(final.df$SNP==x)] <- bmodel$coefficients[[2]]
  final.df$p_t[which(final.df$SNP==x)] <- test_time[[5]][2]
  final.df$p_anc[which(final.df$SNP==x)] <- test_anc[[5]][2]
  final.df$n[which(final.df$SNP==x)] <- nobs(bmodel)
  
  # get partial r2
  anc_r2 <- rsq.partial(bmodel, bmodel_time, adj=F)
  time_r2 <- rsq.partial(bmodel, bmodel_anc, adj=F)
  
  # fill in table
  final.df$t_r2[which(final.df$SNP==x)] <- time_r2$partial.rsq
  final.df$anc_r2[which(final.df$SNP==x)] <- anc_r2$partial.rsq
  
  rm(bmodel)
  rm(bmodel_anc)
  rm(bmodel_time)
}

# remove SNPs with warnings
final.df <- final.df[which(!(final.df$SNP %in% bad_snps)), ]

# save SNPs with warnings to separate file
write.table(bad_snps, bad_path, sep="\t", row.names=F, quote=F )

# export table with model coefficients
write.table(final.df, export_path, sep="\t", row.names=F, quote=F)
