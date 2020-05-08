## runPRSControl.R
#` Get summary statistics for randomized control sets PRS
#` Note: should only be run in pipeline

extract <- function(m) {
  # extracts date summary stats from linear model
  sm <- summary(m)
  m.df <- data.frame(sm$coefficients)
  date_stats <- m.df[2, ]
  return(date_stats)
}

args = commandArgs(trailingOnly = TRUE)
score_path <- args[1]
anno_path <- args[2]
pc_path <- args[3]
use_beta <- args[4] # use effect size estimates in PRS calculation
loc_path <- args[5]
out_dir <- args[6] # the out folder where 4 subdirectories located
set <- args[7] # used for file naming
source_path <- args[8]
use_anc <- args[9]

source(paste(source_path, 'add1khgLocation.R', sep='/'))

# load data
score.df = read.csv(score_path, as.is = TRUE, header = TRUE)
anno.df = read.csv(anno_path, as.is = TRUE, na.strings="..", header = TRUE)
pc.df = read.table(pc_path, as.is=T)
loc.df = read.csv(loc_path, as.is=T, header=T)

score.df$Sample <- lapply(score.df$ID, function(x) {
  strsplit(x, "_")[[1]][1] })
score.df$Sample <- unlist(score.df$Sample)

# merge data frames
colnames(pc.df)[1] <- 'ID'

# choose column for merging based on match proportion
match <- length(which(anno.df$Instance.ID %in% score.df$ID))
if (match < (nrow(anno.df)-nrow(anno.df)/2)) { #arbitrary for poor matching
  score.anno.df <- merge(score.df, anno.df, by = 'ID', all.x=T)
  fix_sample <- which(score.anno.df$Instance.ID != score.anno.df$Sample & 
                        score.anno.df$Date>0)
  score.anno.df$Sample[fix_sample] <- score.anno.df$Instance.ID[fix_sample]
  final.df <- merge(score.anno.df, pc.df, by.y = 'ID', by.x = 'Sample')
} else {
  score.anno.df <- merge(score.df, anno.df, by.x = 'ID', by.y = 'Instance.ID', 
                         all.x=T)
  fix_sample <- which(score.anno.df$ID != score.anno.df$Sample & 
                        score.anno.df$Date>0)
  score.anno.df$Sample[fix_sample] <- score.anno.df$ID[fix_sample]
  final.df <- merge(score.anno.df, pc.df, by.x='Sample', by.y='ID')
}

colnames(final.df)[which(colnames(final.df)=='V12')] <- 'Population_pca'
colnames(final.df)[which(colnames(final.df)=='Latitude')] <- 'Lat.'
colnames(final.df)[which(colnames(final.df)=='Longitude')] <- 'Long.'

# make time for 1KGP samples 0 
final.df$Date[which(final.df$Population_pca != 'ANC')] <- 0
# add longitude and latitude to 1KHG samples 
final.df <- add1khgLocation(final.df, loc.df)
# assign ANC population to ancient samples
final.df$Population_pca[final.df$ID %in% anno.df$ID] <- 'ANC'
final.df$Population_pca[final.df$ID %in% anno.df$Instance.ID] <- 'ANC'

# delete data with missing values 
final.df <- final.df[!is.na(final.df$Date), ]
# subset data frames into ancient only and ancient + 1KGP Europeans
anc.df <- final.df[which(final.df$Population_pca=='ANC'), ]
eur.df <- final.df[which(final.df$Population_pca %in% 
                           c('ANC','CEU','TSI','FIN','GBR','IBS')), ]

if(use_beta == T | use_beta == "true") {
  ## Run GLM for only ancient samples
  n_anc <- nrow(anc.df)
  # plain regression
  anc_model <- glm(Score ~ Date, family = "binomial", weights=Weight, data=anc.df)
  ap <- extract(anc_model)
  anc_plain_row <- c(ap[1], ap[2], ap[4], n_anc)
  write.table(anc_plain_row, file=paste(out_dir, 'anc/plain', set, sep='/'),
              col.names=F, row.names=F, quote=F)
  
  # using PCs
  anc_model_pc <- glm(Score ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
                 + V10 + V11 + Lat. + Long., weights=Weight,
                 family = "binomial",
                 data = anc.df)
  apc <- extract(anc_model_pc)
  anc_pc_row <- c(apc[1], apc[2], apc[4], n_anc)
  write.table(anc_pc_row, file=paste(out_dir, 'anc/pc', set, sep='/'),
              col.names=F, row.names=F, quote=F)
  
  if(use_anc == T | use_anc == "true"){
    # Partition dataset by ADMIXTURE proportions
    # Early Farmer
    ef.df <- anc.df[anc.df$ANA>0.6, ]
    filter_ef <- ef.df$ID[which(ef.df$Date<5000 & ef.df$YAM>0.3)]
    ef.df <- ef.df[!(ef.df$ID %in% filter_ef), ]
    # Hunter gatherer
    hg.df <- anc.df[which(anc.df$HG>0.6), ]
    # Steppe/Yamnaya
    filter_sp <- which(anc.df$Date<=5000 & anc.df$YAM>0.3)
    sp.df <- anc.df[filter_sp, ]
    
    # farmer model
    n_ef <- nrow(ef.df)
    model.ef <- glm(Score ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 
                    + V10 + V11 + Lat. + Long., 
                    family="binomial", weights = Weight, data=ef.df)
    ef <- extract(model.ef)
    ef_row <- c(ef[1], ef[2], ef[4], n_ef)
    write.table(ef_row, file=paste(out_dir, 'anc/ef', set, sep='/'),
                col.names=F, row.names=F, quote=F)
    # steppe 
    n_sp <- nrow(sp.df)
    model.sp <- glm(Score ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 
                    + V10 + V11 + Lat. + Long., 
                    family="binomial", weights = Weight, data=sp.df)
    sp <- extract(model.sp)
    sp_row <- c(sp[1], sp[2], sp[4], n_sp)
    write.table(sp_row, file=paste(out_dir, 'anc/yam', set, sep='/'),
                col.names=F, row.names=F, quote=F)
    # hunter gatherer
    n_hg <- nrow(hg.df)
    model.hg <- glm(Score ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 
                    + V10 + V11 + Lat. + Long.,
                    family="binomial", weights = Weight, data=hg.df)
    hg <- extract(model.hg)
    hg_row <- c(hg[1], hg[2], hg[4], n_hg)
    write.table(hg_row, file=paste(out_dir, 'anc/hg', set, sep='/'),
                col.names=F, row.names=F, quote=F)
  }
  # using admixture ancestry proportions
  # anc_model_ad <- glm(Score ~ Date + ANA + YAM + Latitude + Longitude,
  #                       family = "gaussian",
  #                       data = anc.df)
  # aad <- extract(anc_model_ad)
  # anc_ad_row <- c(aad[1], aad[2], aad[4], n_anc)
  # write.table(anc_ad_row, file=paste(out_dir, 'anc/admix', set, sep='/'),
  #             col.names=F, row.names=F, quote=F)
  
  ## Run GLM for 1000 Genomes + ancient samples
  # n_eur <- nrow(eur.df)
  # plain regression
  # eur_model <- glm(Score ~ Date, family = "gaussian", data = eur.df)
  # ep <- extract(eur_model)
  # eur_plain_row <- c(ep[1], ep[2], ep[4], n_eur)
  # write.table(eur_plain_row, file=paste(out_dir, 'anc-1kg/plain', set, sep='/'),
  #             col.names=F, row.names=F, quote=F)
  # use PCs 
  # model_pc_eur <- glm(Score ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9
  #                     + V10 + V11 + Latitude + Longitude,
  #                     family = "gaussian",
  #                     data = eur.df)
  # epc <- extract(model_pc_eur)
  # eur_pc_row <- c(epc[1], epc[2], epc[4], n_eur)
  # write.table(eur_pc_row, file=paste(out_dir, 'anc-1kg/pc', set, sep='/'),
  #             col.names=F, row.names=F, quote=F)
} else {
  ## Run GLM for just ancient
  anc.df$NumLight <- anc.df$TotalNumSNPs - anc.df$NumDark
  values_anc <- cbind(anc.df$NumDark, anc.df$NumLight)
  colnames(values_anc) <- c('dark', 'light')
  n_anc <- nrow(anc.df)

  # plain regression with no covariates
  anc_model <- glm(values_anc ~ Date, 
                   family = "binomial", data = anc.df)
  ap <- extract(anc_model)
  anc_plain_row <- c(ap[1], ap[2], ap[4], n_anc)
  write.table(anc_plain_row, file=paste(out_dir, 'anc/plain', set, sep='/'),
              col.names=F, row.names=F, quote=F)

  # using PCs
  anc_model_pc <- glm(values_anc ~ Date + V2 + V3 + V4 + V5 + V6 + V7 +
                  V8 + V9 + V10 + V11 + Lat. + Long.,
                  family = "binomial", data = anc.df)
  apc <- extract(anc_model_pc)
  anc_pc_row <- c(apc[1], apc[2], apc[4], n_anc)
  write.table(anc_pc_row, file=paste(out_dir, 'anc/pc', set, sep='/'),
              col.names=F, row.names=F, quote=F)

  if(use_anc == T | use_anc == "true"){
    # Partition dataset by ADMIXTURE proportions
    # Early Farmer
    ef.df <- anc.df[anc.df$ANA>0.6, ]
    filter_ef <- ef.df$ID[which(ef.df$Date<5000 & ef.df$YAM>0.3)]
    ef.df <- ef.df[!(ef.df$ID %in% filter_ef), ]
    # Hunter gatherer
    hg.df <- anc.df[which(anc.df$HG>0.6), ]
    # Steppe/Yamnaya
    filter_sp <- which(anc.df$Date<=5000 & anc.df$YAM>0.3)
    sp.df <- anc.df[filter_sp, ]
    
    # farmer model
    n_ef <- nrow(ef.df)
    values.ef <- cbind(ef.df$NumDark, ef.df$NumLight)
    colnames(values.ef) <- c('dark', 'light')
    model.ef <- glm(values.ef ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 
                    + V10 + V11 + Lat. + Long., 
                    family = "binomial", data=ef.df)
    
    ef <- extract(model.ef)
    ef_row <- c(ef[1], ef[2], ef[4], n_ef)
    write.table(ef_row, file=paste(out_dir, 'anc/ef', set, sep='/'),
                col.names=F, row.names=F, quote=F)
    # steppe 
    n_sp <- nrow(sp.df)
    values.sp <- cbind(sp.df$NumDark, sp.df$NumLight)
    colnames(values.sp) <- c('dark', 'light')
    
    model.sp <- glm(values.sp ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 
                    + V10 + V11 + Lat. + Long., 
                    family="binomial", data=sp.df)
    sp <- extract(model.sp)
    sp_row <- c(sp[1], sp[2], sp[4], n_sp)
    write.table(sp_row, file=paste(out_dir, 'anc/yam', set, sep='/'),
                col.names=F, row.names=F, quote=F)
    # hunter gatherer
    n_hg <- nrow(hg.df)
    values.hg <- cbind(hg.df$NumDark, hg.df$NumLight)
    colnames(values.hg) <- c('dark', 'light')
    
    model.hg <- glm(values.hg ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 
                    + V10 + V11 + Lat. + Long., 
                    family="binomial",data=hg.df)
    hg <- extract(model.hg)
    hg_row <- c(hg[1], hg[2], hg[4], n_hg)
    write.table(hg_row, file=paste(out_dir, 'anc/hg', set, sep='/'),
                col.names=F, row.names=F, quote=F)
  }
  # using admix proportions
  # anc_model_ad <- glm(values_anc ~ Date + ANA + YAM + Latitude + Longitude,
  #                 family = "binomial", data = anc.df)
  # aad <- extract(anc_model_ad)
  # anc_ad_row <- c(aad[1], aad[2], aad[4], n_anc)
  # write.table(anc_ad_row, file=paste(out_dir, 'anc/admix', set, sep='/'),
  #             col.names=F, row.names=F, quote=F)

  ## Run GLM for ancient + 1KGP Europeans
  # eur.df$NumLight <- eur.df$TotalNumSNPs - eur.df$NumDark
  # values_eur <- cbind(eur.df$NumDark, eur.df$NumLight)
  # colnames(values_eur) <- c('dark', 'light')
  # n_eur <- nrow(eur.df)

  # plain regression with no covariates
  # eur_model <- glm(values_eur ~ Date, family = "binomial", data = eur.df)
  # ep <- extract(eur_model)
  # eur_plain_row <- c(ep[1], ep[2], ep[4], n_eur)
  # write.table(eur_plain_row, file=paste(out_dir, 'anc-1kg/plain', set, sep='/'),
  #             col.names=F, row.names=F, quote=F)

  # using PCs
  # model_pc_eur <- glm(values_eur ~ Date + V2 + V3 + V4 + V5 + V6 +
  #                     V7 + V8 + V9 + V10 + V11 + Latitude + Longitude,
  #                     family = "binomial", data = eur.df)
  # epc <- extract(model_pc_eur)
  # eur_pc_row <- c(epc[1], epc[2], epc[4], n_eur)
  # write.table(eur_pc_row, file=paste(out_dir, 'anc-1kg/pc', set, sep='/'),
  #             col.names=F, row.names=F, quote=F)

}
