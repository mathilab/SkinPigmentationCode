## Make plots of time series of PRS controlling for ancestry 
#` Update: weighted score regression now using logistic regression

args = commandArgs(trailingOnly = TRUE)
score_path <- args[1]
anno_path <- args[2]
pc_path <- args[3]
use_beta <- args[4]
avg_score_path <- args[5]
use_anc <- args[6]
out <- args[7]
source_path <- args[8]
loc_path <- args[9]

source(paste(source_path, 'plotPRStime.R', sep='/'))
source(paste(source_path, 'plotPRStimeGroup.R', sep='/'))
source(paste(source_path, 'add1khgLocation.R', sep='/'))
       
# load data 
score.df = read.csv(score_path, as.is = TRUE, header = TRUE)
anno.df = read.csv(anno_path, as.is = TRUE, na.strings="..", header = TRUE)
pc.df = read.table(pc_path, as.is=T)
pop.df = read.table(avg_score_path, as.is=T, header=T)
loc.df = read.csv(loc_path, as.is=T, header=T)

# remove AMR from super population averages
pop.df <- pop.df[which(pop.df$Population!='AMR'), ]
pop.df$Date <- 0

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

# perform classification 
if(use_beta == T | use_beta == "T") {
  # only ancient samples
  model_pc <- glm(Score ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 
                 + V10 + V11 + Lat. + Long., weights = Weight,
                 family = "binomial", 
                 data = anc.df)
  # ancient samples + 1kgp Europeans
  model_pc_eur <- glm(Score ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 
                       + V10 + V11 + Lat. + Long., weights = Weight,
                       family = "binomial", 
                       data = eur.df)
  # save model summary 
  sink(paste(out,'_anc.txt',sep=''))
  print(summary(model_pc))
  sink()
  
  sink(paste(out,'_anc+1kgpEur.txt',sep=''))
  print(summary(model_pc_eur))
  sink()
  # make plot  
  plotPRStime(use_beta,anc.df,eur.df,pop.df,out)
  
  if(use_anc == T | use_anc == "true"){
    # using admixture ancestry components 
    model_ancestry <- glm(Score ~ Date + ANA + YAM + Lat. + Long.,
                    family = "binomial", weights=Weight,
                    data = anc.df)
    # save model summary
    sink(paste(out,'_ancestry.txt',sep=''))
    print(summary(model_ancestry))
    sink()
    # make plot and model summaries for each ancestral group
    plotPRStimeGroup(anc.df,pop.df,use_beta,out)
  }
} else { 
  # run GLM for ancient
  anc.df$NumLight <- anc.df$TotalNumSNPs - anc.df$NumDark
  values_anc <- cbind(anc.df$NumDark, anc.df$NumLight)
  colnames(values_anc) <- c('dark', 'light')
  
  model_pc <- glm(values_anc ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + 
                  V8 + V9 + V10 + V11 + Lat. + Long.,
                  family = "binomial", data = anc.df)

  # run GLM for ancient + 1KGP Europeans
  eur.df$NumLight <- eur.df$TotalNumSNPs - eur.df$NumDark
  values_eur <- cbind(eur.df$NumDark, eur.df$NumLight)
  colnames(values_eur) <- c('dark', 'light')
  
  model_pc_eur <- glm(values_eur ~ Date + V2 + V3 + V4 + V5 + V6 + 
                      V7 + V8 + V9 + V10 + V11 + Lat. + Long.,
                      family = "binomial", data = eur.df)
  
  # save model summary 
  sink(paste(out,'_anc.txt',sep=''))
  print(summary(model_pc))
  sink()
  
  sink(paste(out,'_anc+1kgpEur.txt',sep=''))
  print(summary(model_pc_eur))
  sink()
  
  # make plot 
  plotPRStime(use_beta,anc.df,eur.df,pop.df,out)
  
  if(use_anc == T | use_anc == "true"){
    # using admixture ancestry components 
    model_ancestry <- glm(values_anc ~ Date + ANA + YAM + Lat. + Long.,
                     family = "binomial", 
                     data = anc.df)
    # save model summary
    sink(paste(out,'_admix.txt',sep=''))
    print(summary(model_ancestry))
    sink()
    # make plot and model summaries for each ancestral group
    plotPRStimeGroup(anc.df,pop.df,use_beta,out)
  }
}
