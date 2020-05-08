## Plot subfigure 4A
#` Note: need to run the UKBio time series pipeline beforehand and 
#` mask out large effect loci

library(ggplot2)

se_flag <- T

prepDf <- function(score.df, anno.df, pc.df) {
  score.df$Sample <- lapply(score.df$ID, function(x) {
    strsplit(x, "_")[[1]][1] })
  score.df$Sample <- unlist(score.df$Sample)
  
  score.anno.df <- merge(score.df, anno.df, by.x = 'ID', by.y = 'Instance.ID', 
                        all.x=T)
  fix_sample <- which(score.anno.df$ID != score.anno.df$Sample & 
                        score.anno.df$Date>0)
  score.anno.df$Sample[fix_sample] <- score.anno.df$ID[fix_sample]
  final.df <- merge(score.anno.df, pc.df, by.x='Sample', by.y='ID')
  
  colnames(final.df)[which(colnames(final.df)=='V12')] <- 'Population_pca'
  colnames(final.df)[which(colnames(final.df)=='Latitude')] <- 'Lat.'
  colnames(final.df)[which(colnames(final.df)=='Longitude')] <- 'Long.'
  
  # make time for 1KGP samples 0 
  final.df$Date[which(final.df$Population_pca != 'ANC')] <- 0
  # assign ANC population to ancient samples
  final.df$Population_pca[final.df$ID %in% anno.df$ID] <- 'ANC'
  final.df$Population_pca[final.df$ID %in% anno.df$Instance.ID] <- 'ANC'
  # delete data with missing values 
  final.df <- final.df[!is.na(final.df$Date), ]
  # subset data frames into ancient only
  final.df[which(final.df$Population_pca=='ANC'), ]
}

## paths to CSV of samples' PRS from getPRSHaploid.py 
score_path <- ''
score_path2 <- ''
score_path3 <- ''
score_path4 <- ''
score_path5 <- ''
score_path6 <- ''

# sample annotation file
anno_path <- 'v37.2_anno_FINAL.csv'
# evec file from smartPCA for the samples
pc_path <- ''
# folder path for output
out <- ''
# fig1 folder path
source_path <- ''

source(paste(source_path, 'plotPRStime.R', sep='/'))
source(paste(source_path, 'plotPRStimeGroup.R', sep='/'))
source(paste(source_path, 'add1khgLocation.R', sep='/'))

# load data 
score.df = read.csv(score_path, as.is = TRUE, header = TRUE)
score.df2 = read.csv(score_path2, as.is = TRUE, header = TRUE)
score.df3 = read.csv(score_path3, as.is = TRUE, header = TRUE)
score.df4 = read.csv(score_path4, as.is = TRUE, header = TRUE)
score.df5 = read.csv(score_path5, as.is = TRUE, header = TRUE)
score.df6 = read.csv(score_path6, as.is = TRUE, header = TRUE)

anno.df = read.csv(anno_path, as.is = TRUE, na.strings="..", header = TRUE)
pc.df = read.table(pc_path, as.is=T)

# merge data frames
colnames(pc.df)[1] <- 'ID'
 
# set up data frames for plotting
plot.df1 <- prepDf(score.df, anno.df, pc.df)
plot.df2 <- prepDf(score.df2, anno.df, pc.df)
plot.df3 <- prepDf(score.df3, anno.df, pc.df)
plot.df4 <- prepDf(score.df4, anno.df, pc.df)
plot.df5 <- prepDf(score.df5, anno.df, pc.df)
plot.df6 <- prepDf(score.df6, anno.df, pc.df)

# these are unresidualized plots
ggplot(plot.df1, aes(x=Date, y=Score)) + 
  theme_classic(base_size=15) + 
  theme(axis.title=element_text(size=18)) +
  scale_x_reverse() + 
  xlab("Date (years BP)") + ylab("Score") +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  stat_smooth(method="glm", aes(weight=Weight), 
              method.args=list(family="binomial"), se=se_flag, fill='blue')  +
  stat_smooth(data=plot.df2, method="glm", aes(weight=Weight), 
              method.args=list(family="binomial"), se=se_flag,
              color='red', fill='red')  +
  stat_smooth(data=plot.df3, method="glm", aes(weight=Weight), 
              method.args=list(family="binomial"), se=se_flag,
              color='magenta', fill='magenta')  +
  stat_smooth(data=plot.df4, method="glm", aes(weight=Weight), 
              method.args=list(family="binomial"), se=se_flag,
              color='green', fill='green')  +
  stat_smooth(data=plot.df5, method="glm", aes(weight=Weight), 
              method.args=list(family="binomial"), se=se_flag,
              color='dimgray', fill='dimgray')  +
  stat_smooth(data=plot.df6, method="glm", aes(weight=Weight), 
              method.args=list(family="binomial"), se=se_flag,
              color='black')  

# plot residualized scores based on time 
runModel <- function(df) {
  model_pc <- glm(Score ~ Date + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 
                  + V10 + V11 + Lat. + Long., weights = Weight,
                  family = "binomial", 
                  data = df)
}

lm1 <- runModel(plot.df1)
lm2 <- runModel(plot.df2)
lm3 <- runModel(plot.df3)
lm4 <- runModel(plot.df4)
lm5 <- runModel(plot.df5)
lm6 <- runModel(plot.df6)

getPredictedScores <- function(lm) {
  intercept <- coefficients(lm)[1]
  beta <- coefficients(lm)[2]
  times <- data.frame(seq(0,4e4,100))
  # calculate scores based on beta
  log_odds <- times * beta + intercept
  odds <- exp(log_odds)
  Score <- odds / (1 + odds)
  pred.lm <- data.frame(cbind(times,Score))
  colnames(pred.lm) <- c('Date','Score')
  pred.lm
  }

lm1.df <- getPredictedScores(lm1)
lm1.df$group <- 1
lm2.df <- getPredictedScores(lm2)
lm2.df$group <- 2
lm3.df <- getPredictedScores(lm3)
lm3.df$group <- 3
lm4.df <- getPredictedScores(lm4)
lm4.df$group <- 4
lm5.df <- getPredictedScores(lm5)
lm5.df$group <- 5
lm6.df <- getPredictedScores(lm6)
lm6.df$group <- 6

plot.df <- rbind(lm1.df, lm2.df, lm3.df, lm4.df, lm5.df, lm6.df)

width=0.8
ggplot(data = plot.df, aes(x = Date, y = Score)) +
  theme_classic(base_size=15) +
  theme(legend.background=element_rect(color='white', size=0.5),
        legend.position=c(0.8,0.8),
        legend.key.size=unit(0.3,'cm'),
        legend.title=element_text(size=10)
        ) + 
  scale_x_reverse() + 
  scale_y_continuous(limits=c(0.7,1), breaks=c(0.7,0.85,1)) +
  stat_smooth(aes(x = Date, y = Score, color=as.factor(group)), se=F) +
  scale_color_manual(name="# SNPs removed", 
                     labels=c('0','1','2','3','4','5','6'),
                     values=c('blue','green','orange','cyan','darkgrey','red')
                     )
  
ggsave(file = paste(out,'fig4a.pdf',sep='/'), width = 4, height = 4)
