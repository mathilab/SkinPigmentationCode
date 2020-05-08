## Plot PRS time series stratified by European ancestral group

plotPRStimeGroup <- function(data.df,avg.df,use_beta,out) {
  # Partition dataset by ADMIXTURE proportions
  # Early Farmer
  ef.df <- data.df[data.df$ANA>0.6, ]
  filter_ef <- ef.df$ID[which(ef.df$Date<5000 & ef.df$YAM>0.3)]
  ef.df <- ef.df[!(ef.df$ID %in% filter_ef), ]
  # Hunter gatherer
  hg.df <- data.df[which(data.df$HG>0.6), ]
  # Steppe/Yamnaya
  filter_sp <- which(data.df$Date<=5000 & data.df$YAM>0.3)
  sp.df <- data.df[filter_sp, ]
  # assign factors for legend 
  data.df$color[which(data.df$HG>0.6)] <- 'HG'
  data.df$color[data.df$ID %in% ef.df$ID] <- 'EF'
  data.df$color[filter_sp] <- 'SP'
  
  if(use_beta==T | use_beta == "T") {
    # run GLM for each group
    # early farmer
    model.ef <- glm(Score ~ Date, family="binomial", weights = Weight, data=ef.df)
    sink(paste(out,'_ef.txt',sep=''))
    print(summary(model.ef))
    sink()
    # steppe 
    model.sp <- glm(Score ~ Date, family="binomial", weights = Weight, data=sp.df)
    sink(paste(out,'_sp.txt',sep=''))
    print(summary(model.sp))
    sink()
    # hunter gatherer
    model.hg <- glm(Score ~ Date, family="binomial", weights = Weight, data=hg.df)
    sink(paste(out,'_hg.txt',sep=''))
    print(summary(model.hg))
    sink()
    # plot groups into one figure
    ggplot(data.df, aes(x=Date, y=Score)) + theme_classic() +
      scale_color_manual(name="Population", 
                         labels=c('EF','HG','SP'),
                         values=c('steelblue1','red','palegreen2')
      ) + 
      theme(legend.position=c(0.13,0.85),
            legend.title=element_text(size=6),
            legend.text=element_text(size=6),
            legend.background=element_rect(size=0.15, linetype="solid", color='grey'),
            legend.key.size = unit(0.4, "cm")
      ) + guides(size=FALSE) +
      scale_x_reverse() + 
      xlab("Date (years BP)") + ylab("Score") +
      scale_y_continuous(breaks=c(0,0.5,1)) +
      geom_point(data=data.df[!(is.na(data.df$color)), ], 
                 mapping=aes(x=Date, y=Score, color=color, size=Weight^7), 
                 pch=16, alpha=0.6) +
      stat_smooth(method="glm", method.args=list(family="binomial"), se=F,
                  aes(weight=Weight), lty='longdash',
                  color=rgb(0,0,0,0.8)) +
      stat_smooth(data=ef.df,mapping=aes(x=Date, y=Score, weight=Weight), 
                  method="glm", method.args=list(family="binomial"), 
                  se=T, inherit.aes=F, color="blue", fill='blue') +
      stat_smooth(data=hg.df,mapping=aes(x=Date, y=Score, weight=Weight), 
                  method="glm", method.args=list(family="binomial"), 
                  se=T, inherit.aes=F, color='red3', fill='red3') +
      stat_smooth(data=sp.df,mapping=aes(x=Date, y=Score, weight=Weight), 
                  method="glm", method.args=list(family="binomial"), 
                  se=T, inherit.aes=F, color='darkgreen', fill='darkgreen') 
    
    ggsave(file = paste(out,'_stratified.pdf',sep=''), width = 4, height = 4) 
  } else {
    # run GLM for each group
    # early farmer 
    values.ef <- cbind(ef.df$NumDark, ef.df$NumLight)
    colnames(values.ef) <- c('dark', 'light')
    model.ef <- glm(values.ef ~ ef.df$Date, family = "binomial")
    sink(paste(out,'_ef.txt',sep=''))
    print(summary(model.ef))
    sink()
    # steppe 
    values.sp <- cbind(sp.df$NumDark, sp.df$NumLight)
    colnames(values.sp) <- c('dark', 'light')
    model.sp <- glm(values.sp ~ sp.df$Date, family = "binomial")
    sink(paste(out,'_sp.txt',sep=''))
    print(summary(model.sp))
    sink()
    # hunter gatherer
    values.hg <- cbind(hg.df$NumDark, hg.df$NumLight)
    colnames(values.hg) <- c('dark', 'light')
    model.hg <- glm(values.hg ~ hg.df$Date, family = "binomial")
    sink(paste(out,'_hg.txt',sep=''))
    print(summary(model.hg))
    sink()

    # plot groups into one figure
    ggplot(data.df, aes(x=Date, y=Score)) + theme_classic() +
      scale_color_manual(name="Population", 
                         labels=c('EF','HG','SP'),
                         values=c('steelblue1','red','palegreen2')
                         ) + 
      theme(legend.position=c(0.13,0.85),
            legend.title=element_text(size=6),
            legend.text=element_text(size=6),
            legend.background=element_rect(size=0.15, linetype="solid", color='grey'),
            legend.key.size = unit(0.4, "cm")
            ) +
      guides(size=F) +          
      scale_x_reverse() + 
      scale_y_continuous(breaks=c(0,0.5,1)) + 
      xlab("Date (years BP)") + ylab("Score") +
      geom_point(data=data.df[!(is.na(data.df$color)), ], 
                 mapping=aes(x=Date, y=Score, size=(TotalNumSNPs)^7, color=color), 
                 pch=16, alpha=0.6) +
      stat_smooth(method="glm", method.args=list(family="binomial"), se=F,
                  aes(weight=TotalNumSNPs), lty='longdash',
                  color=rgb(0,0,0,0.8)) +
      stat_smooth(data=ef.df,mapping=aes(x=Date, y=Score, weight=TotalNumSNPs), 
                  method="glm", method.args=list(family="binomial"), 
                  se=T, inherit.aes=F, color="blue", fill='blue') +
      stat_smooth(data=hg.df,mapping=aes(x=Date, y=Score, weight=TotalNumSNPs), 
                  method="glm", method.args=list(family="binomial"), 
                  se=T, inherit.aes=F, color='red3', fill='red3') +
      stat_smooth(data=sp.df,mapping=aes(x=Date, y=Score, weight=TotalNumSNPs), 
                  method="glm", method.args=list(family="binomial"), 
                  se=T, inherit.aes=F, color='darkgreen', fill='darkgreen')

    ggsave(file = paste(out,'_stratified.pdf',sep=''), width = 4, height = 4) 
  }
}
