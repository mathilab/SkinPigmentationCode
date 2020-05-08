## Plot time series PRS of ancient and present-day samples 
library(ggplot2)

plotPRStime <- function(use_beta,anc.df,eur.df,avg.df,out) {
  lwidth = 1
  a = 0.4
  
  ## Weighted linear regression
  if(use_beta==T | use_beta == "T") {
  # run simple regression
  model_anc <- glm(Score ~ Date, family="binomial", weights = Weight, data=anc.df)
  ord_anc <- order(anc.df$Date)    
  
  model_eur <- glm(Score ~ Date, family="binomial", weights = Weight, data=eur.df)
  ord_eur <- order(eur.df$Date)
  
  # save model summary
  sink(paste(out,'_anc_simple.txt',sep=''))
  print(summary(model_anc))
  sink()
  
  sink(paste(out,'_anc+1kgpEur_simple.txt',sep=''))
  print(summary(model_eur))
  sink()
  
  # plot only ancient humans
  ggplot(anc.df, aes(x=Date, y=Score)) + 
    geom_point(color="blue", alpha=a, pch=16, aes(size=Weight^7)) + 
    theme_classic() + theme(legend.position="none") +
    scale_x_reverse() + 
    stat_smooth(method="glm", aes(weight=Weight), 
                method.args=list(family="binomial"), se=T) +
    xlab("Date (years BP)") + ylab("Score") +
    geom_point(data=avg.df,mapping=aes(x=Date, y=Score), pch=18) + 
    geom_text(data=avg.df, aes(label=Population), size=5, nudge_x = 2000) +
    scale_y_continuous(breaks=c(0,0.5,1)) 
  ggsave(file = paste(out,'_anc.pdf',sep=''), width = 4, height = 4)
  
  # plot ancient humans and 1KGP Europeans
  ggplot(eur.df, aes(x=Date, y=Score)) + 
    geom_point(color="blue", alpha=a, pch=16, aes(size=Weight^7)) + 
    theme_classic() + theme(legend.position="none") +
    scale_x_reverse() + 
    stat_smooth(method="glm",aes(weight=Weight),
                method.args=list(family="binomial"), se=T) +
    xlab("Date (years BP)") + ylab("Score") 
  ggsave(file = paste(out,'_anc+1kgpEur.pdf',sep=''), width = 4, height = 4)

  } else {## Equal weights--logistic regression
    # run simple regression
    anc.df$NumLight <- anc.df$TotalNumSNPs - anc.df$NumDark
    values_anc <- cbind(anc.df$NumDark, anc.df$NumLight)
    colnames(values_anc) <- c('dark', 'light')
    model_anc <- glm(values_anc ~ anc.df$Date, family="binomial")
    
    eur.df$NumLight <- eur.df$TotalNumSNPs - eur.df$NumDark
    values_eur <- cbind(eur.df$NumDark, eur.df$NumLight)
    colnames(values_eur) <- c('dark', 'light')
    model_eur <- glm(values_eur ~ eur.df$Date, family="binomial")
    
    # save model summary
    sink(paste(out,'_anc_simple.txt',sep=''))
    print(summary(model_anc))
    sink()
    
    sink(paste(out,'_anc+1kgpEur_simple.txt',sep=''))
    print(summary(model_eur))
    sink()
    
    # plot only ancient humans
    ggplot(anc.df, aes(x=Date, y=Score)) + 
      geom_point(color="blue", alpha=a, aes(size=TotalNumSNPs^7), pch=16) + 
      theme_classic() + theme(legend.position="none") +
      scale_x_reverse() + 
      # scale_radius(limits = c(0,50)) +
      stat_smooth(method="glm", aes(weight=TotalNumSNPs),
                  method.args=list(family="binomial"), se=T) +
      xlab("Date (years BP)") + ylab("Score") +
      scale_y_continuous(breaks=c(0,0.5,1)) + 
      geom_point(data=avg.df,mapping=aes(x=Date, y=Score), pch=18) + 
      geom_text(data=avg.df, aes(label=Population), size=4, nudge_x = 2000)
    
    ggsave(file = paste(out,'_anc.pdf',sep=''), width = 4, height = 4)
    
    # plot ancient humans and 1KGP Europeans
    ggplot(eur.df, aes(x=Date, y=Score)) + 
      geom_point(color="blue", alpha=a, aes(size=TotalNumSNPs^7), pch=16) + 
      theme_classic() + theme(legend.position="none") +
      scale_x_reverse() + #scale_radius(limits = c(0,50)) +
      stat_smooth(method="glm", aes(weight=TotalNumSNPs),
                  method.args=list(family="binomial"), se=T) +
      xlab("Date (years BP)") + ylab("Score") +
      scale_y_continuous(breaks=c(0,0.5,1))
    ggsave(file = paste(out,'_anc+1kgpEur.pdf',sep=''), width = 4, height = 4)
  }
  
}
