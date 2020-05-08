## Plot Qx statistics

library(ggplot2)

# CSV file with Qx statistic output with columns 'Set','Qx_pval','log10p'
qx_all <- '/Volumes/work/pigment/qst/qx_stats_v3.csv'

out_file_path <- '/Users/Dan/Documents/OneDrive/MathiesonLab/paper-draft/figures/fig4b.pdf'

qx.df <- read.csv(qx_all, as.is=T)

qx.df <- qx.df[!(is.na(qx.df[,1])), ]

qx.df$Set <- factor(qx.df$Set,
                    levels = c('5','6','7','8','9','10','15','20','25','30','40'
                               ,'50','75','100','170'))

qx.df$log10p <- -log10(as.numeric(qx.df$Qx_pval))

# plot just GWAS effect size
qx.df$ordering <- 'gwas'
ggplot(data=qx.df,
       aes(x=Set,y=log10p, group=ordering)) +
  geom_point(aes(x=Set,y=log10p), size=1) +
  geom_line(lwd=0.4) +
  theme_classic(base_size = 9) +
  theme(text=element_text(family='Helvetica'),
        axis.title=element_text(size=13)) +
  labs(x='Number of SNPs',
       y=bquote(''~-log[10]~'('~Q[x]~ 'p-value)'))

ggsave(filename = out_file_path,
       height=2.8, width=3)
