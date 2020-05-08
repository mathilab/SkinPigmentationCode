## Determine percentile of given SNP beta coefficient from regression model
#` in genome-wide distribution of betas
#` Note: generate model coefficents beforehand

library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
ctrl_path <- args[1]
b_input <- as.numeric(args[2])
t_input <- as.numeric(args[3])
p_input <- as.numeric(args[4])
file <- args[5]

# load data
ctrl.df = read.table(ctrl_path, as.is=T, header=F)
colnames(ctrl.df) <- c('beta','se','p','n')

ctrl.df$t <- ctrl.df$beta / ctrl.df$se

# determine percentile for beta, t, p 
perc_b <- ecdf(ctrl.df$beta)
b <- perc_b(b_input)

perc_t <- ecdf(ctrl.df$t)
t <- perc_t(t_input)

perc_p <- ecdf(ctrl.df$p)
p <- perc_p(p_input)

# save to a text file
sink(paste(file,'_percent.txt',sep=''))
cat(paste('beta:',b, '\n', sep=" "))
cat(paste('T-stat:', t, '\n', sep=" "))
cat(paste('p-value:', p, '\n', sep=" "))
sink()

# make histograms 
w <- 2.5
pdf(file=paste(file,'_hist_beta.pdf',sep=''), width = w, height = 2)
ggplot(ctrl.df, aes(x=beta)) + geom_histogram(aes(y=..density..)) +
  geom_density(alpha=.2, fill="#FF6666") + 
  geom_vline(xintercept=b_input)
dev.off()

## S2c (Clean up the x-axis)
# pdf(file=paste(file,'_hist_beta.pdf',sep=''), width = w, height = 2)
# ggplot(ctrl.df, aes(x=beta)) + geom_histogram(aes(y=..density..)) +
#   geom_density(alpha=.2, fill="#FF6666") +
#   geom_vline(xintercept=b_input) +
#   scale_x_continuous(breaks=c(-3E-5,0,4e-5))
# dev.off()

# S2E (Clean up the x-axis)
# pdf(file=paste(file,'_hist_beta.pdf',sep=''), width = w, height = 2)
# ggplot(ctrl.df, aes(x=beta)) + geom_histogram(aes(y=..density..)) +
#   geom_density(alpha=.2, fill="#FF6666") +
#   geom_vline(xintercept=b_input) +
#   scale_x_continuous(breaks=c(-3E-5,0,5e-5))
# dev.off()

# Supplement unweighted UKB 1240K 40kyBP scores
# pdf(file=paste(file,'_hist_beta.pdf',sep=''), width = w, height = 2)
# ggplot(ctrl.df, aes(x=beta)) + geom_histogram(aes(y=..density..)) +
#   geom_density(alpha=.2, fill="#FF6666") +
#   geom_vline(xintercept=b_input) +
#   scale_x_continuous(breaks=c(-1E-5,0,1E-5))
# dev.off()

pdf(file=paste(file,'_hist_t.pdf',sep=''), width = w, height = 2)
ggplot(ctrl.df, aes(x=t)) + geom_histogram(aes(y=..density..))  +
  geom_density(alpha=.2, fill="#FF6666") + 
  geom_vline(xintercept=t_input)
dev.off()

pdf(file=paste(file,'_hist_p.pdf',sep=''), width = w, height = 2)
ggplot(ctrl.df, aes(x=p)) + geom_histogram(aes(y=..density..))  +
  geom_density(alpha=.2, fill="#FF6666") + 
  geom_vline(xintercept=p_input)
dev.off()

# make q-q plot of -log10(p)
pdf(file=paste(file,'_qq_p.pdf',sep=''), width = 4, height = 3.5)
par(mar=c(4,4,1,2))
qqplot(-log10(runif(nrow(ctrl.df))), -log10(ctrl.df$p),
       pch=16, cex=1, 
       xlab='Theoretical', ylab='Empirical')
abline(0,1)
dev.off()

qqplot(ctrl.df$p, runif(nrow(ctrl.df)), pch=16,
       xlab='Theoretical', ylab='Empirical', cex=1,
       cex.lab=1)
abline(0,1)
