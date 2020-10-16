## plotAncTimePval.R
#` Make plot of ancestry p value against time p value
#` Write table with ancestry and time p-values and SNP annotations
args = commandArgs(trailingOnly = T)

# library(wordcloud)
library(ggplot2)
library(ggrepel)

selscan <- args[1] # coefficients for all genes
gene_path <- args[2] # master file
name <- args[3] # file name for output figures and models
no_filter_flag <- as.logical(args[4])
trait_selscan <- args[5] # Just for figure 2a, 2d

genes <- read.csv(gene_path, as.is=T, header=T, na.strings = '')
selscan <- read.table(selscan, as.is=T, header=T)

if (no_filter_flag) {
  background <- selscan
  selscan <- read.table(trait_selscan, as.is=T, header=T)
}

# get coefficients for trait SNPs
genes$SNP <- paste(genes$CHROM, genes$POS, sep = '_')
df2 <- merge(selscan, genes, by='SNP')

# determine direction of change in terms of light allele
df2$alt_light <- df2$Alt == df2$Light
df2$light_beta <- NA
df2$light_beta[which(df2$alt_light==T)] <- df2$beta_t[which(df2$alt_light==T)]
df2$light_beta[which(df2$alt_light==F)] <- -df2$beta_t[which(df2$alt_light==F)]

df2$light_inc <- 0
df2$light_inc[which(df2$light_beta<0)] <- 1
df2$sqrt_light_beta <- abs(df2$light_beta)^(1/2.5)
df2$sqrt_light_beta[which(df2$light_inc==1)] <- -df2$sqrt_light_beta[which(df2$light_inc==1)]

# log10(p-val)
df2$log10p_anc <- log10(df2$p_anc)
df2$log10p_t <- log10(df2$p_t)

# # set bonferonni 0.05 for number of tests
# p_threshold <- -log10(0.05 / nrow(df2))

# find quantiles based on genome-wide distribution
if (no_filter_flag) {
  time_q <- quantile(background$p_t, .05, na.rm=T)
  anc_q <- quantile(background$p_anc, .05, na.rm=T)
} else {
  time_q <- quantile(selscan$p_t, .05, na.rm=T)
  anc_q <- quantile(selscan$p_anc, .05, na.rm=T)
}

# Take top 5 ancestry and time SNPs
anc_th <- quantile(df2$p_anc,0.03, na.rm=T)
time_th <- quantile(df2$p_t,0.03, na.rm=T)

if (nrow(df2)<20) {
  top_snps <- which(!(is.na(df2$p_t)))
} else {
  top_snps <- which(df2$p_anc < anc_th | df2$p_t < time_th)
}

## Plot ancestry and time p-values figure
b <- ggplot(df2, aes(x=-log10(p_t), y=-log10(p_anc))) +
  theme_classic(base_size = 20) +
  labs(x=bquote('Time ['~-log[10]~'(p-value)]'), y=bquote('Ancestry ['~-log[10]~'(p-value)]')) +
  theme(text=element_text(family='Helvetica'))

b + geom_point(aes(size=1.5, colour=sqrt_light_beta)) +
  theme(legend.position='none')  +
  geom_label_repel(data=df2[top_snps, ],
                   aes(label=NearestGene), size=5,
                   force=5, segment.alpha = 0.7, min.segment.length = 0.2,
                   arrow=arrow(type='open', length = unit(0.2, "cm")),
                   fill=rgb(1,1,1,0.5)) +
  geom_vline(aes(xintercept = -log10(time_q)), linetype='dashed', colour='red') +
  geom_hline(aes(yintercept = -log10(anc_q)), linetype='dashed', colour='red') +
  scale_color_gradient2(low='red', mid='grey', high='blue')

file=paste(name, ".pdf", sep="")
ggsave(filename=file, width = 4.5, height = 5)

 # remove rows with NA
filter_na <- is.na(df2$p_t)
out.df <- df2[!(filter_na), ]

## Export table with summary stats and SNP annotations
write.table(out.df, file=paste(name, ".csv", sep=""), quote=F, row.names=F,
            sep=',')


## Examine distribution of trait SNPs relative to genome-wide
# Remove pigmentation snps from selscan
background <- selscan[!(selscan$SNP %in% df2$SNP), ]
# time term
pdf(file=time_out ,width = 3, height = 3)
par(mar=c(4,4,1,2))
qqplot(-log10(background$p_t), -log10(df2$p_t),
       pch = 16, frame = T, cex=0.2, col='black',
       ylab='Skin pigmentation SNPs', xlab='Genome-wide'#, xlim=c(0,15)
)
abline(0,1, col=rgb(0, 0, 0, max=255,alpha=100))
dev.off()

wilcox.test(x=background$p_t, y=df2$p_t)

# ancestry term
pdf(file=anc_out ,width = 3, height = 3)
par(mar=c(4,4,1,2))
qqplot(-log10(background$p_anc), -log10(df2$p_anc),
       pch = 16, frame = T, cex=0.2, col='black',
       ylab='Skin pigmentation SNPs', xlab='Genome-wide'
)
abline(0,1, col=rgb(0, 0, 0, max=255,alpha=100))
dev.off()

wilcox.test(x=background$p_anc, y=df2$p_anc)
