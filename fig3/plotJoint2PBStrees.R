## Plot 2 PBS distributions against each other from different trees
#` Note: This is run manually. Must generate PBS distributions beforehand
#` for the branches to be plotted against each other.

library(wordcloud)
library(ggplot2)
library(ggrepel)
library(data.table)

## Main paths (EDIT USER-SPECIFIED PATHS)
## __________________________________________________________________________
# Master file with UKB SNPs
master_path <- '/Users/danju/desktop/pigment/ukbio/runs2/p_5e-8_1240k_T/clump_true_r2_0.05/ukb_master_v37_15ky_ukb_skin_r2_0.05_p_5e-8.csv'
# Ancestry and time p-vals for individual SNP regressions
pval_path <- '/Users/danju/Desktop/pigment/ukbio/runs2/p_5e-8_1240k_T/clump_true_r2_0.05/anc-time/v37_15ky_ukb_skin_r2_0.05_p_5e-8.csv'
base_dir <- '/Users/danju/Desktop/pigment'

## __________________________________________________________________________

plotFig3 <- function(gen_path,trait_path,genes_path,fq_path,g1_name,g2_name,out_path) {
  gen.df <- read.table(gen_path, as.is=T, header=T)
  trait.df <- read.table(trait_path, as.is=T, header=T)
  fq.df <- fread(fq_path, header=T)
  genes.df <- read.table(genes_path, as.is=T)
  pval.df <- read.csv(pval_path, as.is=T)

  master.df <- read.csv(master_path, as.is=T, header=T)
  master.df$ID <- paste(master.df$CHROM, master.df$POS, sep='_')
  # Add gene names to trait df if names are known
  trait.df2 <- merge(trait.df, master.df[ , c('ID','NearestGene','Alt','Light')], by='ID')

  # Add genes to gen.df
  gen.df$nearest <- genes.df$V7

  # Take top trait loci in the groups
  anc_pbs_th <- quantile(gen.df$pbs_anc, 0.99, na.rm=T)
  eur_pbs_th <- quantile(gen.df$pbs_eur, 0.99, na.rm=T)

  top_snps <- trait.df2$ID[which(trait.df2$pbs_anc > anc_pbs_th | trait.df2$pbs_eur > eur_pbs_th)]

  # Get ancestry driven SNPs (pass Bonferoni corrected p-value)
  anc_filt <- which(pval.df$p_anc < 0.05/nrow(pval.df))
  anc_snps <- pval.df$SNP[anc_filt]

  # Get time driven SNPs (pass Bonferoni corrected p-value)
  time_filt <- which(pval.df$p_t < 0.05/nrow(pval.df))
  time_snps <- pval.df$SNP[time_filt]

  # Get 99.9, 99, 80 quantiles for genome-wide PBS for both groups
  g2_q1 <- quantile(gen.df$pbs_eur, .999, na.rm=T)
  g2_q2 <- quantile(gen.df$pbs_eur, .99, na.rm=T)
  g2_q3 <- quantile(gen.df$pbs_eur, .8, na.rm=T)

  g1_q1 <- quantile(gen.df$pbs_anc, .999, na.rm=T)
  g1_q2 <- quantile(gen.df$pbs_anc, .99, na.rm=T)
  g1_q3 <- quantile(gen.df$pbs_anc, .8, na.rm=T)

  ## Calculate difference in light allele frequency in European and Ancient
  # trait.fq.df <- merge(trait.df2, fq.df[ , c('ID','anc_fq','eur_fq')], by='ID')
  # trait.fq.df$g1_light_fq[trait.fq.df$Alt==trait.fq.df$Light] <- trait.fq.df$anc_fq[trait.fq.df$Alt==trait.fq.df$Light]
  # trait.fq.df$g1_light_fq[trait.fq.df$Alt!=trait.fq.df$Light] <- 1 - trait.fq.df$anc_fq[trait.fq.df$Alt!=trait.fq.df$Light]
  # trait.fq.df$g2_light_fq[trait.fq.df$Alt==trait.fq.df$Light] <- trait.fq.df$eur_fq[trait.fq.df$Alt==trait.fq.df$Light]
  # trait.fq.df$g2_light_fq[trait.fq.df$Alt!=trait.fq.df$Light] <- 1 - trait.fq.df$eur_fq[trait.fq.df$Alt!=trait.fq.df$Light]
  # trait.fq.df$light_diff <- trait.fq.df$g1_light_fq - trait.fq.df$g2_light_fq

  ## Calculate difference in light allele frequency in European and YRI
  # trait.fq.df <- merge(trait.df2, fq.df[ , c('ID','out_fq','eur_fq')], by='ID')
  # trait.fq.df$g1_light_fq[trait.fq.df$Alt==trait.fq.df$Light] <- trait.fq.df$out_fq[trait.fq.df$Alt==trait.fq.df$Light]
  # trait.fq.df$g1_light_fq[trait.fq.df$Alt!=trait.fq.df$Light] <- 1 - trait.fq.df$out_fq[trait.fq.df$Alt!=trait.fq.df$Light]
  # trait.fq.df$g2_light_fq[trait.fq.df$Alt==trait.fq.df$Light] <- trait.fq.df$eur_fq[trait.fq.df$Alt==trait.fq.df$Light]
  # trait.fq.df$g2_light_fq[trait.fq.df$Alt!=trait.fq.df$Light] <- 1 - trait.fq.df$eur_fq[trait.fq.df$Alt!=trait.fq.df$Light]
  # trait.fq.df$light_diff <- trait.fq.df$g1_light_fq - trait.fq.df$g2_light_fq

  ## Calculate difference in light allele frequency in Ancient and YRI (Used in paper)
  trait.fq.df <- merge(trait.df2, fq.df[ , c('ID','out_fq','anc_fq')], by='ID')
  trait.fq.df$g1_light_fq[trait.fq.df$Alt==trait.fq.df$Light] <- trait.fq.df$out_fq[trait.fq.df$Alt==trait.fq.df$Light]
  trait.fq.df$g1_light_fq[trait.fq.df$Alt!=trait.fq.df$Light] <- 1 - trait.fq.df$out_fq[trait.fq.df$Alt!=trait.fq.df$Light]
  trait.fq.df$g2_light_fq[trait.fq.df$Alt==trait.fq.df$Light] <- trait.fq.df$anc_fq[trait.fq.df$Alt==trait.fq.df$Light]
  trait.fq.df$g2_light_fq[trait.fq.df$Alt!=trait.fq.df$Light] <- 1 - trait.fq.df$anc_fq[trait.fq.df$Alt!=trait.fq.df$Light]
  trait.fq.df$light_diff <- trait.fq.df$g1_light_fq - trait.fq.df$g2_light_fq

  # Get top 5 ancestry and time SNPs
  anc_th <- quantile(pval.df$p_anc,0.03, na.rm=T)
  time_th <- quantile(pval.df$p_t,0.03, na.rm=T)

  time_series_snps <- pval.df$SNP[which(pval.df$p_anc < anc_th | pval.df$p_t < time_th)]

  name_snps <- c(time_series_snps, top_snps)

  # Make joint distribution plot
  b <- ggplot(gen.df, aes(x=pbs_eur, y=pbs_anc)) +
    theme_classic(base_size = 15) +
    labs(x=paste(g2_name, 'PBS', sep=' '),
         y=paste(g1_name, 'PBS', sep=' ')) +
    theme(text=element_text(family='Helvetica', size = 18))

  b + geom_point(aes(alpha=0.02),colour="grey",size=1, pch=16) +
    theme(legend.position='none') +
    geom_point(data=trait.fq.df, aes(x=pbs_eur, y=pbs_anc, fill=light_diff),
               size=3, shape=22, alpha=0.8) +
    scale_fill_gradient2(low='red', mid='white', high='darkblue') +
    geom_vline(xintercept = g2_q2, linetype="dashed", color = "orange") + #99%
    geom_vline(xintercept = g2_q1, linetype="dashed", color = "red") +#99.9%
    geom_hline(yintercept = g1_q1, linetype="dashed", color = "red") + #99.9%
    geom_hline(yintercept = g1_q2, linetype="dashed", color = "orange") + #99%
    geom_label_repel(data=trait.df2[trait.df2$ID %in% name_snps, ],
                     aes(label=NearestGene), size=4, segment.size=0.2,
                     force=25, segment.alpha = 0.7, min.segment.length = 0.1,
                     arrow=arrow(type='open', length = unit(0.2, "cm")),
                     fill=rgb(1,1,1,0.5)) +
    geom_point(data=trait.fq.df[trait.fq.df$ID %in% anc_snps, ],
               aes(x=pbs_eur, y=pbs_anc),
               shape=4, size=1) +
    geom_point(data=trait.fq.df[trait.fq.df$ID %in% time_snps, ],
               aes(x=pbs_eur, y=pbs_anc),
               shape=3, size=1) +
    xlim(-0.1, 1.1) + ylim(c(-0.1, 1))

  ggsave(out_path, width = 4.5, height = 4.5)


  ## Test if trait SNPs differ from background
  # wilcox.test(x=gen.df$pbs_anc, y=trait.df$pbs_anc)
  #
  # pdf(file=out ,width = 3, height = 3)
  # par(mar=c(4,4,1,2))
  # qqplot(gen.df$pbs_anc, trait.df$pbs_anc,
  #        pch = 16, frame = T, cex=0.5, col='black',
  #        ylab='Skin pigmentation SNPs', xlab='Genome-wide',
  #        xlim=c(0,1), ylim=c(0,1)
  # )
  # abline(0,1, col=rgb(0, 0, 0, max=255,alpha=100))
  # dev.off()
}

## Make plots
# EF
gen_path <- '/Users/danju/Desktop/pigment/pbs/ef-joint/pbs_table.tsv'
trait_path <- '/Users/danju/Desktop/pigment/pbs/ef-joint/trait_pbs.tsv'
genes_path <- '/Users/danju/Desktop/pigment/pbs/ef-joint/closest.tsv'
fq_path <- '/Users/danju/Desktop/pigment/pbs/ef-joint/all_fq.tsv'
g1_name <- 'EF'
g2_name <- 'GBR'
out_path <- '/Users/danju/OneDrive/MathiesonLab/remote/fig3/ef_gbr_pbs.pdf'
plotFig3(gen_path,trait_path,genes_path,fq_path,g1_name,g2_name,out_path)

# HG
gen_path <- '/Users/danju/Desktop/pigment/pbs/hg-joint/pbs_table.tsv'
trait_path <- '/Users/danju/Desktop/pigment/pbs/hg-joint/trait_pbs.tsv'
genes_path <- '/Users/danju/Desktop/pigment/pbs/hg-joint/closest.tsv'
fq_path <- '/Users/danju/Desktop/pigment/pbs/hg-joint/all_fq.tsv'
g1_name <- 'HG'
g2_name <- 'GBR'
out_path <- '/Users/danju/OneDrive/MathiesonLab/remote/fig3/hg_gbr_pbs.pdf'
plotFig3(gen_path,trait_path,genes_path,fq_path,g1_name,g2_name,out_path)

# SP
gen_path <- '/Users/danju/Desktop/pigment/pbs/sp-joint/pbs_table.tsv'
trait_path <- '/Users/danju/Desktop/pigment/pbs/sp-joint/trait_pbs.tsv'
genes_path <- '/Users/danju/Desktop/pigment/pbs/sp-joint/closest.tsv'
fq_path <- '/Users/danju/Desktop/pigment/pbs/sp-joint/all_fq.tsv'
g1_name <- 'SP'
g2_name <- 'GBR'
out_path <- '/Users/danju/OneDrive/MathiesonLab/remote/fig3/sp_gbr_pbs.pdf'
plotFig3(gen_path,trait_path,genes_path,fq_path,g1_name,g2_name,out_path)

## Supplemental figures
# # EF - 20 SNP windows
# gen_path <- '/Users/danju/Desktop/pigment/pbs/ef-joint-20/pbs_table.tsv'
# trait_path <- '/Users/danju/Desktop/pigment/pbs/ef-joint-20/trait_pbs.tsv'
# genes_path <- '/Users/danju/Desktop/pigment/pbs/ef-joint-20/closest.tsv'
# fq_path <- '/Users/danju/Desktop/pigment/pbs/ef-joint-20/all_fq.tsv'
# g1_name <- 'EF'
# g2_name <- 'GBR'
# out_path <- '/Users/danju/desktop/ef_gbr_pbs_20.pdf'
# # SP - 20 SNP windows
# gen_path <- '/Users/danju/Desktop/pigment/pbs/sp-joint-20/pbs_table.tsv'
# trait_path <- '/Users/danju/Desktop/pigment/pbs/sp-joint-20/trait_pbs.tsv'
# genes_path <- '/Users/danju/Desktop/pigment/pbs/sp-joint-20/closest.tsv'
# fq_path <- '/Users/danju/Desktop/pigment/pbs/sp-joint-20/all_fq.tsv'
# g1_name <- 'SP'
# g2_name <- 'GBR'
# out_path <- '/Users/danju/desktop/sp_gbr_pbs_20.pdf'
# # HG - 20 SNP windows
# gen_path <- '/Users/danju/Desktop/pigment/pbs/hg-joint-20/pbs_table.tsv'
# trait_path <- '/Users/danju/Desktop/pigment/pbs/hg-joint-20/trait_pbs.tsv'
# genes_path <- '/Users/danju/Desktop/pigment/pbs/hg-joint-20/closest.tsv'
# fq_path <- '/Users/danju/Desktop/pigment/pbs/hg-joint-20/all_fq.tsv'
# g1_name <- 'HG'
# g2_name <- 'GBR'
# out_path <- '/Users/danju/desktop/hg_gbr_pbs_20.pdf'
