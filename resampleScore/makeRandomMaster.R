## makeRandomMaster.R
## Make master file for random SNPs genome-wide 

library(data.table)

args = commandArgs(trailingOnly = TRUE)
control_path <- args[1]
master_path <- args[2]
all_path <- args[3]
label_path <- args[4]
use_beta <- args[6]

# read data
con <- scan(control_path, what = character())
master.df <- read.csv(master_path, header=T, as.is=T)
label <- scan(label_path, what = character())
all.df <- fread(all_path, stringsAsFactors = F)

# give header to all.df
colnames(all.df) <- c('random','rref','ralt','ranc','rder','rderfq')
# get rid of "snp" in first entry of label
label <- label[2:length(label)]

# combine label and random set for merging
control_mat <- rbind(label, con)
c.df <- data.frame(t(control_mat), stringsAsFactors = F)
colnames(c.df) <- c('snp','random')

# merge with all.df
master.df$snp <- paste(master.df$CHROM, master.df$POS, sep='_')
cm.df <- merge(c.df, master.df, by='snp')
cma.df <- merge(cm.df, all.df, by='random')

# assign random snps dark/light allele
i <- which(cma.df$Dark == cma.df$Der)
j <- which(cma.df$Dark != cma.df$Der)

cma.df$rdark <- NA
cma.df$rlight <- NA

cma.df$rdark[i] <- cma.df$rder[i]
cma.df$rdark[j] <- cma.df$ranc[j]
  
cma.df$rlight[i] <- cma.df$ranc[i]
cma.df$rlight[j] <- cma.df$rder[j]

# obtain summary statistics from UKB 
if (use_beta==T | use_beta=="true" | use_beta=='T') {
  cma.df$rbeta <- abs(cma.df$Beta)
  k <- which(cma.df$rdark == cma.df$ralt) 
  l <- which(cma.df$rdark != cma.df$ralt)
  cma.df$rbeta[k] <- cma.df$Beta[k]
  cma.df$rbeta[l] <- -cma.df$rbeta[l]
  cma.df$p <- NA
  cma.df$pop <- NA
} else { # put missing values
  cma.df$rbeta <- 0
  cma.df$SE <- 0
  cma.df$p <- NA
  cma.df$pop <- NA
}

# fill in extraneous columns
cma.df$minor <- NA
cma.df$major <- NA
cma.df$near <- NA
cma.df$source <- 'random'

# separate chr_pos
variant <- strsplit(cma.df$random, split = "_")
cma.df$chr <- lapply(variant, function(x){x[1]})
cma.df$pos <- lapply(variant, function(x){x[2]})
cma.df$chr <- unlist(cma.df$chr)
cma.df$pos <- unlist(cma.df$pos)
final.df <- cma.df[ , c('chr','pos','random','rdark','rlight','rref','ralt','ranc',
                        'rder','minor','major','rbeta','SE','p','near','pop',
                        'source')]

col_names <- c('CHROM', 'POS', 'RSID', 'Dark', 'Light','Ref', 'Alt', 'Anc', 
               'Der', 'Minor', 'Major', 'Beta', 'SE',	'p', 'NearestGene', 'Pop',
               'Source')

colnames(final.df) <- col_names
  
write.table(final.df, args[5], sep=',', quote=F, row.names=F)
