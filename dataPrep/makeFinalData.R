## Prepare final dataset for analysis from Reich Lab v37

library(ggplot2)

######## Paths ########
dir <- '/Users/danju/Desktop/pigment/ancientEuroVCF/v37/v37.2.1240K'
# folder containing eigenstrat files
archaic <- '/Users/danju/Desktop/pigment/misc/v37_archaic_instance_id.txt'
# file listing archaic sample Instance ID to exclude
map_save_path <- '/Users/danju/Desktop/sfig1.pdf'
# export path for map of samples included/excluded
#######################

## Prep ancient datasets
setwd(dir)
#` read data
anno.df <- read.csv('v37.2.1240K.clean4.anno.csv', as.is=T)
arch_samples <- scan(archaic, what = "character")
#` fix column classes in anno.df
anno.df$lat <- as.numeric(anno.df$Lat.)
anno.df$long <- as.numeric(anno.df$Long.)
anno.df$date <- as.numeric(anno.df$Average.of.95.4..date.range.in.calBP..defined.as.1950.CE.)
anno.df$Coverage <- as.numeric(anno.df$Coverage)
#` exclude present-day samples
samples.df <- anno.df[anno.df$date > 0, ]
#` exclude archaic individuals
samples.df <- samples.df[!samples.df$Instance.ID %in% arch_samples, ]
#` exclude low coverage (<0.1x)
samples.df <- samples.df[samples.df$Coverage>=0.1, ]
## Remove 390k data type samples
samples.df <- samples.df[samples.df$Data.type != '390K', ]
#` remove samples outside of broad Eurasia zone
lat_l = 0
lat_h = 75
lon_l = -15
lon_h = 107
samples.df <- samples.df[samples.df$lat < lat_h & samples.df$lat > lat_l & samples.df$long>lon_l & samples.df$long < lon_h, ]
samples.df <- samples.df[!is.na(samples.df$Instance.ID), ]
#` subset data for samples included in analysis
lat_l2 = 35
lat_h2 = 75
lon_l2 = -15
lon_h2 = 60
include.df <- samples.df[samples.df$lat < lat_h2 & samples.df$lat > lat_l2 & samples.df$long>lon_l2 & samples.df$long < lon_h2, ]
#` remove samples in Middle-East and West of Urals
botai <- which(include.df$lon>55 & include.df$lat<50)
include.df <- include.df[!(include.df$Instance.ID %in% include.df$Instance.ID[botai]), ]
include.df <- include.df[include.df$Country != 'Iran', ]
#` subset data for samples exlcuded in analysis
exclude.df <- samples.df[!(samples.df$Instance.ID %in% include.df$Instance.ID), ]
#` pick map region for West Eurasia
world.df <- map_data('world')
eur <- world.df[world.df$lat < lat_h & world.df$lat > lat_l & world.df$long>lon_l & world.df$long < lon_h, ]
#` Make map of included/excluded samples
ggplot() + theme_minimal() + theme(legend.position = 'None') +
  geom_polygon(data = eur, aes(x=long, y = lat, group = group), fill="grey", alpha=0.3) +
  coord_map() + geom_polygon(data = eur, aes(x = long, y = lat, group=group), fill='grey', alpha=0.3, color = "white", lwd=0.5) + geom_point(data = include.df, aes(x = long, y = lat), size=1.5, shape=1, alpha=0.7) +
  geom_point(data = exclude.df, aes(x = long, y = lat), size=1.5, shape=13, alpha=0.7) +
  scale_x_continuous(breaks = round(seq(min(eur$long), max(eur$long), by = 5),1)) +
  scale_y_continuous(breaks = round(seq(min(eur$lat), max(eur$lat), by = 5),1))

ggsave(map_save_path)



#` Save filtered dataset
write.table(include.df, file="v37_1240k_anno_analysis.csv", row.names = F, quote=F, sep=',')
## Prepare eig files to check for duplication
#` Eigen tools must be installed
dir.create("work")

## Change ind file populations to Ignore
ind.df <- read.table('v37.2.1240K.ind', as.is=T)
ind.df$V3[!(ind.df$V1 %in% include.df$Instance.ID)] <- 'Ignore'
write.table(ind.df, file="v37.2.1240K_ignore.ind", row.names = F, quote=F, col.names = F)

#` Make par file
sink(file = "work/eig2ped.txt")
cat("genotypename:  v37.2.1240K.geno\n")
cat("snpname: v37.2.1240K.snp\n")
cat("indivname: v37.2.1240K_ignore.ind\n")
cat("outputformat:  PED\n")
cat("genotypeoutname: work/ancient.ped\n")
cat("snpoutname:  work/ancient.snp\n")
cat("indivoutname:  work/ancient.ind\n")
sink()

#` Convert format from eig to PLINK .bed
system("convertf -p work/eig2ped.txt")
system("awk '{print $1,\"\t\",$2,\"\t\",$3,\"\t\",$4}' work/ancient.snp > work/ancient.map")

#` Calculate IBS pairwise
system("plink --file work/ancient --distance square gz ibs")
system("plink --file work/ancient --missing gz")

#` Check and remove duplicates
ibs <- read.table("plink.mibs.gz", as.is=TRUE, header=FALSE)
ids <- read.table("plink.mibs.id", as.is=TRUE)
missing <- read.table("plink.imiss.gz", as.is=TRUE, header=TRUE)
rownames(missing) <- missing$IID

ibs <- as.matrix(ibs)
rownames(ibs) <- colnames(ibs) <- ids[,2]

## Remove duplicate samples
duplicates <- data.frame(KEEP=character(), REMOVE=character(), IBS=numeric(), MISS.K=numeric(), MISS.R=numeric(), stringsAsFactors=FALSE)
di <- 1
max_ibs <- c()
max_name1 <- c()
max_name2 <- c()
keep_vec <- c()
mean_vec <- c()
sd_vec <- c()
## Iain's previous code to remove duplicates (use this as a reference)
for(i in 1:(NROW(ibs)-1)){
  this.row <- ibs[i,(i+1):NCOL(ibs)]
  mm <- mean(ibs[i,])
  cutoff <- 1-0.55*(1-mm)
  dup <- which(this.row>cutoff)+i

  # get max IBS to identify a cutoff point
  max_ibs <- append(max_ibs, max(this.row))
  max_name1 <- append(max_name1, ids$V2[(i+which(ibs[i,(i+1):NCOL(ibs)]>=max(this.row)))])
  max_name2 <- append(max_name2, ids$V2[i])
  mean_vec <- append(mean_vec, mean(ibs[i,]))
  sd_vec <- append(sd_vec, sd(ibs[i,]))
  if(missing$F_MISS[missing$IID==max_name1[i]] <= missing$F_MISS[missing$IID==max_name2[i]]) {
    keep_vec <- append(keep_vec, 1)
    } else {
    keep_vec <- append(keep_vec, 0)
    }

  if(length(dup)>1)
  {
    cat(paste0("Multiple duplicates ", rownames(ibs)[i], " (", length(dup), ")\n"))
  }else if(length(dup)==1)
  {
    for(d in dup){
      m1 <- missing[rownames(ibs)[i],"F_MISS"]
      m2 <- missing[rownames(ibs)[d],"F_MISS"]
      if(m1<=m2){
        duplicates[di,] <-  c(rownames(ibs)[i], colnames(ibs)[d], ibs[i,d], m1, m2)
      }else{
        duplicates[di,] <-  c( colnames(ibs)[d],rownames(ibs)[i], ibs[i,d], m2, m1)
      }
      di=di+1
    }
  }
}

## Check names of sampled with high IBS values with Iain's previous code
max.df <- data.frame(max_ibs, max_name1, max_name2, keep_vec, mean_vec, sd_vec,
                     stringsAsFactors = F)

# create a normalized max IBS
max.df$max_ibs_norm <- max.df$max_ibs / max.df$mean_vec
max.df$ibs_z <- (max.df$max_ibs - max.df$mean_vec) / max.df$sd_vec

## Find all samples with a familial relation based on annotation
family_terms <- c('brother', 'child', 'parent', 'mother', 'daughter', 'son',
                  'sister', 'father', 'parents')
fam_samples <- grep(paste(family_terms,collapse="|"), anno.df$Group.ID )
anno_fam.df <- anno.df[fam_samples, ]
fam_names <- anno.df$Instance.ID[fam_samples]
max_fam.df <- max.df[max.df$max_name2 %in% fam_names | max.df$max_name1 %in% fam_names, ]

## Look at I1903 family
fam1 <- c('I1903', 'I10619', 'I10624' , 'I10620', 'I10621')
fam1.df <- anno.df[anno.df$Instance.ID %in% fam1, ]
# The related individuals not in there

# Histogram of max IBS
hist(max.df$max_ibs_norm, breaks=20)
hist(max_fam.df$max_ibs_norm)
hist(max.df$ibs_z, breaks=20)
hist(max_fam.df$ibs_z, breaks=12)

# Assign which sample to keep/remove based on missing sites
max.df$keep[max.df$keep_vec==1] <- max.df$max_name1[max.df$keep_vec==1]
max.df$keep[max.df$keep_vec==0] <- max.df$max_name2[max.df$keep_vec==0]
max.df$remove[max.df$keep_vec==1] <- max.df$max_name2[max.df$keep_vec==1]
max.df$remove[max.df$keep_vec==0] <- max.df$max_name1[max.df$keep_vec==0]

maxkeep.df <- max.df[max.df$ibs_z>7, ]
maxkeep.df$remove %in% duplicates$REMOVE
dup_names_remove <- maxkeep.df$remove
## Add Loschbour to exclusion list b/c no pseudo-haploid
dup_names_remove <- append(dup_names_remove, 'Loschbour_published.DG')
## Switch Stuttgart_published.DG out to be replaced by I0018
dup_names_remove <- dup_names_remove[dup_names_remove != 'I0018']
dup_names_remove <- append(dup_names_remove, 'Stuttgart_published.DG')
## export list of samples to remove
file_name_duplicates <- '/Users/danju/desktop/pigment/misc/v37_duplicates.txt'
write.table(dup_names_remove, file=file_name_duplicates, row.names = F,
            quote = F, col.names = F)

## Choose I0018 over Stuttgart_published.DG b/c pseudo-haploid
grep('dg', maxkeep.df$keep, ignore.case=T) # which samples are diploid
overlap <- maxkeep.df$keep %in% duplicates$KEEP
length(which(overlap==T))

############################################################################
## I manually went through the max_fam.df to choose which samples to exclude
#` based on annotations and coverage.  <v37_family_exclude_iid.txt>
############################################################################

# # Look for potentially familially related samples
# fam.df <- max.df[max.df$max_ibs>0.78, ]
# fam.df <- fam.df[!(fam.df$keep%in%maxkeep.df$keep), ]
# # Filter annotation file
# col_keep <- c('Instance.ID', 'Master.ID', 'Data.type', 'Publication', 'Average.of.95.4..date.range.in.calBP..defined.as.1950.CE.', 'Group.ID', 'Location', 'Country', 'Lat.', 'Long.', 'Coverage')
# family_anno.df <- anno.df[anno.df$Instance.ID %in% fam.df$remove, col_keep]
# colnames(fam.df)[6] <- 'Instance.ID'
# family_anno.df <- merge(family_anno.df, fam.df[ , c(1,5,6)] , by = 'Instance.ID')

## Create final dataset!
exclude_family <- '/Users/danju/desktop/pigment/misc/v37_family_exclude_iid.txt'

duplicate_list <- scan(file_name_duplicates, what = 'character')
family_list <- scan(exclude_family, what = 'character')
exclude_samples <- append(duplicate_list, family_list)

## remove duplicates and family members from include.df
final.df <- include.df[!(include.df$Instance.ID %in% exclude_samples), ]

final_filename <- '/Users/danju/desktop/pigment/ancientEuroVCF/v37/v37.2.1240K/final_ancient_capture_anno.tsv'
write.table(final.df, file=final_filename, row.names=F, quote=F, sep='\t')

## Filter the eigen files with convertf
## Change ind file populations to Ignore
ind.df <- read.table('v37.2.1240K.ind', as.is=T)
ind.df$V3[!(ind.df$V1 %in% final.df$Instance.ID)] <- 'Ignore'
write.table(ind.df, file="v37.2.1240K_ignore2.ind", row.names = F, quote=F, col.names = F)

#` Make par file
sink(file = "work2/eig2ped.txt")
cat("genotypename:  v37.2.1240K.geno\n")
cat("snpname: v37.2.1240K.snp\n")
cat("indivname: v37.2.1240K_ignore2.ind\n")
cat("outputformat:  PED\n")
cat("genotypeoutname: work2/ancient.ped\n")
cat("snpoutname:  work2/ancient.snp\n")
cat("indivoutname:  work2/ancient.ind\n")
sink()

system("convertf -p work2/eig2ped.txt")

## Convert back to eigen format
#` Make par file
sink(file = "final/ped2eig.txt")
cat("genotypename:  work2/ancient.ped\n")
cat("snpname: work2/ancient.pedsnp\n")
cat("indivname: work2/ancient.pedind\n")
cat("outputformat:  EIGENSTRAT\n")
cat("genotypeoutname: final/v37.2_1240k_filtered.geno\n")
cat("snpoutname:  final/v37.2_1240k_filtered.snp\n")
cat("indivoutname:  final/v37.2_1240k_filtered.ind\n")
cat("familynames:     NO")
sink()

system("convertf -p final/ped2eig.txt")

## Merge v37 with new data from: Villalba-MoucoCurrentBiology2019 + OlaldeNature2019
setwd("final/")
#` Make par file for first merge
add1_path <- '/Users/danju/Desktop/pigment/ancientEuroVCF/Villalba-MoucoCurrentBiology2019/'
sink(file = "wd1/merge1.txt")
cat("geno1: v37.2_1240k_filtered.geno\n")
cat("snp1:  v37.2_1240k_filtered.snp\n")
cat("ind1:  v37.2_1240k_filtered.ind\n")
cat(paste("geno2:\t", add1_path, "new11_1240k_pulldown_autosomes.geno\n", sep=''))
cat(paste("snp2:\t", add1_path, "new11_1240k_pulldown_autosomes.snp\n", sep=''))
cat(paste("ind2:\t", add1_path, "new11_1240k_pulldown_autosomes.ind\n", sep=''))
cat("genooutfilename: wd1/merge1.geno\n")
cat("snpoutfilename:  wd1/merge1.snp\n")
cat("indoutfilename:  wd1/merge1.ind\n")
cat("outputformat:  packedancestrymap\n")
cat("strandcheck: NO")
sink()
#` merge v37 and Iberian HG dataset
system("mergeit -p wd1/merge1.txt")

## check lost SNPs
# v37snp.df <- read.table("v37.2_1240k_filtered.snp", as.is=T)
# vm.df <- read.table(paste(add1_path, "new11_1240k_pulldown_autosomes.snp", sep=''), as.is=T)
# v37snp.df <- v37snp.df[v37snp.df$V2 != 23 & v37snp.df$V2 != 24, ]
# merge1.df <- read.table('wd1/merge1.snp', as.is=T)
# missingsnps.df <- v37snp.df[!(v37snp.df$V1 %in% merge1.df$V1), ]
# missingsnps.df <- merge(missingsnps.df, vm.df, by='V1')

#` Make par file for first merge
add2_path <- '/Users/danju/Desktop/pigment/ancientEuroVCF/OlaldeNature2019/'
#` Filter based on coverage, contamination, and kinship
#` Manually create a list of samples to exclude
olalde_exlude <- scan(paste(add2_path, 'remove.txt', sep=''), what='character')
olalde_ind.df <- read.table(paste(add2_path, 'Olalde_et_al_genotypes.ind', sep=''), as.is=T)

olalde_ind.df$V3[olalde_ind.df$V1 %in% olalde_exlude] <- 'Ignore'
write.table(olalde_ind.df, file=paste(add2_path, 'Olalde_et_al_genotypes_ignore.ind', sep=''),
            row.names = F, quote=F, col.names = F)

sink(file = "final_merge.txt")
cat("geno1: wd1/merge1.geno\n")
cat("snp1:  wd1/merge1.snp\n")
cat("ind1:  wd1/merge1.ind\n")
cat(paste("geno2:\t", add2_path, "Olalde_et_al_genotypes.geno\n", sep=''))
cat(paste("snp2:\t", add2_path, "Olalde_et_al_genotypes.snp\n", sep=''))
cat(paste("ind2:\t", add2_path, "Olalde_et_al_genotypes_ignore.ind\n", sep=''))
cat("genooutfilename: v37.2_1240k_FINAL.geno\n")
cat("snpoutfilename:  v37.2_1240k_FINAL.snp\n")
cat("indoutfilename:  v37.2_1240k_FINAL.ind\n")
cat("outputformat:  EIGENSTRAT\n")
cat("strandcheck: NO")
sink()
#` merge v37 and Olalde Iberian dataset
system("mergeit -p final_merge.txt")

## Convert EIG to VCF
## DO THIS ON PMACS HPC with eigenstrat2vcf.py
## Convert EIGENSTRAT to VCF; python 2.7.5
# python /home/danju/packages/gdc/eigenstrat2vcf.py -r v37.2_1240k_FINAL | bgzip > v37.2_1240k_FINAL.vcf.gz


## Create final annotation files with specific columns and the new datasets
#` Filter annotations for v37 dataset
v37anno.df <- anno.df[anno.df$Instance.ID %in% final.df$Instance.ID, ]
#` Retain variables for statistical analysis and rename
keep_col <- c(1,2,7,8,9,11,13,14,15,16,17,18,20,21)
v37anno.df <- v37anno.df[ , keep_col]
colnames(v37anno.df)[5] <- 'Date'
colnames(v37anno.df)[12] <- 'Ychrom'
v37_anno_save <- '/Users/danju/desktop/pigment/metadata/v37.2_anno_clean.tsv'
write.table(v37anno.df, file=v37_anno_save, row.names=F, quote=F, sep='\t')
## Manually add Villalba-MoucoCurrentBiology2019 + OlaldeNature2019 individuals
final_ind_path <- '/Users/danju/Desktop/pigment/ancientEuroVCF/v37/v37.2.1240K/final/v37.2_1240k_FINAL.ind'
final_ind.df <- read.table(final_ind_path, as.is=T)
## Remove samples from anno that were excluded in genotype files
edit_anno_path <- '/Users/danju/Desktop/pigment/metadata/v37.2_olalde_vm_anno.csv'
edit.anno.df <- read.csv(edit_anno_path, as.is=T)
edit.anno.df2 <- edit.anno.df[edit.anno.df$Instance.ID %in% final_ind.df$V1, ]
## Exclude duplicates from anno
edit.anno.df2 <- edit.anno.df2[!(duplicated(edit.anno.df2$Instance.ID)), ]
#` Save this version of annotation and manually add locations to VM and dates to
#` Olalde samples to finish the final annotation set
anno_save <- '/Users/danju/Desktop/pigment/metadata/v37.2_anno_clean2.tsv'
write.table(edit.anno.df2, file=anno_save, quote=F, row.names=F, sep='\t')
## For VM paper have to estimate longitude and latitude based on paper
#` Get rough date estimates for Olalde samples
olalde_date_path <- '/Users/danju/desktop/pigment/misc/olalde_ID_date.csv'
olalde.date.df <- read.csv(olalde_date_path, as.is=T, header=F)
olalde_dates <- rep(NA, length(olalde.date.df$V2))
for (i in 1:length(olalde.date.df$V2)) {
  date_range <- strsplit(olalde.date.df$V2[i], split = ' ')[[1]][[1]]
  try(
    {date1 <- as.numeric(strsplit(date_range, split = '–|-')[[1]][[1]])
  date2 <- as.numeric(strsplit(date_range, split = '–|-')[[1]][[2]])
  olalde_dates[i] <- (date1 + date2)/2 + 1950}
  )
}
olalde.date.df$date <- unlist(olalde_dates)
#` Fix row 191 'I10866' date manually
olalde.date.df$date[191] <- (43-51) / 2 + 1950
#' Save and use vlookup in Excel
olalde_date_save <- '/Users/danju/Desktop/pigment/misc/olalde_dates.tsv'
write.table(olalde.date.df[ ,c(1,2,4)], file=olalde_date_save, quote=F, row.names=F, sep='\t')
#` Manually remove two low coverage samples from VM (BAL003 and CHA004)

## Create anno file of samples within the past 15k years BP
anno_path <- '/Users/danju/Desktop/pigment/metadata/v37.2_anno_FINAL.csv'
anno.final.df <- read.csv(anno_path, as.is=T)
anno.15.df <- anno.final.df[anno.final.df$Date <= 15000, ]
anno_15_save <- '/Users/danju/Desktop/pigment/metadata/v37.2_anno_15kybp.tsv'
write.table(anno.15.df, file=anno_15_save, quote=F, row.names=F, sep='\t')

## Subset VCF file to include only individuals in past 15k years BP
setwd('/Users/danju/Desktop/pigment/ancientEuroVCF/v37/v37.2.1240K/final')
anno15.df <- read.csv(file='/Users/danju/Desktop/pigment/metadata/v37.2_anno_15kybp.csv', as.is=T)
subset_id <- anno15.df$Instance.ID
write.table(subset_id, file='instance.id_15kybp.txt', quote=F, row.names=F, col.names=F)
system("bcftools view -S instance.id_15kybp.txt -O z v37.2_1240k_FINAL.vcf.gz > v37.2_1240k_15kybp.vcf.gz")

## Convert VCF to Plink
system("plink --vcf v37.2_1240k_15kybp.vcf.gz --make-bed --double-id --keep-allele-order --out v37.2_1240k_15kybp")

##### Prep Human Origins dataset for PCA #####
## Create Human Origins West Eurasian genotype files
setwd('/Users/danju/Desktop/pigment/ancientEuroVCF/v37/v37.2.1240K')
#` Filter European individuals
ho_ind.df <- read.table('v37.2.1240K_HumanOrigins.ind', as.is=T)
ho_ignore <- ho_ind.df$V1[!(ho_ind.df$V1 %in% ho_samples)]
ho_ind.df$V3[ho_ind.df$V1 %in% ho_ignore] <- 'Ignore'
ho_ignore_save <- 'v37.2.1240K_HumanOrigins_ignore.ind'
write.table(ho_ind.df, file=ho_ignore_save, quote=F, row.names=F, col.names=F)

#` Make par file
sink(file = "eig2eig.txt")
cat("genotypename:  v37.2.1240K_HumanOrigins.geno\n")
cat("snpname: v37.2.1240K_HumanOrigins.snp\n")
cat("indivname: v37.2.1240K_HumanOrigins_ignore.ind\n")
cat("outputformat:  EIGENSTRAT\n")
cat("genotypeoutname: ho/eur_HO.geno\n")
cat("snpoutname:  ho/eur_HO.snp\n")
cat("indivoutname:  ho/eur_HO.ind\n")
sink()

#` Convert format from eig to PLINK .bed
system("convertf -p eig2eig.txt")

## Create individual file based on list of West Eurasian samples
ho_sample_path <- '/Users/danju/Desktop/pigment/misc/HOsamples.txt'
ho_samples <- scan(ho_sample_path, what='character')
ho_sample.df <- ho_ind.df[ho_ind.df$V1 %in% ho_samples, c(1,3)]
ho_sample_save <- '/Users/danju/Desktop/pigment/misc/HO_sample_pop.tsv'
write.table(ho_sample.df, ho_sample_save, quote=F, row.names=F, col.names=F,
            sep='\t')

ho_pop <- row.names(table(ho_sample.df$V3))
ho_pop_save <- '/Users/danju/Desktop/pigment/misc/HO_pop.tsv'
write.table(ho_pop, ho_pop_save, quote=F, row.names=F, col.names=F, sep='\t')
## runPCA_run5.txt lists files used for PCA run
