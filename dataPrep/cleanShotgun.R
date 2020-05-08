## Prepare dataset of indivduals prepared by shotgun sequencing 

######## Paths ########
archaic <- '/Users/danju/Desktop/pigment/misc/v37_archaic_instance_id.txt'
# file listing archaic sample Instance ID to exclude
# folder containing eigenstrat files
v37anno_path <- '/Users/danju/Desktop/pigment/ancientEuroVCF/v37/v37.2.1240K/v37.2.1240K.clean4.anno.csv'
#shotgun_anno <- '/Users/danju/Desktop/pigment/misc/shotgun.anno.csv'
## Need to manually change the IDs in original shotgun.anno file to match the 
#` corresponding individual master.id in v37 annotation dataset. Also, remove
#` the duplicate samples in the anno file: RISE516 from Allentoft, RISE515 from
#` Allentoft, Sope (ie MA826) from Saag
shotgun_anno <- '/Users/danju/Desktop/pigment/misc/shotgun.anno_EDIT.csv'
#######################

v37anno.df <- read.csv(v37anno_path, as.is=T)
sg.anno.df <- read.csv(shotgun_anno, as.is=T)
exclude_arch <- scan(archaic, what='character')

## Remove archaic individuals from v37 anno
v37anno.df <- v37anno.df[!(v37anno.df$Instance.ID %in% exclude_arch), ]
# merge
master.df <- merge(v37anno.df, sg.anno.df, by.x='Master.ID', by.y='ID')

## Remove 1240k individuals
master.sg.df <- master.df[master.df$Data.type!='1240K' , ]

## Find which samples did not match up
sg.missing.df <- sg.anno.df[!(sg.anno.df$ID %in% master.df$Master.ID), ]

## Remove archaic individuals from new shotgun samples list
new.sg.samples.df <- sg.missing.df[!(sg.missing.df$Culture %in% c('Neanderthal','Denisovan')), ]

## Make bamlist for apulldown 1240k sites of new samples 
bam_list <- new.sg.samples.df[ , 1:2]
new_bam_save <- '/Users/danju/Desktop/pigment/misc/new_shotgun_bam1240k'
write.table(bam_list, new_bam_save, quote=F, row.names=F, col.names=F, sep='\t')
## Run apulldown.py on HPC for these non-overlapping samples

## Create 1240k VCF for overlapping shotgun samples from v37.2
ind_file_path <- '/Users/danju/Desktop/pigment/ancientEuroVCF/v37/v37.2.1240K/v37.2.1240K.ind'
v37ind.df <- read.table(ind_file_path, as.is=T)
v37ind.df$V3[!(v37ind.df$V1 %in% master.sg.df$Instance.ID)] <- 'Ignore'
ignore_sg_save <- '/Users/danju/Desktop/pigment/ancientEuroVCF/v37/v37.2.1240K/sg_ignore.ind'
write.table(v37ind.df, ignore_sg_save, quote=F, col.names=F, row.names=F)
#`Remove non-shotgun samples
setwd('/Users/danju/Desktop/pigment/ancientEuroVCF/v37/v37.2.1240K')
#` Make par file 
sink(file = "work-sg/eig2eig.txt")
cat("genotypename:  v37.2.1240K.geno\n")
cat("snpname: v37.2.1240K.snp\n")
cat("indivname: sg_ignore.ind\n")
cat("outputformat:  EIGENSTRAT\n")
cat("genotypeoutname: work-sg/sg499.geno\n")
cat("snpoutname:  work-sg/sg499.snp\n")
cat("indivoutname:  work-sg/sg499.ind\n")
sink()

#` Convert format from eig to PLINK .bed
system("convertf -p work-sg/eig2eig.txt")

## Convert EIG to VCF
## DO THIS ON PMACS HPC with eigenstrat2vcf.py; python 2.7.5
# python /home/danju/packages/gdc/eigenstrat2vcf.py -r sg499 | bgzip > sg499.vcf.gz

## Merge VCF files
setwd('/Users/danju/desktop/pigment/ancientEuroVCF/updated-sg')

system('bcftools merge -O z sg499.vcf.gz v37.2_1240k_new_ind.vcf.gz > sg_1240k_FINAL.vcf.gz')

## Remove sex chromosomes based on 1240k array
setwd('/Users/danju/Desktop/pigment/ancientEuroVCF/updated-sg')
system('bcftools view -v snps -O z -R /Users/danju/Desktop/pigment/misc/1240k_no_sex_snps.tsv sg_1240k_FINAL.vcf.gz > sg_1240k_autosome.vcf.gz')

############ Create anno file for regression analyses ############
anno_path <- '/Users/danju/Desktop/pigment/misc/shotgun.anno_EDIT.csv'
anno.df <- read.csv(anno_path, as.is=T, na.strings = '..') 

# filter by geography
filter <- which(anno.df$Latitude > 35 & anno.df$Longitude < 60 & anno.df$Longitude > -15)
anno.final.df <- anno.df[filter, ]
# remove archaics still in the dataset
anno.final.df <- anno.final.df[anno.final.df$Culture != 'Neanderthal', ]
# add Instance.ID from v37 anno downloaded off Reich Lab website
v37.anno.df <- read.csv('/Users/danju/Desktop/pigment/ancientEuroVCF/v37/v37.2.1240K/v37.2.1240K.clean4.anno.csv', as.is=T)
anno.final.df2 <- merge(anno.final.df, v37.anno.df[ , 1:2], all.x=T, by.x='ID', by.y='Master.ID')

## Fill in Instance IDs for samples with NA
anno.final.df2$Instance.ID[is.na(anno.final.df2$Instance.ID)] <- anno.final.df2$ID[is.na(anno.final.df2$Instance.ID)]

## Manually remove duplicates based on Instance.ID
rm_duplicates <- c('I1504','I5044','I4139','I5035','I1502','I5024','I4132','I5022','I5025','I5020','I5017','Kunila2','I5021','I1497','I4629','I4627','I1495','I1378','I4628','I1500','I1496','I1499','I1505','I1498','I1506','I0018','I1508','I0017','I1507','I4626','I0585','I0001','Loschbour_snpAD.DG','I0061','I4632','I5408','I1819','Kostenki14','MA968','Matojo_d','ATP2_d')

anno.final.df3 <- anno.final.df2[!(anno.final.df2$Instance.ID %in% rm_duplicates), ]

## Remove Oase1 because not shotgun sequenced
anno.final.df3 <- anno.final.df3[anno.final.df3$ID != 'Oase1', ]

## Manually remove 1st degree relatives from regression model annotation set
rm_family <- c() # no first degree relatives

# fixed date for Ukraine_N1 since that was wrong in shotgun.anno

# save file for actual regression model
out_path <- '/Users/danju/Desktop/pigment/metadata/shotgun.anno_REGRESSION.tsv'
write.table(anno.final.df3, file=out_path, quote=F, row.names=F, col.names=T, sep='\t')

## Create anno file for use in just 18 skin SNPs VCF file analysis
#` This needs to be done because the IDs for the genotype file are based on
#` Iain's original shotgun anno since all samples are pulled down together
link_path <- '/Users/danju/Desktop/pigment/misc/iainID_masterID_shotgun.txt'
link.df <- read.table(link_path, as.is=T, head=T)
anno.final.df3$ID[anno.final.df3$ID %in% link.df$Master.ID] <- link.df$iain.ID

save_anno_skin <- '/Users/danju/Desktop/pigment/metadata/sg_anno_for_manual_curated_snps.tsv'
write.table(anno.final.df3, file=save_anno_skin, quote=F, row.names=F, col.names=T, sep='\t')
