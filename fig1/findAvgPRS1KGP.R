## Find average PRS for major population groups in 1000 Genomes Project

args = commandArgs(trailingOnly=T)
score_path <- args[1]
id_path <- args[2]
out_path <- args[3]
source_path <- args[4]

source(paste(source_path, 'label1KGPpop.R', sep='/'))

score.df = read.csv(score_path, as.is=T, header=T)
id.df = read.csv(id_path, as.is=T, header=T)

# remove chromosome from name in score data frame
score.df$Sample <- lapply(score.df$ID, function(x) {
  strsplit(x, "_")[[1]][1] })
  
score.df$Sample <- unlist(score.df$Sample)
score_1kgp.df <- merge(score.df, id.df, by = 'Sample')

# assign super population labels
score_1kgp.df$superpop <- label1KGPpop(score_1kgp.df$Population)

# calculate average score per super population
super_pop <- c('EAS','EUR','AFR','AMR','SAS')
means <- aggregate(score_1kgp.df$Score, list(score_1kgp.df$superpop), mean)
stdev <- aggregate(score_1kgp.df$Score, list(score_1kgp.df$superpop), 
                   function(x) {sd(x)}
                   )
final.df <- data.frame(means$Group.1,means$x,stdev$x)
colnames(final.df) <- c('Population','Score',"St.dev")

write.table(final.df, file=out_path, row.names=F, quote=F)
