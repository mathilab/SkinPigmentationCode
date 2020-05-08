## filterByFq.R 
#` Part of plotPBSJointDiffTrees pipeline
#` Take frequency file and remove sites without enough individuals sampled

args = commandArgs(trailingOnly = T)

fq_path <- args[1]
filter_th <- as.integer(args[2])
out_file <- args[3]

fq.df <- read.csv(fq_path, as.is=T, header=T)

out.df <- fq.df[fq.df$N>=filter_th, ]

write.table(out.df, out_file, row.names=F, col.names=T, quote=F,
            sep=',')