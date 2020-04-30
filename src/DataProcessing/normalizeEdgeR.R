data <- read.table('/hpc/compgen/users/mnieboer/data/pipeline/read_counts/pipeline_readCounts_raw.txt', header=TRUE, sep='\t')
library(edgeR)

d <- data[,2:dim(data)[2]]
rownames(d) <- data[,1]

group = rep('brca', 162)
d <- DGEList(counts = d, group = group)

TMM <- calcNormFactors(d, method = 'TMM')

TMM$samples

#plotMDS(TMM)

normCounts <- cpm(TMM)

outfile2="/hpc/compgen/users/mnieboer/data/pipeline/read_counts/brca_tmm.txt"
write.table(normCounts, file=outfile2, row.names=T,col.names=NA, quote=F,sep="\t")