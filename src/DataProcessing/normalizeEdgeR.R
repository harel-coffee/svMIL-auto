
totalExpression <- read.table('../data/expression/HMF_merged.txt', header=TRUE, sep='\t')
library(edgeR)

d <- totalExpression[,2:dim(totalExpression)[2]]
print(d)
#somehow row names are duplicate, add them back later
#rownames(d) <- totalExpression[,1]

group = rep('HMF', dim(d)[2])
d <- DGEList(counts = d, group = group)

TMM <- calcNormFactors(d, method = 'TMM')

normCounts <- cpm(TMM)

#add back row names as a column here
normCountsWithGeneNames <- data.frame(totalExpression[,1])
normCountsWithGeneNames <- cbind(normCountsWithGeneNames, normCounts)
colnames(normCountsWithGeneNames)[1] <- ''


outfile="/hpc/compgen/users/mnieboer/data/svMIL2/svMIL/data/expression/HMF_TMM.txt"
write.table(normCountsWithGeneNames, file=outfile, row.names=F,col.names=T, quote=F,sep="\t")