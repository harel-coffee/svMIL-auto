
data <- read.table('../data/expression/brcaGtex_merged.txt', header=TRUE, sep='\t')
library(edgeR)

d <- data[,2:dim(data)[2]]
rownames(d) <- data[,1]

group <- as.factor(c(rep('GTEx', 290), rep('BRCA', 169)))

dge <- DGEList(counts = d, group = group)

TMM <- calcNormFactors(dge, method = 'TMM')

TMM$samples


normCounts <- cpm(TMM)



#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

#BiocManager::install("RUVSeq")

###how were the data merged?

data <- read.table('/hpc/compgen/users/mnieboer/data/evidenceStacking/data/expression/brcaGtex_merged.txt', header=TRUE, sep='\t')
d <- data[,2:dim(data)[2]]
rownames(d) <- data[,1]

housekeepingGenes <- read.table('../data/genes/HK_genes.txt', header=FALSE, sep='\t')
housekeepingGenes <- housekeepingGenes[,1]
housekeepingGenes

housekeepingGenes <- gsub(" ", "", housekeepingGenes, fixed=TRUE)

#housekeepingGenes = c('C1orf43', 'CHMP2A', 'EMC7', 'GPI', 'PSMB2', 'PSMB4', 'RAB7A', 'REEP5', 'SNRPD3','VCP', 'VPS29')


#genes <- rownames(d)[grep(paste(appendedHKGenes,collapse="|"), rownames(d))]
genes <- housekeepingGenes

genes = c()
for (rowname in rownames(d)){

	if (rowname %in% housekeepingGenes){
		genes = c(genes, rowname)
	}

}
genes
library(RUVSeq)


x <- as.factor(c(rep('GTEx', 290), rep('BRCA', 169)))

normCounts <- normCounts * 2
normCounts <- ceiling(normCounts)

set <- newSeqExpressionSet(as.matrix(normCounts),phenoData = data.frame(x, row.names=colnames(normCounts)))
#set <- newSeqExpressionSet(as.matrix(d),phenoData = data.frame(x, row.names=colnames(d)))

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
#plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#plotPCA(set, col=colors[x], cex=1.2)

#upper-quartile normalization
#set <- betweenLaneNormalization(set, which="upper")
#set <- betweenLaneNormalization(set, which="full")
#plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#plotPCA(set, col=colors[x], cex=1.2)


set1 <- RUVg(set, genes, k=1)

#plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#plotPCA(set1, col=colors[x], cex=1.2)


#output the normalized counts
normalizedCounts = normCounts(set1)
#normalizedCounts = counts(set1)

#normalizedCounts

outfile2='../data/expression/brcaGtex_ruv_TMM.txt'
write.table(normalizedCounts, file=outfile2, row.names=T,col.names=NA, quote=F,sep="\t")

#split the counts into BRCA and GTEx files and output these as well
gtexCounts <- normalizedCounts[,1:290]

colnames(gtexCounts)
dim(gtexCounts)

brcaCounts <- normalizedCounts[,291:459]
colnames(brcaCounts)
dim(brcaCounts)

outfile3='../data/expression/brca_ruv_TMM.txt'
write.table(brcaCounts, file=outfile3, row.names=T,col.names=NA, quote=F,sep="\t")

outfile4='../data/expression/gtex_ruv_TMM.txt'
write.table(gtexCounts, file=outfile4, row.names=T,col.names=NA, quote=F,sep="\t")