knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache=FALSE)

devtools::load_all(".")

library(SummarizedExperiment)
library(dplyr)
library(DT)
library(TCGAbiolinks)

listSamples <- scan("/Users/mnieboer/Documents/Projects/StructuralVariants/Code/Outliers/sampleNames.txt", what="", sep=",")

listSamples <- c("TCGA-E9-A1NG-11A-52R-A14M-07","TCGA-BH-A1FC-11A-32R-A13Q-07",
                 "TCGA-A7-A13G-11A-51R-A13Q-07","TCGA-BH-A0DK-11A-13R-A089-07",
                 "TCGA-E9-A1RH-11A-34R-A169-07","TCGA-BH-A0AU-01A-11R-A12P-07",
                 "TCGA-C8-A1HJ-01A-11R-A13Q-07","TCGA-A7-A13D-01A-13R-A12P-07",
                 "TCGA-A2-A0CV-01A-31R-A115-07","TCGA-AQ-A0Y5-01A-11R-A14M-07")

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  experimental.strategy = "RNA-Seq",
                  platform = "Illumina HiSeq",
                  file.type = "results",
                  barcode = listSamples, 
                  legacy = TRUE)

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseqSE,"raw_count") # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")

dataGE <- dataBRCA[sample(rownames(dataBRCA),10),sample(colnames(dataBRCA),7)]

knitr::kable(dataGE[1:10,2:3], digits = 2, 
             caption = "Example of a matrix of gene expression (10 genes in rows and 2 samples in columns)",
             row.names = TRUE)

# Downstream analysis using gene expression data  
# TCGA samples from IlluminaHiSeq_RNASeqV2 with type rsem.genes.results
# save(dataBRCA, geneInfo , file = "dataGeneExpression.rda")
library(TCGAbiolinks)

# normalization of genes
dataNorm <- TCGAanalyze_Normalization(tabDF = dataBRCA, geneInfo =  geneInfo)

# quantile filter of genes
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), 
                                   typesample = c("TP"))

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])

dataDEGsFiltLevel$FDR <- format(dataDEGsFiltLevel$FDR, scientific = TRUE)
knitr::kable(dataDEGsFiltLevel[1:10,], digits = 2,
             caption = "Table of DEGs after DEA", row.names = FALSE)

pCutoff = 0.01/geneInd

geneInd = 1
significantIndices = c()
for (fdr in dataDEGsFiltLevel$FDR){
  
  if(as.numeric(fdr) < 0.01){
    fdrs = c(fdrs, fdr)
    print(fdr)
    significantIndices = c(significantIndices, geneInd)
  }
  
  geneInd = geneInd + 1
}



fdrs = dataDEGsFiltLevel$FDR
logF = dataDEGsFiltLevel$logFC
genes = dataDEGsFiltLevel$mRNA
sink("/Users/mnieboer/Documents/Projects/StructuralVariants/Code/IntegrativeOmics/data/DEG/BRCA_logFC.txt")
for (ind in 1:length(fdrs)){
  cat(genes[ind],"\t",logF[ind])  
  cat("\n")
} 

sink()






