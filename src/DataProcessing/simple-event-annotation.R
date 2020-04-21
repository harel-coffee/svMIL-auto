
##Choose installation source depending on which is available...
#library(BiocInstaller)

#biocLite("StructuralVariantAnnotation", lib='/hpc/compgen/users/mnieboer/data/dataProcessing/R')

#install.packages("devtools", repos = "http://cran.us.r-project.org")

#library(devtools)
#install_github("PapenfussLab/StructuralVariantAnnotation")
#withr::with_libpaths(new = "/hpc/cog_bioinf/ridder/users/mnieboer/data/dataProcessing/R", install_github("PapenfussLab/StructuralVariantAnnotation"))
#devtools::install_github("PapenfussLab/StructuralVariantAnnotation", args = c('--library="/hpc/cog_bioinf/ridder/users/mnieboer/data/dataProcessing/R"'))

#install.packages("stringr")
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)
#' Simple SV type classifier
simpleEventType <- function(gr) {
  return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
          ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
           ifelse(strand(gr) == strand(partner(gr)), "INV",
            ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
             "DUP")))))
}

#Read through the folders and get the SV files

files <- list.files(path="/hpc/compgen/users/mnieboer/data/somatics/", pattern="*.sv.ann.vcf.gz$", full.names=TRUE, recursive=TRUE)

print(files)

#Process each file
lapply(files, function(file) {
	cat(file)
	vcf <- readVcf(file, "hg19")
	info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
	row.names=c("SIMPLE_TYPE"),
	Number=c("1"),
	Type=c("String"),
	Description=c("Simple event type annotation based purely on breakend position and orientation."))), "DataFrame"))
	gr <- breakpointRanges(vcf)
	gr <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"] #filter SVs that did not pass
	svtype <- simpleEventType(gr)
	info(vcf)$SIMPLE_TYPE <- NA_character_
	info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
	info(vcf[gr$sourceId])$SVLEN <- gr$svLen
	newName = paste0(file, '.svTypes.passed')
	writeVcf(vcf, newName)
})
