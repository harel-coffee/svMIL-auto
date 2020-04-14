files = dict(
	
	#In SV mode, snvFile can be left empty. For SNV mode, svFile can be left empty. 
	svDir = '../../../../somatics', #for running with HMF data
	causalGenesFile = '../../data/genes/CCGC.tsv', #file with all causal genes from cosmic. 
	nonCausalGenesFile = '../../data/genes/hg19_proteinCodingGenes.bed', #file with all protein-coding genes.
	hiCInteractionsFile = '../../data/hic/regions_regions_rel.csv', #All intrachromosomal regions that interact with each other
	promoterFile = '../../data/promoters/epdnew_hg38ToHg19_9vC8m.bed', #File with promoters in human, not cell-specific
	cpgFile = '../../data/cpg/cpgIslandExt.txt', #All human CpG sites
	tfFile = '../../data/tf/tf_experimentallyValidated.bed', #All validated human TFs	
	hicFile = '../../data/hic/HMEC_groupedTADInteractions.txt',
	rankedGeneScoreDir = "linkedSVGenePairs", #File that the output scores will be written to. The output will be in a folder with the provided UUID under this main results folder
	hg19CoordinatesFile = "../../data/chromosomes/hg19Coordinates.txt",
	snvDir = '../../../../somatics',
	cnvDir = '../../../../somatics',
	geneNameConversionFile = '../../data/genes/allGenesAndIdsHg19.txt', #used with HMF SNVs, converting ENSG identifiers to gene names.
	expressionFile = '/hpc/compgen/users/mnieboer/data/pipeline/read_counts/brca_tmm.txt',
	
	#specific for HMEC
	tadFile = "../../data/tads/HMEC_Lieberman-raw_TADs.bed", #File with TADs specific for breast tissue
	eQTLFile = '../../data/eQTLs/breast_eQTLs.bed', #File with eQTLs specific for breast tissue
	enhancerFile = '../../data/enhancers/hmec_encoderoadmap_elasticnet.117.txt', #File with enhancers specific for normal cell lines
	h3k9me3 = '../../data/histones/ENCFF065FJK_H3K9me3.bed',
	h3k4me3 = '../../data/histones/ENCFF065TIH_H3K4me3.bed',
	h3k27ac = '../../data/histones/ENCFF154XFN_H3K27ac.bed',
	h3k27me3 = '../../data/histones/ENCFF291WFP_H3K27me3.bed',
	h3k4me1 = '../../data/histones/ENCFF336DDM_H3K4me1.bed',
	h3k36me3 = '../../data/histones/ENCFF906MJM_H3K36me3.bed',
	dnaseIFile = '../../data/dnase/ENCFF336OGZ_dnase_hmec.bed',
	chromHmmFile = '../../data/chromhmm/GSE57498_HMEC_ChromHMM.bed',
	rnaPolFile = '../../data/rnapol/ENCFF433ZKP.bed',
	superEnhancerFile = '../../data/enhancers/se_20200212_HMEC.bed',
	ctcfFile = '../../data/ctcf/ENCFF288RFS.bed'
	
)

general = dict(
	
	source = 'HMF',
	cancerType = 'BRCA', #used to annotate the SVs, but is not required.
	shuffleTads = False #Should TAD positions be shuffled
)

interactions = dict( #settings specific for Hi-C interaction data
	binSize = 5000 
	
	
)
