files = dict(

	svDir = '../../../somatics', #for running with HMF data
	causalGenesFile = '../data/genes/CCGC.tsv', #file with all causal genes from cosmic.
	nonCausalGenesFile = '../data/genes/hg19_proteinCodingGenes.bed', #file with all protein-coding genes.
	promoterFile = '../data/promoters/epdnew_hg38ToHg19_9vC8m.bed', #File with promoters in human, not cell-specific
	cpgFile = '../data/cpg/cpgIslandExt.txt', #All human CpG sites
	tfFile = '../data/tf/tf_experimentallyValidated.bed_clustered.bed', #All validated human TFs
	rankedGeneScoreDir = "linkedSVGenePairs", #File that the output scores will be written to. The output will be in a folder with the provided UUID under this main results folder
	hg19CoordinatesFile = "../data/chromosomes/hg19Coordinates.txt",
	snvDir = '../../../somatics',
	cnvDir = '../../../somatics',
	geneNameConversionFile = '../data/genes/allGenesAndIdsHg19.txt', #used with HMF SNVs and expression, converting ENSG identifiers to gene names.
	expressionFile = '../data/expression/brca_ruv_TMM.txt',
	#expressionFile = '/hpc/compgen/users/mnieboer/data/pipeline/read_counts/brca_tmm.txt',
	gtexExpressionFile = '../data/expression/gtex_ruv_TMM.txt',

	#specific for HMEC
	tadFile = "../data/tads/hg19.TADs/GM12878_Lieberman-raw_TADs.txt", #File with TADs specific for breast tissue
	eQTLFile = '../data/eQTLs/gm12878/blood_eQTLs.bed_clustered.bed', #File with eQTLs specific for breast tissue
	enhancerFile = '../data/enhancers/gm12878/encoderoadmap_elasticnet.114.csv', #File with enhancers specific for normal cell lines
	h3k9me3 = '../data/histones/hmec/ENCFF065FJK_H3K9me3.bed',
	h3k4me3 = '../data/histones/gm12878/ENCFF295GNH_h3k4me3.bed_clustered.bed',
	h3k27ac = '../data/histones/gm12878/ENCFF816AHV_h3k27ac.bed_clustered.bed',
	h3k27me3 = '../data/histones/gm12878/ENCFF247VUO_h3k27me3.bed_clustered.bed',
	h3k4me1 = '../data/histones/gm12878/ENCFF921LKB_h3k4me1.bed_clustered.bed',
	h3k36me3 = '../data/histones/hmec/ENCFF906MJM_H3K36me3.bed',
	dnaseIFile = '../data/dnase/gm12878/ENCFF097LEF.bed_clustered.bed',
	chromHmmFile = '../data/chromhmm/hmec/GSE57498_HMEC_ChromHMM.bed',
	rnaPolFile = '../data/rnapol/gm12878/ENCFF120VUT.bed_clustered.bed',
	superEnhancerFile = '../data/superEnhancers/gm12878/GM12878.bed',
	ctcfFile = '../data/ctcf/gm12878/ENCFF833FTF.bed_clustered.bed',
	hicFile = '../data/hic/HMEC_groupedTADInteractions.txt'

)

general = dict(

	source = 'HMF',
	cancerType = 'BRCA', #used to annotate the SVs, but is not required.
	shuffleTads = False, #Should TAD positions be shuffled
	crd = False,
	gtexControl = False, #Should we use GTEx expression as a normal control, or use the non-disrupted TADs in other patients as control?
	geneENSG = False, #True if in the expression data (e.g. from HMF) the gene names are ENSG IDs. Otherwise, use False if these are regular gene names.
	bagThreshold = 700 #Subsample the bags randomly if there are more than this amount. To save computational load.
)
