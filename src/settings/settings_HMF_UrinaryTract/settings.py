files = dict(

	#to run on the PCAWG data, we need both the directory of the SVs, and also the metadata to link file identifiers to cancer types
	svDir = '/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/',
	snvDir = '/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/',
	cnvDir = '/hpc/cuppen/shared_resources/HMF_data/DR-104/data/somatics/',
	metadataHMF = '/hpc/cuppen/shared_resources/HMF_data/DR-104/metadata/metadata.tsv',
	expressionDir = '/hpc/cuppen/shared_resources/HMF_data/DR-104/data/isofix/',
	normalizedExpressionFile = '../data/expression/HMF_TMM.txt',

	causalGenesFile = '../data/genes/CCGC.tsv', #file with all causal genes from cosmic.
	nonCausalGenesFile = '../data/genes/hg19_proteinCodingGenes.bed', #file with all protein-coding genes.
	promoterFile = '../data/promoters/epdnew_hg38ToHg19_9vC8m.bed', #File with promoters in human, not cell-specific
	cpgFile = '../data/cpg/cpgIslandExt.txt', #All human CpG sites
	tfFile = '../data/tf/tf_experimentallyValidated.bed_clustered.bed', #All validated human TFs
	rankedGeneScoreDir = "linkedSVGenePairs", #File that the output scores will be written to. The output will be in a folder with the provided UUID under this main results folder
	hg19CoordinatesFile = "../data/chromosomes/hg19Coordinates.txt",
	geneNameConversionFile = '../data/genes/allGenesAndIdsHg19.txt', #used with HMF SNVs and expression, converting ENSG identifiers to gene names.

	#specific for HMEC
	tadFile = "../data/tads/urinaryTract/Bladder_Donor-BL1-raw_TADs.txt",
	eQTLFile = '../data/eQTLs/urinaryTract/blood_eQTLs.bed_clustered.bed',
	enhancerFile = '../data/enhancers/urinaryTract/fantom5_elasticnet.5.csv',
	h3k4me3 = '../data/histones/urinaryTract/ENCFF599MLG_h3k4me3.bed_clustered.bed',
	h3k27ac = '../data/histones/urinaryTract/ENCFF701IRJ_h3k27ac.bed_clustered.bed',
	h3k27me3 = '../data/histones/gm12878/ENCFF247VUO_h3k27me3.bed_clustered.bed', #this one is missing, use blood
	h3k4me1 = '../data/histones/urinaryTract/ENCFF539KBP_h3k4me1.bed_clustered.bed',
	dnaseIFile = '../data/dnase/urinaryTract/ENCFF391VLP.bed',
	chromHmmFile = '../data/chromhmm/hmec/GSE57498_HMEC_ChromHMM.bed', #use hmec, missing for blood also
	rnaPolFile = '../data/rnapol/gm12878/ENCFF120VUT.bed', #same here, missing data
	superEnhancerFile = '../data/superEnhancers/urinaryTract/Bladder.bed',
	ctcfFile = '../data/ctcf/gm12878/ENCFF833FTF.bed' #same here

)

general = dict(

	source = 'HMF',
	cancerType = 'Urinary tract', #is used to get the right cancer type from the data, use the cancer type name used in the metadata.
	shuffleTads = False, #Should TAD positions be shuffled
	crd = False,
	gtexControl = False, #Should we use GTEx expression as a normal control, or use the non-disrupted TADs in other patients as control?
	geneENSG = False, #True if in the expression data the gene names are ENSG IDs. Otherwise, use False if these are regular gene names.
	bagThreshold = 700 #Subsample the bags randomly if there are more than this amount. To save computational load.
)
