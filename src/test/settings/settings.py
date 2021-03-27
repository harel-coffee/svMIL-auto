"""
	This is an example of the settings file used by svMIL2. General settings are needed
	and paths to specific files. The files need to be changed according to where the dataset
	is located, and the regulatory data should be updated based on which tissue types to use.

"""

general = dict(
	source = 'Local', #Use to tell svMIL to run on non-HMF data
	cancerType = "Test", #Set to the cancer type that you want to run on. Should match the cancer type name in the metadata file
	shuffleTads = False, #Will randomly shuffle TAD positions if True
	crd = False, #deprecated
	gtexControl = False, #deprecated
	geneENSG = False, #deprecated
	bagThreshold = 700, #If there are more than 700 bags, randomly subsample for scalability.
)
files = dict(
	svDir = "test/somatics/", #Where are the SVs per patient?
	snvDir = "test/somatics/", #Where are the SNVs per patient?
	cnvDir = "test/somatics/", #Where are the CNVs per patient?
	metadataHMF = "test/metadata/metadata.tsv", #Metadata file with patient IDs and cancer types
	expressionDir = "test/expression/", #Expression folder to check if patients have expression data or not
	normalizedExpressionFile = "test/expression/normalized_expression.txt", #Normalized expression file used by the model


	### Regulatory data, here using hmec as an example.
	#Change these paths to the data of the cancer type that you are using.
	causalGenesFile = "../data/genes/CCGC.tsv",
	nonCausalGenesFile = "../data/genes/hg19_proteinCodingGenes.bed",
	promoterFile = "../data/promoters/epdnew_hg38ToHg19_9vC8m.bed",
	cpgFile = "../data/cpg/cpgIslandExt.txt",
	tfFile = "../data/tf/tf_experimentallyValidated.bed_clustered.bed",
	chromHmmFile = "../data/chromhmm/hmec/GSE57498_HMEC_ChromHMM.bed",
	rankedGeneScoreDir = "linkedSVGenePairs",
	hg19CoordinatesFile = "../data/chromosomes/hg19Coordinates.txt",
	geneNameConversionFile = "../data/genes/allGenesAndIdsHg19.txt",
	tadFile = "../data/tads/hmec/HMEC_Lieberman-raw_TADs.bed",
	eQTLFile = "../data/eQTLs/hmec/hmec_eQTLs.bed_clustered.bed",
	enhancerFile = "../data/enhancers/hmec/hmec_encoderoadmap_elasticnet.117.txt",
	h3k4me3 = "../data/h3k4me3/hmec/ENCFF065TIH_h3k4me3.bed_clustered.bed",
	h3k27me3 = "../data/h3k27me3/hmec/ENCFF291WFP_h3k27me3.bed_clustered.bed",
	h3k27ac = "../data/h3k27ac/hmec/ENCFF154XFN_h3k27ac.bed_clustered.bed",
	h3k4me1 = "../data/h3k4me1/hmec/ENCFF336DDM_h3k4me1.bed_clustered.bed",
	dnaseIFile = "../data/dnase/hmec/ENCFF336OGZ.bed",
	rnaPolFile = "../data/rnapol/hmec/ENCFF433ZKP.bed",
	superEnhancerFile = "../data/superEnhancers/hmec/se_20200212_HMEC.bed",
	ctcfFile = "../data/ctcf/hmec/ENCFF288RFS.bed",
)
