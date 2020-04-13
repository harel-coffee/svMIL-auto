files = dict(

	#In SV mode, snvFile can be left empty. For SNV mode, svFile can be left empty.
	svFile = '../../data/svs/liver_pcawg_parsed.txt',
	causalGenesFile = '../../data/genes/CCGC.tsv', #file with all causal genes from cosmic.
	nonCausalGenesFile = '../../data/genes/hg19_proteinCodingGenes.bed', #file with all protein-coding genes.
	hiCInteractionsFile = '../../data/hic/regions_regions_rel.csv', #All intrachromosomal regions that interact with each other
	lncRNAFile = '../../data/lncRNAs/lncRNA.bed', #File with all lncRNAs and positions
	eQTLFile = '../../data/eQTLs/breast_eQTLs.bed', #File with eQTLs specific for breast tissue
	promoterFile = '../../data/promoters/epdnew_hg38ToHg19_9vC8m.bed', #File with promoters in human, not cell-specific
	cpgFile = '../../data/cpg/cpgIslandExt.txt', #All human CpG sites
	tfFile = '../../data/tf/tf_experimentallyValidated.bed', #All validated human TFs
	hicFile = '../../data/hic/HMEC_groupedTADInteractions.txt',
	tadFile = "../../data/tads/HMEC_Lieberman-raw_TADs.bed", #File with TADs specific for breast tissue
	rankedGeneScoreDir = "linkedSVGenePairs", #File that the output scores will be written to. The output will be in a folder with the provided UUID under this main results folder
	hg19CoordinatesFile = "../../data/chromosomes/hg19Coordinates.txt",
	snvDir = '../../data/snvs/snv_mnv/',
	cnvFile = '../../data/cnvs/all_samples.consensus_CN.by_gene.170214.txt',
	geneNameConversionFile = '../../data/genes/allGenesAndIdsHg19.txt', #used with HMF SNVs, converting ENSG identifiers to gene names.
	metaDataFile = '../../data/svs/icgc_metadata.tsv', #metadata to map identifiers
	expressionFile = '../../data/expression/tophat_star_fpkm_uq.v2_aliquot_gl.tsv',

	#specific for liver, mixed with HMEC where data unavailable
	enhancerFile = '../../data/enhancers/liver/encoderoadmap_elasticnet.64.csv',
	h3k9me3 = '../../data/histones/liver/ENCFF845SNW_H3K9me3.bed',
	h3k4me3 = '../../data/histones/liver/ENCFF065NRN_H3K4me3.bed',
	h3k27ac = '../../data/histones/liver/ENCFF752QSK_H3K27ac.bed',
	h3k27me3 = '../../data/histones/liver/ENCFF933NEK_H3K27me3.bed',
	h3k4me1 = '../../data/histones/liver/ENCFF327BBS_H3K4me1.bed',
	h3k36me3 = '../../data/histones/liver/ENCFF744HEP_H3K36me3.bed',
	dnaseIFile = '../../data/dnase/liver/ENCFF286LYP.bed',
	chromHmmFile = '../../data/chromhmm/GSE57498_HMEC_ChromHMM.bed', #HMEC
	rnaPolFile = '../../data/rnapol/ENCFF433ZKP.bed', #HMEC
	superEnhancerFile = '../../data/enhancers/se_20200411_hepg2.bed',
	ctcfFile = '../../data/ctcf/liver/ENCFF690BYG.bed'

)

general = dict(

	source = 'PCAWG',
	mode = 'SV', #Options are: SV, SNV or SV+SNV
	cancerType = 'LIVER',
	nonCoding = True, #Should we count gains/losses for genes that are affected by non-coding SVS?
	coding = False, #Should we include gains and losses of non-coding elements, or only genes that are directly affected by SVs?
	snvs = False, #Should we include or exclude genes that have SNVs?
	cnvs = False, #should we include or exclude genes that have CNVs?
	gains = True, #Enable looking at gains of features
	losses = True, #enable looking at losses of features
	tads = True, #Include TADs in the ranking yes/no, only to rank by how often TAD boundaries themselves are disrupted by SVs. (we are getting very dependent on TADs, so perhaps force this as always required)
	eQTLs = True, #Include eQTLs in the ranking yes/no
	enhancers = True, #Include enhancers in the ranking yes/no
	promoters = True, #Include promoters in the ranking yes/no
	cpgIslands = True, #Include CpG islands in the ranking yes/no
	transcriptionFactors = True, #Include TFs in the ranking yes/no
	hiC = True, #Include HiC interactions in the ranking yes/no
	histones = True, #Include the histone marks in the ranking yes/no (maybe these need to be split later)
	dnaseI = True, #Include DNAse I hypersensitivity sites yes/no
	chromHMM = True,
	rnaPol = True,
	superEnhancers = True,
	ctcfSites = True,
	boundaryStrength = True,
	methylation = False,
	gainOfInteractions = True, #This depends on TADs and interactions.
	shuffleTads = False, #Should TAD positions be shuffled
	lncRNA = False #use lncRNAs instead of eQTLs
)

interactions = dict( #settings specific for Hi-C interaction data
	binSize = 5000


)
