files = dict(
	
	#In SV mode, snvFile can be left empty. For SNV mode, svFile can be left empty. 
	#svFile = '../../data/svs/brca_tcga_parsed_05022019.txt', #TCGA BRCA SVs
	#svFile = 'codingEffectSVs.txt',
	#svFile = '/Users/mnieboer/Documents/Projects/BoundaryMutations/data/ucec_svs_parsed_12082019.txt',
	svFile = '../../data/svs/GRCh37_hg19_variants_2020-02-25_filtered.txt', #germline SV test
	svDir = '../../../../somatics', #for running with HMF data
	snvFile = '../../data/snvs/cosmicNCV.txt', #all SNVs
	causalGenesFile = '../../data/genes/CCGC.tsv', #file with all causal genes from cosmic. 
	nonCausalGenesFile = '../../data/genes/hg19_proteinCodingGenes.bed', #file with all protein-coding genes.
	#excludedSVs = 'codingEffectSVs.txt', #move this file to the right place
	excludedSVs = 'Output/RankedGenes/12112019/BRCA/coding_geneSVPairs.txt__codingEffectSVs2.txt',
	heatDiffusionFile = '../../data/hic/diffusionScores.txt', # HEAT DIFFUSION IS BROKEN File with the heat diffusion scores for each region in the Hi-C network
	hiCInteractionsFile = '../../data/hic/regions_regions_rel.csv', #All intrachromosomal regions that interact with each other
	lncRNAFile = '../../data/lncRNAs/lncRNA.bed', #File with all lncRNAs and positions
	eQTLFile = '../../data/eQTLs/breast_eQTLs.bed', #File with eQTLs specific for breast tissue
	promoterFile = '../../data/promoters/epdnew_hg38ToHg19_9vC8m.bed', #File with promoters in human, not cell-specific
	cpgFile = '../../data/cpg/cpgIslandExt.txt', #All human CpG sites
	tfFile = '../../data/tf/tf_experimentallyValidated.bed', #All validated human TFs	
	hicFile = '../../data/hic/HMEC_groupedTADInteractions.txt',
	#methylationFile = '../../data/methylation/BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt',
	tadFile = "../../data/tads/HMEC_Lieberman-raw_TADs.bed", #File with TADs specific for breast tissue
	rankedGeneScoreDir = "Output/RankedGenes", #File that the output scores will be written to. The output will be in a folder with the provided UUID under this main results folder
	hg19CoordinatesFile = "../../data/chromosomes/hg19Coordinates.txt",
	#snvDir = '../../data/snvs/gdac.broadinstitute.org_BRCA.Mutation_Packager_Calls.Level_3.2016012800.0.0/',
	snvDir = '../../../../somatics',
	cnvDir = '../../../../somatics',
	geneNameConversionFile = '../../data/genes/allGenesAndIdsHg19.txt', #used with HMF SNVs, converting ENSG identifiers to gene names.
	
	#specific for HMEC
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
	
	#settings specific for mcf7
	# enhancerFile = '../../data/enhancers/mcf7/mcf7_fantom5_elasticnet.409.txt',
	# #enhancerFile = '../../data/enhancers/mcf7/enhancer_gene_MCF-7.txt',
	# h3k9me3 = '../../data/histones/mcf7/ENCFF454ZUB_h3k9me3.bed',
	# h3k4me3 = '../../data/histones/mcf7/ENCFF727UPU_h3k4me3.bed',
	# h3k27ac = '../../data/histones/mcf7/ENCFF223CXM_h3k27ac.bed',
	# h3k27me3 = '../../data/histones/mcf7/ENCFF151QZZ_h3k27me3.bed',
	# h3k4me1 = '../../data/histones/mcf7/ENCFF674BKS_h3k4me1.bed',
	# h3k36me3 = '../../data/histones/mcf7/ENCFF001VCW_h3k36me3.bed',
	# dnaseIFile = '../../data/dnase/mcf7/ENCFF846DFL_dnaseI.bed',
	# chromHmmFile = '../../data/chromhmm/mcf7/GSE57498_MCF7_ChromHMM.bed',
	# rnaPolFile = '../../data/rnapol/mcf7/ENCFF002DBQ.bed'
	
)

general = dict(
	
	source = 'HMF',
	mode = 'SV', #Options are: SV, SNV or SV+SNV
	cancerType = 'BRCA',
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
