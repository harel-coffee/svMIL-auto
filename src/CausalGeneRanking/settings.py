files = dict(
	
	#In SV mode, snvFile can be left empty. For SNV mode, svFile can be left empty. 
	svFile = '../../data/TPTNTestSet/brca_tcga_parsed_05022019.txt', #TCGA BRCA SVs
	snvFile = '../../data/SNVs/cosmicNCV.txt', #all SNVs
	causalGenesFile = '../../data/Genes/CCGC.tsv', #file with all causal genes from cosmic. 
	nonCausalGenesFile = '../../data/Genes/hg19_proteinCodingGenes.bed', #file with all protein-coding genes. 
	heatDiffusionFile = '../../data/HiC/diffusionScores.txt', # HEAT DIFFUSION IS BROKEN File with the heat diffusion scores for each region in the Hi-C network
	hiCInteractionsFile = '../../data/HiC/regions_regions_rel.csv', #All intrachromosomal regions that interact with each other
	lncRNAFile = '../../data/lncRNAs/lncRNA.bed', #File with all lncRNAs and positions
	eQTLFile = '../../data/eQTLs/breast_eQTLs.txt', #File with eQTLs specific for breast tissue
	enhancerFile = '../../data/enhancers/enhancer_gene_GM12878.txt', #File with enhancers specific for normal cell lines
	promoterFile = '../../data/promoters/epdnew_hg38ToHg19_9vC8m.bed', #File with promoters in human, not cell-specific
	cpgFile = '../../data/cpg/cpgIslandExt.txt', #All human CpG sites
	tfFile = '../../data/tf/allTFs.bed', #All human TFs from gm12878
	hicFile = '../../data/HiC/HMEC_interactions.txt',
	tadFile = "../../data/tads/HMEC_Lieberman-raw_TADs.bed", #File with TADs specific for breast tissue
	rankedGeneScoreDir = "./RankedGenes" #File that the output scores will be written to. The output will be in a folder with the provided UUID under this main results folder
)

general = dict(
	
	mode = 'SV', #Options are: SV, SNV or SV+SNV
	cancerType = 'BRCA', #Use to specify which cancer type the data should be filtered by
	tads = True, #Include TADs in the ranking yes/no, only to rank by how often TAD boundaries themselves are disrupted by SVs. (we are getting very dependent on TADs, so perhaps force this as always required)
	eQTLs = True, #Include eQTLs in the ranking yes/no
	enhancers = True, #Include enhancers in the ranking yes/no
	promoters = True, #Include promoters in the ranking yes/no
	cpgIslands = True, #Include CpG islands in the ranking yes/no
	transcriptionFactors = False, #Include TFs in the ranking yes/no
	hiC = False, #Include HiC interactions in the ranking yes/no
	gainOfInteractions = True, #This depends on TADs and interactions.
	shuffleTads = False, #Should TAD positions be shuffled
	lncRNA = False #use lncRNAs instead of eQTLs
)

interactions = dict( #settings specific for Hi-C interaction data
	binSize = 5000 
	
	
)
