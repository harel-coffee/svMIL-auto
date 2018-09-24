files = dict(
	
	#In SV mode, snvFile can be left empty. For SNV mode, svFile can be left empty. 
	svFile = '../../data/TPTNTestSet/TP.txt', #All SVs
	snvFile = '../../data/SNVs/cosmicNCV.txt', #all SNVs
	causalGenesFile = '../../data/Genes/Census_allTue Apr 10 14_56_44 2018.tsv', #file with all causal genes from cosmic. 
	nonCausalGenesFile = '../../data/Genes/nonCosmicSubset.txt', #File with a subset of non-causal genes that are not in cosmic. 
	heatDiffusionFile = '../../data/HiC/diffusionScores.txt', # HEAT DIFFUSION IS BROKEN File with the heat diffusion scores for each region in the Hi-C network
	hiCInteractionsFile = '../../data/HiC/regions_regions_rel.csv' #All intrachromosomal regions that interact with each other
)

general = dict(
	
	mode = 'SV', #Options are: SV, SNV or SV+SNV
	tads = False, #Include TADs in the ranking yes/no, only to rank by how often TAD boundaries themselves are disrupted by SVs. 
	eQTLs = False, #Include eQTLs in the ranking yes/no
	interactionChains = False, ##THIS IS FOR HEAT DIFFUSION BUT DOES NOT WORK, SETTING SHOULD HAVE A BETTER NAME #Include genomic 3D interactions in the ranking yes/no
	gainOfInteractions = True #This depends on TADs and interactions. 
)

interactions = dict( #settings specific for Hi-C interaction data
	binSize = 5000 
	
	
)
