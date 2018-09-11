files = dict(
	
	#In SV mode, snvFile can be left empty. For SNV mode, svFile can be left empty. 
	svFile = '../../data/TPTNTestSet/TP.txt',
	snvFile = '../../data/SNVs/cosmicNCV.txt',
	causalGenesFile = '../../data/Genes/Census_allTue Apr 10 14_56_44 2018.tsv',
	heatDiffusionFile = '../../data/HiC/diffusionScores.txt' #File with the heat diffusion scores for each region in the Hi-C network
)


general = dict(
	
	mode = 'SV', #Options are: SV, SNV or SV+SNV
	tads = False, #Include TADs in the ranking yes/no
	eQTLs = False, #Include eQTLs in the ranking yes/no
	interactions = True #Include genomic 3D interactions in the ranking yes/no
)

