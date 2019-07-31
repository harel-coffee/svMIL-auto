"""
	Script to run all the naive & rule-based approaches, collect the output genes, and then see if there exists enrichment for COSMIC or DEG genes in these sets 

"""
	
import sys
import numpy as np
import matplotlib.pyplot as plt
import settings
from inputParser import InputParser
from neighborhoodDefiner import NeighborhoodDefiner
from geneRanking import GeneRanking

#1. Filter SVs based on if these cause DEGs or affect COSMIC genes in the coding way

def getSVsWithCodingEffects():
	
	codingPairs = np.loadtxt(sys.argv[1], dtype='object')
	degPairs = np.load(sys.argv[1] + "_nonCodingPairDEGs.npy", allow_pickle=True)
	cosmicGenesFile = sys.argv[2]
	
	codingSVGenes = dict()
	for pair in codingPairs[:,0]:
		
		splitPair = pair.split("_")
		svEntries = splitPair[1:]
		sv = "_".join(svEntries)
		
		if sv not in codingSVGenes:
			codingSVGenes[sv] = []
		codingSVGenes[sv].append(splitPair[0])
		
	#get the COSMIC genes
	
	cosmicGenes = []
	with open(cosmicGenesFile, 'rb') as f:
		lineCount = 0
		for line in f:
			if lineCount == 0:
				lineCount += 1
				continue
			
			splitLine = line.split("\t")
			
			geneName = splitLine[0]
			cosmicGenes.append(geneName)

	codingEffectSVs = []
	for sv in codingSVGenes:
		
		degCount = 0
		cosmicCount = 0
		degAndCosmicCount = 0
		degOrCosmicCount = 0
		
		for gene in codingSVGenes[sv]:
			pair = gene + "_" + sv
			
			if gene in cosmicGenes or pair in degPairs[:,0]:
				if sv not in codingEffectSVs:
					codingEffectSVs.append(sv)

	return codingEffectSVs
	
codingEffectSVs = getSVsWithCodingEffects()	

#2. Find all genes within a window of the filtered SVs
def findAffectedGenesWithinWindow():
	
	#Read all SVs and filter these for the coding effect SVs
	somaticSVs = InputParser().getSVsFromFile(sys.argv[3], "all", codingEffectSVs)
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

	#Specify a window to look around SVs. 2 mb seems good as I recall from quite some time ago. But this threshold may change
	window = 2000000
	
	affectedGenes = []
	#For every SV, look at 2 mb to the left of the left breakpoint, and look at 2 mb to the right of the right breakpoint.
	for sv in somaticSVs:
		
	
		#Get all genes on chr1
		chr1GeneSubset = genes[np.where(genes[:,0] == sv[0])]
		chr2GeneSubset = genes[np.where(genes[:,0] == sv[3])]
		
		#The only genes we consider are the ones that are starting or ending within the window, and are not within the SV.
		#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
		startMatches = (chr1GeneSubset[:,1] <= sv[1]) * (chr1GeneSubset[:,2] >= (sv[1] - window)) * (chr1GeneSubset[:,2] <= sv[1])
		
		#The reverse for the end of the SV.
		endMatches = (chr2GeneSubset[:,2] >= sv[5]) * (chr2GeneSubset[:,1] <= (sv[5] + window)) * (chr2GeneSubset[:,1] >= sv[5])
		
		matchingGenesStart = chr1GeneSubset[startMatches] #genes must be either matching on the left or right.
		matchingGenesEnd = chr2GeneSubset[endMatches] #genes must be either matching on the left or right.
		matchingGenes = np.concatenate((matchingGenesStart, matchingGenesEnd), axis=0)
		for gene in matchingGenes:
			if gene[3].name not in affectedGenes:
				affectedGenes.append(gene[3].name)
				
	
	return affectedGenes
	
#affectedGenes = findAffectedGenesWithinWindow()
#print "Number of affected genes in windowed approach: ", len(affectedGenes)

	


#3. Find all genes within the TAD of the filtered SVs
def findAffectedGenesByTadDisruptions(codingEffectSVs):
		
	somaticSVs = InputParser().getSVsFromFile(sys.argv[3], "all", codingEffectSVs)
	tads = InputParser().getTADsFromFile(sys.argv[4])
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
	
	#Remove translocations, these do not disrupt boundaries
	filteredSomaticSVs = []
	for sv in somaticSVs:
		if sv[0] == sv[3]:
			filteredSomaticSVs.append(sv)
	filteredSomaticSVs = np.array(filteredSomaticSVs, dtype="object")
	
	#For each TAD, determine which SVs start or end within the TAD, and cover the boundary
	affectedGenes = []
	nonCodingSamples = dict()
	for tad in tads:
		
		#match on chromosome
		svSubset = filteredSomaticSVs[filteredSomaticSVs[:,0] == tad[0]]
		geneChrSubset = genes[genes[:,0] == tad[0]]
		
		#if the start of the SV is within the TAD, and the end is outside of the TAD
		startingSVs = svSubset[(svSubset[:,1] >= tad[1]) * (svSubset[:,1] <= tad[2]) * (svSubset[:,5] > tad[2])]
		
		#if the end of the SV is within the TAD, and the start is outside of the TAD
		endingSVs = svSubset[(svSubset[:,5] >= tad[1]) * (svSubset[:,5] <= tad[2]) * (svSubset[:,1] < tad[1])]
		
		if len(startingSVs) > 0 or len(endingSVs) > 0:
			
			#all genes that end after the TAD start, but start before the TAD end
			matchingGenes = geneChrSubset[(geneChrSubset[:,2] >= tad[1]) * (geneChrSubset[:,1] <= tad[2])]
			
			#matchingEndGenes = (geneChrSubset[:,2] >= tad[1]) * (geneChrSubset[:,2] <= tad[2])
		
			#matchingGenes = geneChrSubset[matchingStartGenes * matchingEndGenes]
			for gene in matchingGenes:
			
				svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
				if gene[3].name not in affectedGenes:
					affectedGenes.append(gene[3].name)
				if gene[3].name not in nonCodingSamples:
					nonCodingSamples[gene[3].name] = []
				nonCodingSamples[gene[3].name].append(sv[7])
	
	return affectedGenes

tadAffectedGenes = findAffectedGenesByTadDisruptions(codingEffectSVs)
print "TAD affected genes: ", len(tadAffectedGenes)
exit()
#4. Run the rule-based method on the filtered SVs

#here I wil bypass main for simplicity, but could be much neater I think. the issue is passing the exlucded SVs.
def getGenesWithRuleBasedApproach():
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	causalGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
	
	svData = InputParser().getSVsFromFile(sys.argv[3], "all", codingEffectSVs)
	
	NeighborhoodDefiner(causalGenes, svData, None, 'SV') #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
	
	#3. Do ranking of the genes and report the causal SVs
	print "Ranking the genes for the variants"
	geneRanking = GeneRanking(causalGenes[:,3], svData, 'SV', 'none')
	
	#Read the genes from the ranking
	affectedGenes = []
	for cancerType in geneRanking.scores:

		cancerTypeScores = geneRanking.scores[cancerType]
	
		for row in range(0, cancerTypeScores.shape[0]):
			gene = cancerTypeScores[row][0]
			geneName = gene.name
			
			score = np.sum(cancerTypeScores[row][2:27]) #filter for genes that have at least 1 gain/loss
			if score > 0:
				affectedGenes.append(geneName)
	
	return affectedGenes

ruleBasedAffectedGenes = getGenesWithRuleBasedApproach()
print "rule-based affected genes: ", len(ruleBasedAffectedGenes)

#5. Compare the resulting genes between the two sets. 