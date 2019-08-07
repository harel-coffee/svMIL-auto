"""
	1. Shuffle the SVs, but excluding the ones that have coding effects
	2. Run the 3 non-coding based methods
	3. Compute how many of the genes are in COSMIC, and are BC genes
	4. Compute how many of the SV-gene pairs are DEG.
	5. Output a number for these 3 categories to a file corresponding to the iteration.

"""


import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import settings
from inputParser import InputParser
from neighborhoodDefiner import NeighborhoodDefiner
from geneRanking import GeneRanking
from genomicShuffler import GenomicShuffler
from scipy.stats import chi2_contingency
import pylab as plt
from matplotlib_venn import venn3, venn3_circles

permutationRound = sys.argv[7]
# 
# #1. Filter SVs based on if these cause DEGs or affect COSMIC genes in the coding way
# 
# def getSVsWithCodingEffects():
# 	
# 	codingPairs = np.loadtxt(sys.argv[1], dtype='object')
# 	degPairs = np.load(sys.argv[2], allow_pickle=True)
# 	cosmicGenesFile = sys.argv[3]
# 	
# 	codingSVGenes = dict()
# 	for pair in codingPairs:
# 		
# 		splitPair = pair.split("_")
# 		svEntries = splitPair[1:]
# 		sv = "_".join(svEntries)
# 		
# 		if sv not in codingSVGenes:
# 			codingSVGenes[sv] = []
# 		codingSVGenes[sv].append(splitPair[0])
# 	
# 		
# 	#get the COSMIC genes
# 	
# 	cosmicGenes = []
# 	with open(cosmicGenesFile, 'rb') as f:
# 		lineCount = 0
# 		for line in f:
# 			if lineCount == 0:
# 				lineCount += 1
# 				continue
# 			
# 			splitLine = line.split("\t")
# 			
# 			geneName = splitLine[0]
# 			cosmicGenes.append(geneName)
# 
# 	codingEffectSVs = []
# 	genesAffectedByFilteredSVs = []
# 	genesAffectedByCodingSVs = []
# 	for sv in codingSVGenes:
# 		
# 		degCount = 0
# 		cosmicCount = 0
# 		degAndCosmicCount = 0
# 		degOrCosmicCount = 0
# 		
# 		for gene in codingSVGenes[sv]:
# 			pair = gene + "_" + sv
# 			
# 			if gene not in genesAffectedByCodingSVs:
# 				genesAffectedByCodingSVs.append(gene)
# 			
# 			if gene in cosmicGenes or pair in degPairs[:,0]:
# 				if sv not in codingEffectSVs:
# 					codingEffectSVs.append(sv)
# 				if gene not in genesAffectedByFilteredSVs: #look at all genes that are linked to the SVs. 
# 					genesAffectedByFilteredSVs.append(gene)
# 					
# 	print "Number of genes affected in the coding way: ", len(genesAffectedByCodingSVs)
# 	print "Number of genes affected in the coding way that are DEG or COSMIC: ", len(genesAffectedByFilteredSVs)
# 	return codingEffectSVs
# 	
# codingEffectSVs = getSVsWithCodingEffects()
# print "Number of SVs filtered out with coding effects: ", len(codingEffectSVs)
# np.savetxt('codingEffectSVs.txt', codingEffectSVs, delimiter='\t', fmt='%s')
codingEffectSVs = np.loadtxt('codingEffectSVs.txt', dtype='object')
# 
# 
# #2. Find all genes within a window of the filtered SVs
# def findAffectedGenesWithinWindow():
# 	
# 	#Read all SVs and filter these for the coding effect SVs
# 	somaticSVs = InputParser().getSVsFromFile(sys.argv[4], "all", codingEffectSVs)
# 	#Shuffle the SVs
# 	genomicShuffler = GenomicShuffler()
# 	somaticSVs = genomicShuffler.shuffleSVs(somaticSVs)
# 	
# 	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
# 	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
# 	
# 	#Combine the genes into one set. 
# 	genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
# 
# 	#Specify a window to look around SVs. 2 mb seems good as I recall from quite some time ago. But this threshold may change
# 	window = 2000000
# 	
# 	affectedGenes = []
# 	svGenePairs = []
# 	#For every SV, look at 2 mb to the left of the left breakpoint, and look at 2 mb to the right of the right breakpoint.
# 	for sv in somaticSVs:
# 		
# 		#Split translocations and non-translocations
# 		if sv[0] == sv[3]:
# 
# 			#Get all genes on chr1
# 			chr1GeneSubset = genes[np.where(genes[:,0] == sv[0])]
# 			
# 			#The only genes we consider are the ones that are starting or ending within the window, and are not within the SV.
# 			#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
# 			startMatches = (chr1GeneSubset[:,1] <= sv[1]) * (chr1GeneSubset[:,2] >= (sv[1] - window)) * (chr1GeneSubset[:,2] <= sv[1])
# 			
# 			#The reverse for the end of the SV.
# 			endMatches = (chr1GeneSubset[:,2] >= sv[5]) * (chr1GeneSubset[:,1] <= (sv[5] + window)) * (chr1GeneSubset[:,1] >= sv[5])
# 			
# 			matchingGenesStart = chr1GeneSubset[startMatches] #genes must be either matching on the left or right.
# 			matchingGenesEnd = chr1GeneSubset[endMatches] #genes must be either matching on the left or right.
# 			matchingGenes = np.concatenate((matchingGenesStart, matchingGenesEnd), axis=0)
# 			for gene in matchingGenes:
# 				if gene[3].name not in affectedGenes:
# 					affectedGenes.append(gene[3].name)
# 					svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
# 				svGenePairs.append(gene[3].name + "_" + svStr)
# 	
# 		else:
# 			#look at both sides of the translocation on both chromosomes
# 			
# 			#First for chr 1
# 			chr1GeneSubset = genes[np.where(genes[:,0] == sv[0])]
# 			
# 			#The only genes we consider are the ones that are starting or ending within the window, and are not within the SV.
# 			#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
# 			startMatches = (chr1GeneSubset[:,1] <= sv[1]) * (chr1GeneSubset[:,2] >= (sv[1] - window)) * (chr1GeneSubset[:,2] <= sv[1])
# 			
# 			#The reverse for the end of the SV.
# 			endMatches = (chr1GeneSubset[:,2] >= sv[1]) * (chr1GeneSubset[:,1] <= (sv[1] + window)) * (chr1GeneSubset[:,1] >= sv[1])
# 			
# 			matchingGenesStart = chr1GeneSubset[startMatches] #genes must be either matching on the left or right.
# 			matchingGenesEnd = chr1GeneSubset[endMatches] #genes must be either matching on the left or right.
# 			matchingGenes = np.concatenate((matchingGenesStart, matchingGenesEnd), axis=0)
# 			for gene in matchingGenes:
# 				if gene[3].name not in affectedGenes:
# 					affectedGenes.append(gene[3].name)
# 					svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
# 				svGenePairs.append(gene[3].name + "_" + svStr)
# 			
# 			#Repeat for chr2
# 			chr2GeneSubset = genes[np.where(genes[:,0] == sv[3])]
# 			
# 			#The only genes we consider are the ones that are starting or ending within the window, and are not within the SV.
# 			#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
# 			startMatches = (chr2GeneSubset[:,1] <= sv[5]) * (chr2GeneSubset[:,2] >= (sv[5] - window)) * (chr2GeneSubset[:,2] <= sv[5])
# 			
# 			#The reverse for the end of the SV.
# 			endMatches = (chr2GeneSubset[:,2] >= sv[5]) * (chr2GeneSubset[:,1] <= (sv[5] + window)) * (chr2GeneSubset[:,1] >= sv[5])
# 			
# 			matchingGenesStart = chr2GeneSubset[startMatches] #genes must be either matching on the left or right.
# 			matchingGenesEnd = chr2GeneSubset[endMatches] #genes must be either matching on the left or right.
# 			matchingGenes = np.concatenate((matchingGenesStart, matchingGenesEnd), axis=0)
# 			for gene in matchingGenes:
# 				if gene[3].name not in affectedGenes:
# 					affectedGenes.append(gene[3].name)
# 					svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
# 				svGenePairs.append(gene[3].name + "_" + svStr)
# 			
# 			
# 	return affectedGenes, svGenePairs
# 	
# affectedGenesWindowed, svGenePairsWindowed = findAffectedGenesWithinWindow()
# np.savetxt("Output/windowedSVs.txt", svGenePairsWindowed, delimiter='\t', fmt='%s')
# print "Number of affected genes in windowed approach: ", len(affectedGenesWindowed)
# 
# #3. Find all genes within the TAD of the filtered SVs
# ###Alternative to look at disrupted TADs, not just boundaries
# def findAffectedGenesByTadDisruptions(codingEffectSVs):
# 		
# 	somaticSVs = InputParser().getSVsFromFile(sys.argv[4], "all", codingEffectSVs)
# 	#Shuffle the SVs
# 	genomicShuffler = GenomicShuffler()
# 	somaticSVs = genomicShuffler.shuffleSVs(somaticSVs)
# 	
# 	tads = InputParser().getTADsFromFile(sys.argv[5])
# 	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
# 	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
# 	
# 	#Combine the genes into one set. 
# 	genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
# 	
# 	#For each TAD, determine which SVs start or end within the TAD
# 	affectedGenes = []
# 	svGenePairs = []
# 	nonCodingSamples = dict()
# 	
# 	for sv in somaticSVs:
# 		svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
# 		
# 		#for non-translocations, look at all genes on the left and right of the SV, but that are within the TAD boundaries.
# 		if sv[0] == sv[3]:
# 			
# 			chr1GeneSubset = genes[genes[:,0] == sv[0]]
# 			tadChr1Subset = tads[tads[:,0] == sv[0]]
# 			
# 			#find the TAD that the left side starts in
# 			leftTAD = tadChr1Subset[(tadChr1Subset[:,1] <= sv[1]) * (tadChr1Subset[:,2] >= sv[1])]
# 			
# 			#find the TAD that the right side ends in
# 			rightTAD = tadChr1Subset[(tadChr1Subset[:,1] <= sv[5]) * (tadChr1Subset[:,2] >= sv[5])]
# 
# 			if len(leftTAD) == 0 and len(rightTAD) == 0:
# 				continue
# 
# 			leftMatchingGenes = []
# 			rightMatchingGenes = []
# 			if len(leftTAD) > 0:
# 				leftTAD = leftTAD[0]
# 				
# 				#Get all genes in these TADs that are to the left of the SV
# 				#Make sure to not include genes that are inside the SV
# 				leftEndMatchingGenes = (chr1GeneSubset[:,2] <= sv[1]) * (chr1GeneSubset[:,2] >= leftTAD[1])
# 				leftMatchingGenes = chr1GeneSubset[leftEndMatchingGenes]
# 			if len(rightTAD) > 0:
# 				rightTAD = rightTAD[0]
# 					
# 				#Get all genes in the TAD that are on the right of the SV
# 				rightStartMatchingGenes = (chr1GeneSubset[:,1] >= sv[5]) * (chr1GeneSubset[:,1] <= rightTAD[2])
# 				rightMatchingGenes = chr1GeneSubset[rightStartMatchingGenes]
# 
# 			#Get the unique genes
# 			matchingGenes = []
# 			if len(leftMatchingGenes) > 0 and len(rightMatchingGenes) > 0:
# 				matchingGenes = np.concatenate((leftMatchingGenes, rightMatchingGenes))
# 			if len(leftMatchingGenes) > 0 and len(rightMatchingGenes) == 0:
# 				matchingGenes = leftMatchingGenes
# 			if len(leftMatchingGenes) == 0 and len(rightMatchingGenes) > 0:
# 				matchingGenes = rightMatchingGenes
# 				
# 			if len(matchingGenes) > 0:
# 				
# 				for gene in matchingGenes:
# 				
# 					if gene[3].name not in affectedGenes:
# 						affectedGenes.append(gene[3].name)
# 					if gene[3].name not in nonCodingSamples:
# 						nonCodingSamples[gene[3].name] = []
# 					nonCodingSamples[gene[3].name].append(sv[7])
# 					
# 					svGenePairs.append(gene[3].name + "_" + svStr)
# 		
# 		#repeat for translocations, but look at genes on both sides of the translocation on both chromosomes
# 
# 		else:
# 			chr1GeneSubset = genes[genes[:,0] == sv[0]]
# 			tadChr1Subset = tads[tads[:,0] == sv[0]]
# 			
# 			#find the TAD that the left side starts in
# 			leftTAD = tadChr1Subset[(tadChr1Subset[:,1] <= sv[1]) * (tadChr1Subset[:,2] >= sv[1])]
# 			
# 			if len(leftTAD) >0:
# 				
# 				leftTAD = leftTAD[0]
# 				
# 				#Get all genes in these TADs that are to the left of the SV
# 				
# 				leftEndMatchingGenes = (chr1GeneSubset[:,2] <= sv[1]) * (chr1GeneSubset[:,2] >= leftTAD[1])
# 				leftMatchingGenes = chr1GeneSubset[leftEndMatchingGenes]
# 				
# 				#Get all genes in the TAD that are on the right of the SV
# 				rightStartMatchingGenes = (chr1GeneSubset[:,1] >= sv[1]) * (chr1GeneSubset[:,1] <= leftTAD[2])
# 				rightMatchingGenes = chr1GeneSubset[rightStartMatchingGenes]
# 				
# 				matchingGenes = np.concatenate((leftMatchingGenes, rightMatchingGenes))
# 				if len(matchingGenes) > 0:
# 					
# 					for gene in matchingGenes:
# 					
# 						if gene[3].name not in affectedGenes:
# 							affectedGenes.append(gene[3].name)
# 						if gene[3].name not in nonCodingSamples:
# 							nonCodingSamples[gene[3].name] = []
# 						nonCodingSamples[gene[3].name].append(sv[7])
# 						
# 						svGenePairs.append(gene[3].name + "_" + svStr)
# 				
# 			
# 			#Repeat for the TAD on chr2
# 			chr2GeneSubset = genes[genes[:,0] == sv[3]]
# 			tadChr2Subset = tads[tads[:,0] == sv[3]]
# 			
# 			#find the TAD that the left side starts in
# 			rightTAD = tadChr2Subset[(tadChr2Subset[:,1] <= sv[5]) * (tadChr2Subset[:,2] >= sv[5])]
# 			
# 			if len(rightTAD) > 0:
# 					
# 				rightTAD = rightTAD[0]
# 				
# 				#Get all genes in these TADs that are to the left of the SV
# 				leftEndMatchingGenes = (chr2GeneSubset[:,2] <= sv[5]) * (chr2GeneSubset[:,2] >= rightTAD[1])
# 				leftMatchingGenes = chr2GeneSubset[leftEndMatchingGenes]
# 				
# 				#Get all genes in the TAD that are on the right of the SV
# 				rightStartMatchingGenes = (chr2GeneSubset[:,1] >= sv[5]) * (chr2GeneSubset[:,1] <= rightTAD[2])
# 				rightMatchingGenes = chr2GeneSubset[rightStartMatchingGenes]
# 				
# 				matchingGenes = np.concatenate((leftMatchingGenes, rightMatchingGenes))
# 				if len(matchingGenes) > 0:
# 					
# 					for gene in matchingGenes:
# 					
# 						if gene[3].name not in affectedGenes:
# 							affectedGenes.append(gene[3].name)
# 						if gene[3].name not in nonCodingSamples:
# 							nonCodingSamples[gene[3].name] = []
# 						nonCodingSamples[gene[3].name].append(sv[7])
# 						
# 						svGenePairs.append(gene[3].name + "_" + svStr)
# 
# 	return affectedGenes, svGenePairs
# 
# tadAffectedGenes, tadSVGenePairs = findAffectedGenesByTadDisruptions(codingEffectSVs)
# np.savetxt("Output/tadSVs.txt", tadSVGenePairs, delimiter='\t', fmt='%s')
# 
# print "TAD affected genes: ", len(tadAffectedGenes)


#4. Run the rule-based method on the filtered SVs

#here I wil bypass main for simplicity, but could be much neater I think. the issue is passing the exlucded SVs.
def getGenesWithRuleBasedApproach():
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	causalGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
	
	svData = InputParser().getSVsFromFile(sys.argv[4], "all", codingEffectSVs)
	#Shuffle the SVs
	genomicShuffler = GenomicShuffler()
	svData = genomicShuffler.shuffleSVs(svData)
	
	NeighborhoodDefiner(causalGenes, svData, None, 'SV') #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
	
	#3. Do ranking of the genes and report the causal SVs
	print "Ranking the genes for the variants"
	geneRanking = GeneRanking(causalGenes[:,3], svData, 'SV', 'naive_shuffled', permutationRound)
	
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
ruleSvGenePairs = np.loadtxt('Output/RankedGenes/naive_shuffled/BRCA/nonCoding_geneSVPairs.txt_' + permutationRound, dtype='object')

#output the sv-gene pairs for this permutation
#np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/svGenePairsWindowed.txt_' + permutationRound, svGenePairsWindowed, delimiter='\t', fmt='%s')
#np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/tadSVGenePairs.txt_' + permutationRound, tadSVGenePairs, delimiter='\t', fmt='%s')
np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/ruleSvGenePairs.txt_' + permutationRound, ruleSvGenePairs[:,0], delimiter='\t', fmt='%s')

#Collect the number of COSMIC and BC genes

#For each set, how many of the genes are in COSMIC?
#get the COSMIC genes

cosmicGenesFile = sys.argv[3]
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

# windowedGenesCosmic = []
# for gene in affectedGenesWindowed:
# 	if gene in cosmicGenes:
# 		windowedGenesCosmic.append(gene)
# 
# tadGenesCosmic = []
# for gene in tadAffectedGenes:
# 	if gene in cosmicGenes:
# 		tadGenesCosmic.append(gene)

ruleGenesCosmic = []
for gene in ruleBasedAffectedGenes:
	if gene in cosmicGenes:
		ruleGenesCosmic.append(gene)

#Get the breast cancer specific genes
breastCancerGenesFile = sys.argv[6]
breastCancerGenes = []
with open(breastCancerGenesFile, 'r') as f:
	
	for line in f:
		
		line = line.strip()
		
		breastCancerGenes.append(line)


# windowedGenesBc = []
# for gene in affectedGenesWindowed:
# 	if gene in breastCancerGenes:
# 		windowedGenesBc.append(gene)
# 
# tadGenesBc = []
# for gene in tadAffectedGenes:
# 	if gene in breastCancerGenes:
# 		tadGenesBc.append(gene)

ruleGenesBc = []
for gene in ruleBasedAffectedGenes:
	if gene in breastCancerGenes:
		ruleGenesBc.append(gene)

print len(ruleGenesCosmic)
print ruleGenesBc
exit()


#Collect the DEG genes for each SV-gene pair combination
#The DEGs  will need to be re-computed for each shuffled iteration
windowExprCall = "python computeSVGenePairExpression_oneSet.py Output/RankedGenes/naive_shuffled/BRCA/svGenePairsWindowed.txt_" + permutationRound + " " + sys.argv[8] 
os.system(windowExprCall)
tadExprCall = "python computeSVGenePairExpression_oneSet.py Output/RankedGenes/naive_shuffled/BRCA/tadSVGenePairs.txt_" + permutationRound + " " + sys.argv[8] 
os.system(tadExprCall)
rulesExprCall = "python computeSVGenePairExpression_oneSet.py Output/RankedGenes/naive_shuffled/BRCA/ruleSvGenePairs.txt_" + permutationRound + " " + sys.argv[8] 
os.system(rulesExprCall)

#Read the DEG pairs and determine how many genes are DEG in total
#svGenePairsWindowed = np.loadtxt("Output/windowedSVs.txt", dtype='object')
windowSVsDegPairs = np.load("Output/RankedGenes/naive_shuffled/BRCA/svGenePairsWindowed.txt_" + permutationRound + "_degPairs.npy", allow_pickle=True)
tadSVsDegPairs = np.load("Output/RankedGenes/naive_shuffled/BRCA/tadSVGenePairs.txt_" + permutationRound + "_degPairs.npy", allow_pickle=True)
ruleSVsDegPairs = np.load("Output/RankedGenes/naive_shuffled/BRCA/ruleSvGenePairs.txt_" + permutationRound + "_degPairs.npy", allow_pickle=True)

print "checking for DEG overlap"
windowedDegGenes = []
for pair in svGenePairsWindowed:
	if pair in windowSVsDegPairs[:,0]:
		splitPair = pair.split("_")
		if splitPair[0] not in windowedDegGenes:
			windowedDegGenes.append(splitPair[0])

tadDegGenes = []
for pair in tadSVGenePairs:
	if pair in tadSVsDegPairs[:,0]:
		splitPair = pair.split("_")
		if splitPair[0] not in tadDegGenes:
			tadDegGenes.append(splitPair[0])

ruleDegGenes = []
for pair in ruleSvGenePairs[:,0]:
	if pair in ruleSVsDegPairs[:,0]:
		splitPair = pair.split("_")
		if splitPair[0] not in ruleDegGenes:
			ruleDegGenes.append(splitPair[0])

#Write all counts to output files

np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/Counts/windowedCosmic.txt_' + permutationRound, np.array([len(windowedGenesCosmic)]), fmt='%s')
np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/Counts/tadCosmic.txt_' + permutationRound, np.array([len(tadGenesCosmic)]), fmt='%s')
np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/Counts/rulesCosmic.txt_' + permutationRound, np.array([len(ruleGenesCosmic)]), fmt='%s')

np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/Counts/windowedBc.txt_' + permutationRound, np.array([len(windowedGenesBc)]), fmt='%s')
np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/Counts/tadBc.txt_' + permutationRound, np.array([len(tadGenesBc)]), fmt='%s')
np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/Counts/rulesBc.txt_' + permutationRound, np.array([len(ruleGenesBc)]), fmt='%s')

np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/Counts/windowedDeg.txt_' + permutationRound, np.array([len(windowedDegGenes)]), fmt='%s')
np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/Counts/tadDeg.txt_' + permutationRound, np.array([len(tadDegGenes)]), fmt='%s')
np.savetxt('Output/RankedGenes/naive_shuffled/BRCA/Counts/rulesDeg.txt_' + permutationRound, np.array([len(ruleDegGenes)]), fmt='%s')
