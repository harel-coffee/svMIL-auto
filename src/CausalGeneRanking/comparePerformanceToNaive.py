"""
	Script to run all the naive & rule-based approaches, collect the output genes, and then see if there exists enrichment for COSMIC or DEG genes in these sets 

"""
	
from __future__ import absolute_import
from __future__ import print_function
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import settings
from inputParser import InputParser
from neighborhoodDefiner import NeighborhoodDefiner
from geneRanking import GeneRanking
from scipy.stats import chi2_contingency
import pylab as plt
from matplotlib_venn import venn3, venn3_circles
from scipy import stats

#Get the coding DEG pairs, these can be used later to fitler out SV-gene pairs that are DEG due to nearby coding
codingDegPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')
splitCodingDegPairs = dict()
for pair in codingDegPairs[:,0]:
	splitPair = pair.split("_")
	newPair = splitPair[0] + "_" + splitPair[len(splitPair)-1]
	splitCodingDegPairs[newPair] = 0
# 
# codingDegPairs = np.loadtxt(sys.argv[1], dtype='object')
# splitCodingDegPairs = dict()
# for pair in codingDegPairs:
# 	splitPair = pair.split("_")
# 	newPair = splitPair[0] + "_" + splitPair[len(splitPair)-1]
# 	splitCodingDegPairs[newPair] = 0

# 1. Filter SVs based on if these cause DEGs or affect COSMIC genes in the coding way
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
# 	codingEffectPairs = []
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
# 				if pair not in codingEffectPairs:
# 					codingEffectPairs.append(pair)
# 					
# 	#print "Number of genes affected in the coding way: ", len(genesAffectedByCodingSVs)
# 	#print "Number of genes affected in the coding way that are DEG or COSMIC: ", len(genesAffectedByFilteredSVs)
# 	return codingEffectSVs, codingEffectPairs
# 	
# codingEffectSVs, codingEffectPairs = getSVsWithCodingEffects()
# #print "Number of SVs filtered out with coding effects: ", len(codingEffectSVs)
# np.savetxt('codingEffectSVs.txt', codingEffectSVs, delimiter='\t', fmt='%s')
# np.savetxt('codingEffectPairs.txt', codingEffectPairs, delimiter='\t', fmt='%s')
codingEffectSVs = np.loadtxt('codingEffectSVs.txt', dtype='object')

#2. Find all genes within a window of the filtered SVs
def findAffectedGenesWithinWindow():
	
	#Read all SVs and filter these for the coding effect SVs
	somaticSVs = InputParser().getSVsFromFile(sys.argv[4], "all", [])
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	
	#Fiter out coding effect SVs
	svData = []
	for sv in somaticSVs:
		svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
		if svStr not in codingEffectSVs:
			svData.append(sv)
	somaticSVs = np.array(svData, dtype='object')
	
	#Combine the genes into one set. 
	genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

	#Specify a window to look around SVs. 2 mb seems good as I recall from quite some time ago. But this threshold may change
	window = 2000000
	
	affectedGenes = []
	svGenePairs = []
	#For every SV, look at 2 mb to the left of the left breakpoint, and look at 2 mb to the right of the right breakpoint.
	for sv in somaticSVs:
		svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
		
		#Split translocations and non-translocations
		if sv[0] == sv[3]:

			#Get all genes on chr1
			chr1GeneSubset = genes[np.where(genes[:,0] == sv[0])]
			
			#The only genes we consider are the ones that are starting or ending within the window, and are not within the SV.
			#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
			startMatches = (chr1GeneSubset[:,1] <= sv[1]) * (chr1GeneSubset[:,2] >= (sv[1] - window)) * (chr1GeneSubset[:,2] <= sv[1])
			
			#The reverse for the end of the SV.
			endMatches = (chr1GeneSubset[:,2] >= sv[5]) * (chr1GeneSubset[:,1] <= (sv[5] + window)) * (chr1GeneSubset[:,1] >= sv[5])
			
			matchingGenesStart = chr1GeneSubset[startMatches] #genes must be either matching on the left or right.
			matchingGenesEnd = chr1GeneSubset[endMatches] #genes must be either matching on the left or right.
			matchingGenes = np.concatenate((matchingGenesStart, matchingGenesEnd), axis=0)
			for gene in matchingGenes:
				if gene[3].name not in affectedGenes:
					affectedGenes.append(gene[3].name)
				
				newPair = gene[3].name + "_" + sv[7]
				if newPair not in splitCodingDegPairs:
					svGenePairs.append(gene[3].name + "_" + svStr)
	
		else:
			#look at both sides of the translocation on both chromosomes
			
			#First for chr 1
			chr1GeneSubset = genes[np.where(genes[:,0] == sv[0])]
			
			#The only genes we consider are the ones that are starting or ending within the window, and are not within the SV.
			#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
			startMatches = (chr1GeneSubset[:,1] <= sv[1]) * (chr1GeneSubset[:,2] >= (sv[1] - window)) * (chr1GeneSubset[:,2] <= sv[1])
			
			#The reverse for the end of the SV.
			endMatches = (chr1GeneSubset[:,2] >= sv[1]) * (chr1GeneSubset[:,1] <= (sv[1] + window)) * (chr1GeneSubset[:,1] >= sv[1])
			
			matchingGenesStart = chr1GeneSubset[startMatches] #genes must be either matching on the left or right.
			matchingGenesEnd = chr1GeneSubset[endMatches] #genes must be either matching on the left or right.
			matchingGenes = np.concatenate((matchingGenesStart, matchingGenesEnd), axis=0)
			
		
			for gene in matchingGenes:
				if gene[3].name not in affectedGenes:
					affectedGenes.append(gene[3].name)
					
				newPair = gene[3].name + "_" + sv[7]
				
				if newPair not in splitCodingDegPairs:
					svGenePairs.append(gene[3].name + "_" + svStr)
				
				
			#Repeat for chr2
			chr2GeneSubset = genes[np.where(genes[:,0] == sv[3])]
			
			#The only genes we consider are the ones that are starting or ending within the window, and are not within the SV.
			#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
			startMatches = (chr2GeneSubset[:,1] <= sv[5]) * (chr2GeneSubset[:,2] >= (sv[5] - window)) * (chr2GeneSubset[:,2] <= sv[5])
			
			#The reverse for the end of the SV.
			endMatches = (chr2GeneSubset[:,2] >= sv[5]) * (chr2GeneSubset[:,1] <= (sv[5] + window)) * (chr2GeneSubset[:,1] >= sv[5])
			
			matchingGenesStart = chr2GeneSubset[startMatches] #genes must be either matching on the left or right.
			matchingGenesEnd = chr2GeneSubset[endMatches] #genes must be either matching on the left or right.
			matchingGenes = np.concatenate((matchingGenesStart, matchingGenesEnd), axis=0)
			
			
			for gene in matchingGenes:
				if gene[3].name not in affectedGenes:
					affectedGenes.append(gene[3].name)
				newPair = gene[3].name + "_" + sv[7]
				
				if newPair not in splitCodingDegPairs:
				
					svGenePairs.append(gene[3].name + "_" + svStr)
			
	return affectedGenes, svGenePairs
	
affectedGenesWindowed, svGenePairsWindowed = findAffectedGenesWithinWindow()
np.savetxt('svGenePairsWindowed.txt', svGenePairsWindowed, delimiter='\t', fmt='%s')
print("Number of affected genes in windowed approach: ", len(affectedGenesWindowed))
print("Number of SV-gene pairs in windowed approach: ", len(svGenePairsWindowed))


###Alternative to look at disrupted TADs, not just boundaries
def findAffectedGenesByTadDisruptions(codingEffectSVs):
		
	somaticSVs = InputParser().getSVsFromFile(sys.argv[4], "all", [])
	tads = InputParser().getTADsFromFile(sys.argv[5])
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
	
	#Fiter out coding effect SVs
	svData = []
	for sv in somaticSVs:
		svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
		if svStr not in codingEffectSVs:
			svData.append(sv)
	somaticSVs = np.array(svData, dtype='object')
	
	#For each TAD, determine which SVs start or end within the TAD
	affectedGenes = []
	svGenePairs = []
	nonCodingSamples = dict()
	
	for sv in somaticSVs:
		svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
		
		#for non-translocations, look at all genes on the left and right of the SV, but that are within the TAD boundaries.
		if sv[0] == sv[3]:
			
			chr1GeneSubset = genes[genes[:,0] == sv[0]]
			tadChr1Subset = tads[tads[:,0] == sv[0]]
			
			#find the TAD that the left side starts in
			leftTAD = tadChr1Subset[(tadChr1Subset[:,1] <= sv[1]) * (tadChr1Subset[:,2] >= sv[1])]
			
			#find the TAD that the right side ends in
			rightTAD = tadChr1Subset[(tadChr1Subset[:,1] <= sv[5]) * (tadChr1Subset[:,2] >= sv[5])]

			if len(leftTAD) == 0 and len(rightTAD) == 0:
				continue

			leftMatchingGenes = []
			rightMatchingGenes = []
			if len(leftTAD) > 0:
				leftTAD = leftTAD[0]
				
				#Get all genes in these TADs that are to the left of the SV
				#Make sure to not include genes that are inside the SV
				leftEndMatchingGenes = (chr1GeneSubset[:,2] <= sv[1]) * (chr1GeneSubset[:,2] >= leftTAD[1])
				leftMatchingGenes = chr1GeneSubset[leftEndMatchingGenes]
			if len(rightTAD) > 0:
				rightTAD = rightTAD[0]
					
				#Get all genes in the TAD that are on the right of the SV
				rightStartMatchingGenes = (chr1GeneSubset[:,1] >= sv[5]) * (chr1GeneSubset[:,1] <= rightTAD[2])
				rightMatchingGenes = chr1GeneSubset[rightStartMatchingGenes]

			#Get the unique genes
			if len(leftMatchingGenes) > 0 and len(rightMatchingGenes) > 0:
				matchingGenes = np.concatenate((leftMatchingGenes, rightMatchingGenes))
			if len(leftMatchingGenes) > 0 and len(rightMatchingGenes) == 0:
				matchingGenes = leftMatchingGenes
			if len(leftMatchingGenes) == 0 and len(rightMatchingGenes) > 0:
				matchingGenes = rightMatchingGenes
				
			if len(matchingGenes) > 0:
				
				for gene in matchingGenes:
				
					if gene[3].name not in affectedGenes:
						affectedGenes.append(gene[3].name)
					if gene[3].name not in nonCodingSamples:
						nonCodingSamples[gene[3].name] = []
					nonCodingSamples[gene[3].name].append(sv[7])
					
					newPair = gene[3].name + "_" + sv[7]
					if newPair not in splitCodingDegPairs:
						svGenePairs.append(gene[3].name + "_" + svStr)
		
		#repeat for translocations, but look at genes on both sides of the translocation on both chromosomes

		else:
			chr1GeneSubset = genes[genes[:,0] == sv[0]]
			tadChr1Subset = tads[tads[:,0] == sv[0]]
			
			#find the TAD that the left side starts in
			leftTAD = tadChr1Subset[(tadChr1Subset[:,1] <= sv[1]) * (tadChr1Subset[:,2] >= sv[1])]
			
			if len(leftTAD) >0:
				
				leftTAD = leftTAD[0]
				
				#Get all genes in these TADs that are to the left of the SV
				
				leftEndMatchingGenes = (chr1GeneSubset[:,2] <= sv[1]) * (chr1GeneSubset[:,2] >= leftTAD[1])
				leftMatchingGenes = chr1GeneSubset[leftEndMatchingGenes]
				
				#Get all genes in the TAD that are on the right of the SV
				rightStartMatchingGenes = (chr1GeneSubset[:,1] >= sv[1]) * (chr1GeneSubset[:,1] <= leftTAD[2])
				rightMatchingGenes = chr1GeneSubset[rightStartMatchingGenes]
				
				matchingGenes = np.concatenate((leftMatchingGenes, rightMatchingGenes))
				if len(matchingGenes) > 0:
					
					for gene in matchingGenes:
					
						if gene[3].name not in affectedGenes:
							affectedGenes.append(gene[3].name)
						if gene[3].name not in nonCodingSamples:
							nonCodingSamples[gene[3].name] = []
						nonCodingSamples[gene[3].name].append(sv[7])
						
						newPair = gene[3].name + "_" + sv[7]
						if newPair not in splitCodingDegPairs:
							svGenePairs.append(gene[3].name + "_" + svStr)
				
			
			#Repeat for the TAD on chr2
			chr2GeneSubset = genes[genes[:,0] == sv[3]]
			tadChr2Subset = tads[tads[:,0] == sv[3]]
			
			#find the TAD that the left side starts in
			rightTAD = tadChr2Subset[(tadChr2Subset[:,1] <= sv[5]) * (tadChr2Subset[:,2] >= sv[5])]
			
			if len(rightTAD) > 0:
					
				rightTAD = rightTAD[0]
				
				#Get all genes in these TADs that are to the left of the SV
				leftEndMatchingGenes = (chr2GeneSubset[:,2] <= sv[5]) * (chr2GeneSubset[:,2] >= rightTAD[1])
				leftMatchingGenes = chr2GeneSubset[leftEndMatchingGenes]
				
				#Get all genes in the TAD that are on the right of the SV
				rightStartMatchingGenes = (chr2GeneSubset[:,1] >= sv[5]) * (chr2GeneSubset[:,1] <= rightTAD[2])
				rightMatchingGenes = chr2GeneSubset[rightStartMatchingGenes]
				
				matchingGenes = np.concatenate((leftMatchingGenes, rightMatchingGenes))
				if len(matchingGenes) > 0:
					
					for gene in matchingGenes:
					
						if gene[3].name not in affectedGenes:
							affectedGenes.append(gene[3].name)
						if gene[3].name not in nonCodingSamples:
							nonCodingSamples[gene[3].name] = []
						nonCodingSamples[gene[3].name].append(sv[7])
						
						newPair = gene[3].name + "_" + sv[7]
						if newPair not in splitCodingDegPairs:
							svGenePairs.append(gene[3].name + "_" + svStr)

	return affectedGenes, svGenePairs

tadAffectedGenes, tadSVGenePairs = findAffectedGenesByTadDisruptions(codingEffectSVs)
np.savetxt('tadSVGenePairs.txt', tadSVGenePairs, delimiter='\t', fmt='%s')

print("TAD affected genes: ", len(tadAffectedGenes))
print("Number of SV-gene pairs TAD-based: ", len(tadSVGenePairs))

#4. Run the rule-based method on the filtered SVs

#here I wil bypass main for simplicity, but could be much neater I think. the issue is passing the exlucded SVs.
def getGenesWithRuleBasedApproach():
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	causalGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
	
	svData = InputParser().getSVsFromFile(sys.argv[4], "all", [])
	
	NeighborhoodDefiner(causalGenes, svData, None, 'SV', codingEffectSVs) #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
	
	# somaticSVs = []
	# for sv in svData:
	# 	svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
	# 	if svStr not in codingEffectSVs:
	# 		somaticSVs.append(sv)
	# somaticSVs = np.array(somaticSVs, dtype='object')
	# 
	#3. Do ranking of the genes and report the causal SVs
	print("Ranking the genes for the variants")
	geneRanking = GeneRanking(causalGenes[:,3], svData, 'SV', 'naive', 'none')
	
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
# 
# ruleBasedAffectedGenes = getGenesWithRuleBasedApproach()
# print("rule-based affected genes: ", len(ruleBasedAffectedGenes))
# ruleSvGenePairs = np.loadtxt('Output/RankedGenes/naive/BRCA//nonCoding_geneSVPairs.txt_none', dtype='object')
# 
# np.savetxt('ruleSvGenePairs_withFeatures.txt', ruleSvGenePairs, delimiter='\t', fmt='%s')
# np.savetxt('ruleSvGenePairs.txt', ruleSvGenePairs[:,0], delimiter='\t', fmt='%s')
# 
# #Save all genes in memory to prevent re-computing every time
# np.savetxt('affectedGenesWindowed.txt', affectedGenesWindowed, delimiter='\t', fmt='%s')
# np.savetxt('tadAffectedGenes.txt', tadAffectedGenes, delimiter='\t', fmt='%s')
# np.savetxt('ruleBasedAffectedGenes.txt', ruleBasedAffectedGenes, delimiter='\t', fmt='%s')
# 
# print("Number of sv-gene pairs windowed: ", len(svGenePairsWindowed))
# print("Number of sv-gene pairs tads: ", len(tadSVGenePairs))
# print("Number of sv-gene pairs rules: ", ruleSvGenePairs.shape)



affectedGenesWindowed = np.loadtxt('affectedGenesWindowed.txt', dtype='object')
tadAffectedGenes = np.loadtxt('tadAffectedGenes.txt', dtype='object')
ruleBasedAffectedGenes = np.loadtxt('ruleBasedAffectedGenes.txt', dtype='object')

svGenePairsWindowed = np.loadtxt('svGenePairsWindowed.txt', dtype='object')
tadSVGenePairs = np.loadtxt('tadSVGenePairs.txt', dtype='object')
ruleSvGenePairs = np.loadtxt('ruleSvGenePairs.txt', dtype='object')

#5. Compare the resulting genes between the approaches

#Which genes are different between the TAD and rule based approach?
print("rule-based unique genes not in TAD:")
diffGenes = np.setdiff1d(ruleBasedAffectedGenes, tadAffectedGenes)
for gene in diffGenes:
	print(gene)
print("TAD unique genes not in window:")
diffGenes = np.setdiff1d(tadAffectedGenes, affectedGenesWindowed)
for gene in diffGenes:
	print(gene)
print("rule unique genes not in window:")
diffGenes = np.setdiff1d(ruleBasedAffectedGenes, affectedGenesWindowed)
for gene in diffGenes:
	print(gene)

#Compute the DEGs for each pair set
windowExprCall = "python computeSVGenePairExpression_oneSet.py svGenePairsWindowed.txt" + " " + sys.argv[2] + " " + sys.argv[7] + ' False'
os.system(windowExprCall)
tadExprCall = "python computeSVGenePairExpression_oneSet.py tadSVGenePairs.txt" + " " + sys.argv[2] + " " + sys.argv[7] + ' False'
os.system(tadExprCall)
rulesExprCall = "python computeSVGenePairExpression_oneSet.py ruleSvGenePairs.txt" + " " + sys.argv[2] + " " + sys.argv[7] + ' False'
os.system(rulesExprCall)

# Read the DEG pairs and determine how many genes are DEG in total
windowSVsDegPairs = np.load("svGenePairsWindowed.txt_degPairs.npy", allow_pickle=True, encoding='latin1')
tadSVsDegPairs = np.load("tadSVGenePairs.txt_degPairs.npy", allow_pickle=True, encoding='latin1')
ruleSVsDegPairs = np.load("ruleSvGenePairs.txt_degPairs.npy", allow_pickle=True, encoding='latin1')

print("No of deg pairs windowed: ", windowSVsDegPairs.shape)
print("No of deg pairs tad: ", tadSVsDegPairs.shape)
print("No of deg pairs rule: ", ruleSVsDegPairs.shape)


#get the COSMIC genes

cosmicGenesFile = sys.argv[3]
cosmicGenes = []
with open(cosmicGenesFile, 'r') as f:
	lineCount = 0
	for line in f:
		if lineCount == 0:
			lineCount += 1
			continue
		
		splitLine = line.split("\t")
		
		geneName = splitLine[0]
		cosmicGenes.append(geneName)

windowedGenesCosmic = []
for gene in affectedGenesWindowed:
	if gene in cosmicGenes:
		windowedGenesCosmic.append(gene)

tadGenesCosmic = []
for gene in tadAffectedGenes:
	if gene in cosmicGenes:
		tadGenesCosmic.append(gene)

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

windowedGenesBc = []
for gene in affectedGenesWindowed:
	if gene in breastCancerGenes:
		windowedGenesBc.append(gene)

tadGenesBc = []
for gene in tadAffectedGenes:
	if gene in breastCancerGenes:
		tadGenesBc.append(gene)

ruleGenesBc = []
for gene in ruleBasedAffectedGenes:
	if gene in breastCancerGenes:
		ruleGenesBc.append(gene)

#keep the deg sv-gene pairs, which are cosmic, which are BC, and which are cosmic+deg and bc+deg
windowedCosmicDegPairs = []
windowedBcDegPairs = []
windowedCosmicPairs = []
windowedBcPairs = []

windowedDegGenes = []
windowedCosmicDegGenes = []
windowedBcDegGenes = []
for pair in svGenePairsWindowed:
	splitPair = pair.split("_")
	#check for cosmic or bc
	if splitPair[0] in cosmicGenes:
		windowedCosmicPairs.append(pair)
	if splitPair[0] in breastCancerGenes:
		windowedBcPairs.append(pair)
	
	#combinations with degs
	if pair in windowSVsDegPairs[:,0]:
		if splitPair[0] not in windowedDegGenes:
			windowedDegGenes.append(splitPair[0])
			
			if splitPair[0] in cosmicGenes:
				windowedCosmicDegPairs.append(pair)
				if splitPair[0] not in windowedCosmicDegGenes:
					windowedCosmicDegGenes.append(splitPair[0])
			if splitPair[0] in breastCancerGenes:
				windowedBcDegPairs.append(pair)
				if splitPair[0] not in windowedBcDegGenes:
					windowedBcDegGenes.append(splitPair[0])

print("Windowed no of sv-gene pairs DEG: ", windowSVsDegPairs.shape[0])
print("Windowed no of sv-gene pairs DEG and COSMIC: ", len(windowedCosmicDegPairs))
print("Windowed no of sv-gene pairs DEG and bc: ", len(windowedBcDegPairs))
print("Windowed no of sv-gene pairs COSMIC: ", len(windowedCosmicPairs))
print("Windowed no of sv-gene pairs BC: ", len(windowedBcPairs))

print("Windowed no of unique DEG genes: ", len(windowedDegGenes))
print("Windowed no of unique COSMIC genes: ", len(windowedCosmicDegGenes))
print("Windowed no of unique DEG genes: ", len(windowedBcDegGenes))

tadCosmicDegPairs = []
tadBcDegPairs = []
tadCosmicPairs = []
tadBcPairs = []

tadDegGenes = []
tadCosmicDegGenes = []
tadBcDegGenes = []
for pair in tadSVGenePairs:
	splitPair = pair.split("_")
	
	if splitPair[0] in cosmicGenes:
		tadCosmicPairs.append(pair)
	
	if splitPair[0] in breastCancerGenes:
		tadBcPairs.append(pair)
	
	if pair in tadSVsDegPairs[:,0]:
		
		if splitPair[0] not in tadDegGenes:
			tadDegGenes.append(splitPair[0])
			
		if splitPair[0] in cosmicGenes:
			tadCosmicDegPairs.append(pair)
			if splitPair[0] not in tadCosmicDegGenes:
				tadCosmicDegGenes.append(splitPair[0])
		if splitPair[0] in breastCancerGenes:
			tadBcDegPairs.append(pair)
			if splitPair[0] not in tadBcDegGenes:
				tadBcDegGenes.append(splitPair[0])

print("TAD no of sv-gene pairs DEG: ", tadSVsDegPairs.shape[0])
print("TAD no of sv-gene pairs DEG and COSMIC: ", len(tadCosmicDegPairs))
print("TAD no of sv-gene pairs DEG and bc: ", len(tadBcDegPairs))
print("TAD no of sv-gene pairs COSMIC: ", len(tadCosmicPairs))
print("TAD no of sv-gene pairs BC: ", len(tadBcPairs))

print("TAD no of unique DEG genes: ", len(tadDegGenes))
print("TAD no of unique COSMIC genes: ", len(tadCosmicDegGenes))
print("TAD no of unique DEG genes: ", len(tadBcDegGenes))

ruleCosmicDegPairs = []
ruleBcDegPairs = []
ruleCosmicPairs = []
ruleBcPairs = []

ruleDegGenes = []
ruleCosmicDegGenes = []
ruleBcDegGenes = []
for pair in ruleSvGenePairs:
	splitPair = pair.split("_")
	
	if splitPair[0] in cosmicGenes:
		ruleCosmicPairs.append(pair)
	if splitPair[0] in breastCancerGenes:
		ruleBcPairs.append(pair)
	
	if pair in ruleSVsDegPairs[:,0]:
		
		if splitPair[0] not in ruleDegGenes:
			ruleDegGenes.append(splitPair[0])
			
		if splitPair[0] in cosmicGenes:
			ruleCosmicDegPairs.append(pair)
			if splitPair[0] not in ruleCosmicDegGenes:
				ruleCosmicDegGenes.append(splitPair[0])
		if splitPair[0] in breastCancerGenes:
			ruleBcDegPairs.append(pair)
			if splitPair[0] not in ruleBcDegGenes:
				ruleBcDegGenes.append(splitPair[0])

print("Rules no of sv-gene pairs DEG: ", ruleSVsDegPairs.shape[0])
print("Rules no of sv-gene pairs DEG and COSMIC: ", len(ruleCosmicDegPairs))
print("Rules no of sv-gene pairs DEG and bc: ", len(ruleBcDegPairs))
print("Rules no of sv-gene pairs COSMIC: ", len(ruleCosmicPairs))
print("Rules no of sv-gene pairs BC: ", len(ruleBcPairs))	

print("Rules no of unique DEG genes: ", len(ruleDegGenes))
print("Rules no of unique COSMIC genes: ", len(ruleCosmicDegGenes))
print("Rules no of unique DEG genes: ", len(ruleBcDegGenes))


#Now load all the shuffled files and get the counts. Then do a t-test for the real numbers

def getAllCounts(files):
	
	#go through the files and get the number
	counts = []
	for currentFile in files:
		count = np.loadtxt(currentFile)
		counts.append(count)
	
	return counts	


import glob

shuffledPath = sys.argv[8]
# 
windowedCosmicDegCounts = getAllCounts(glob.glob(shuffledPath + 'windowedCosmicDegPairs.txt*'))
tadCosmicDegCounts = getAllCounts(glob.glob(shuffledPath + 'tadCosmicDegPairs.txt*'))
rulesCosmicDegCounts = getAllCounts(glob.glob(shuffledPath + 'rulesCosmicDegPairs.txt*'))

windowedBcDegCounts = getAllCounts(glob.glob(shuffledPath + 'windowedBcDegPairs.txt*'))
tadBcDegCounts = getAllCounts(glob.glob(shuffledPath + 'tadBcDegPairs.txt*'))
rulesBcDegCounts = getAllCounts(glob.glob(shuffledPath + 'rulesBcDegPairs.txt*'))

windowedDegCounts = getAllCounts(glob.glob(shuffledPath + 'windowedDegPairs.txt*'))
tadDegCounts = getAllCounts(glob.glob(shuffledPath + 'tadDegPairs.txt*'))
rulesDegCounts = getAllCounts(glob.glob(shuffledPath + 'rulesDegPairs.txt*'))


#Do t-tests and get the significance

#DEGs
z = (windowSVsDegPairs.shape[0] - np.mean(windowedDegCounts)) / float(np.std(windowedDegCounts))
windowDegPValue = stats.norm.sf(abs(z))*2
z = (tadSVsDegPairs.shape[0] - np.mean(tadDegCounts)) / float(np.std(tadDegCounts))
tadDegPValue = stats.norm.sf(abs(z))*2
z = (ruleSVsDegPairs.shape[0] - np.mean(rulesDegCounts)) / float(np.std(rulesDegCounts))
rulesDegPValue = stats.norm.sf(abs(z))*2

print("Windowed p-value for DEG genes: ", windowDegPValue)
print("TAD p-value for DEG genes: ", tadDegPValue)
print("Rules p-value for DEG genes: ", rulesDegPValue)

z = (len(windowedCosmicDegPairs) - np.mean(windowedCosmicDegCounts)) / float(np.std(windowedCosmicDegCounts))
windowCosmicPValue = stats.norm.sf(abs(z))*2
z = (len(tadCosmicDegPairs) - np.mean(tadCosmicDegCounts)) / float(np.std(tadCosmicDegCounts))
tadCosmicPValue = stats.norm.sf(abs(z))*2
z = (len(ruleCosmicDegPairs) - np.mean(rulesCosmicDegCounts)) / float(np.std(rulesCosmicDegCounts))
rulesCosmicPValue = stats.norm.sf(abs(z))*2

print("Windowed p-value for COSMIC+DEG sv-gene pairs: ", windowCosmicPValue)
print("TAD p-value for COSMIC+DEG sv-gene pairs: ", tadCosmicPValue)
print("Rules p-value for COSMIC+DEG sv-gene pairs: ", rulesCosmicPValue)

#BC
z = (len(windowedBcDegPairs) - np.mean(windowedBcDegCounts)) / float(np.std(windowedBcDegCounts))
windowBcPValue = stats.norm.sf(abs(z))*2
z = (len(tadBcDegPairs) - np.mean(tadBcDegCounts)) / float(np.std(tadBcDegCounts))
tadBcPValue = stats.norm.sf(abs(z))*2
z = (len(ruleBcDegPairs) - np.mean(rulesBcDegCounts)) / float(np.std(rulesBcDegCounts))
rulesBcPValue = stats.norm.sf(abs(z))*2

print("Windowed p-value for BC+DEG pairs: ", windowBcPValue)
print("TAD p-value for BC+DEG pairs: ", tadBcPValue)
print("Rules p-value for BC+DEG pairs: ", rulesBcPValue)

#Repeat but for genes, and not sv-gene pairs

print("STATS FOR GENES ONLY")
# 
windowedCosmicDegCountsGenes = getAllCounts(glob.glob(shuffledPath + 'windowedCosmicDeg.txt*'))
tadCosmicDegCountsGenes = getAllCounts(glob.glob(shuffledPath + 'tadCosmicDeg.txt*'))
rulesCosmicDegCountsGenes = getAllCounts(glob.glob(shuffledPath + 'rulesCosmicDeg.txt*'))

windowedBcDegCountsGenes = getAllCounts(glob.glob(shuffledPath + 'windowedBcDeg.txt*'))
tadBcDegCountsGenes = getAllCounts(glob.glob(shuffledPath + 'tadBcDeg.txt*'))
rulesBcDegCountsGenes = getAllCounts(glob.glob(shuffledPath + 'rulesBcDeg.txt*'))

windowedDegCountsGenes = getAllCounts(glob.glob(shuffledPath + 'windowedDeg.txt*'))
tadDegCountsGenes = getAllCounts(glob.glob(shuffledPath + 'tadDeg.txt*'))
rulesDegCountsGenes = getAllCounts(glob.glob(shuffledPath + 'rulesDeg.txt*'))


#Do t-tests and get the significance

#DEGs
z = (len(windowedDegGenes) - np.mean(windowedDegCountsGenes)) / float(np.std(windowedDegCountsGenes))
windowDegPValue = stats.norm.sf(abs(z))*2
z = (len(tadDegGenes) - np.mean(tadDegCountsGenes)) / float(np.std(tadDegCountsGenes))
tadDegPValue = stats.norm.sf(abs(z))*2
z = (len(ruleDegGenes) - np.mean(rulesDegCountsGenes)) / float(np.std(rulesDegCountsGenes))
rulesDegPValue = stats.norm.sf(abs(z))*2

print("Windowed p-value for DEG genes: ", windowDegPValue)
print("TAD p-value for DEG genes: ", tadDegPValue)
print("Rules p-value for DEG genes: ", rulesDegPValue)

z = (len(windowedCosmicDegGenes) - np.mean(windowedCosmicDegCountsGenes)) / float(np.std(windowedCosmicDegCountsGenes))
windowCosmicPValue = stats.norm.sf(abs(z))*2
z = (len(tadCosmicDegGenes) - np.mean(tadCosmicDegCountsGenes)) / float(np.std(tadCosmicDegCountsGenes))
tadCosmicPValue = stats.norm.sf(abs(z))*2
z = (len(ruleCosmicDegGenes) - np.mean(rulesCosmicDegCountsGenes)) / float(np.std(rulesCosmicDegCountsGenes))
rulesCosmicPValue = stats.norm.sf(abs(z))*2

print("Windowed p-value for COSMIC+DEG sv-gene pairs: ", windowCosmicPValue)
print("TAD p-value for COSMIC+DEG sv-gene pairs: ", tadCosmicPValue)
print("Rules p-value for COSMIC+DEG sv-gene pairs: ", rulesCosmicPValue)

#BC
z = (len(windowedBcDegGenes) - np.mean(windowedBcDegCountsGenes)) / float(np.std(windowedBcDegCountsGenes))
windowBcPValue = stats.norm.sf(abs(z))*2
z = (len(tadBcDegGenes) - np.mean(tadBcDegCountsGenes)) / float(np.std(tadBcDegCountsGenes))
tadBcPValue = stats.norm.sf(abs(z))*2
z = (len(ruleBcDegGenes) - np.mean(rulesBcDegCountsGenes)) / float(np.std(rulesBcDegCountsGenes))
rulesBcPValue = stats.norm.sf(abs(z))*2

print("Windowed p-value for BC+DEG pairs: ", windowBcPValue)
print("TAD p-value for BC+DEG pairs: ", tadBcPValue)
print("Rules p-value for BC+DEG pairs: ", rulesBcPValue)





