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

# 1. Filter SVs based on if these cause DEGs or affect COSMIC genes in the coding way

def getSVsWithCodingEffects():
	
	codingPairs = np.loadtxt(sys.argv[1], dtype='object')
	degPairs = np.load(sys.argv[2], allow_pickle=True)
	cosmicGenesFile = sys.argv[3]
	
	codingSVGenes = dict()
	for pair in codingPairs:
		
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
	genesAffectedByFilteredSVs = []
	genesAffectedByCodingSVs = []
	codingEffectPairs = []
	for sv in codingSVGenes:
		
		degCount = 0
		cosmicCount = 0
		degAndCosmicCount = 0
		degOrCosmicCount = 0
		
		for gene in codingSVGenes[sv]:
			pair = gene + "_" + sv
			
			if gene not in genesAffectedByCodingSVs:
				genesAffectedByCodingSVs.append(gene)
			
			if gene in cosmicGenes or pair in degPairs[:,0]:
				if sv not in codingEffectSVs:
					codingEffectSVs.append(sv)
				if gene not in genesAffectedByFilteredSVs: #look at all genes that are linked to the SVs. 
					genesAffectedByFilteredSVs.append(gene)
					
				if pair not in codingEffectPairs:
					codingEffectPairs.append(pair)
					
	#print "Number of genes affected in the coding way: ", len(genesAffectedByCodingSVs)
	#print "Number of genes affected in the coding way that are DEG or COSMIC: ", len(genesAffectedByFilteredSVs)
	return codingEffectSVs, codingEffectPairs
	
codingEffectSVs, codingEffectPairs = getSVsWithCodingEffects()
#print "Number of SVs filtered out with coding effects: ", len(codingEffectSVs)
np.savetxt('codingEffectSVs.txt', codingEffectSVs, delimiter='\t', fmt='%s')
np.savetxt('codingEffectPairs.txt', codingEffectPairs, delimiter='\t', fmt='%s')


#2. Find all genes within a window of the filtered SVs
def findAffectedGenesWithinWindow():
	
	#Read all SVs and filter these for the coding effect SVs
	somaticSVs = InputParser().getSVsFromFile(sys.argv[4], "all", codingEffectSVs)
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

	#Specify a window to look around SVs. 2 mb seems good as I recall from quite some time ago. But this threshold may change
	window = 2000000
	
	affectedGenes = []
	svGenePairs = []
	#For every SV, look at 2 mb to the left of the left breakpoint, and look at 2 mb to the right of the right breakpoint.
	for sv in somaticSVs:
		
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
					svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
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
					svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
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
					svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
				svGenePairs.append(gene[3].name + "_" + svStr)
			
			
	return affectedGenes, svGenePairs
	
affectedGenesWindowed, svGenePairsWindowed = findAffectedGenesWithinWindow()
np.savetxt("Output/windowedSVs.txt", svGenePairsWindowed, delimiter='\t', fmt='%s')
print("Number of affected genes in windowed approach: ", len(affectedGenesWindowed))


###Alternative to look at disrupted TADs, not just boundaries
def findAffectedGenesByTadDisruptions(codingEffectSVs):
		
	somaticSVs = InputParser().getSVsFromFile(sys.argv[4], "all", codingEffectSVs)
	tads = InputParser().getTADsFromFile(sys.argv[5])
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
	
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
						
						svGenePairs.append(gene[3].name + "_" + svStr)

	return affectedGenes, svGenePairs

tadAffectedGenes, tadSVGenePairs = findAffectedGenesByTadDisruptions(codingEffectSVs)
np.savetxt("Output/tadSVs.txt", tadSVGenePairs, delimiter='\t', fmt='%s')

print("TAD affected genes: ", len(tadAffectedGenes))


#4. Run the rule-based method on the filtered SVs

#here I wil bypass main for simplicity, but could be much neater I think. the issue is passing the exlucded SVs.
def getGenesWithRuleBasedApproach():
	causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
	
	#Combine the genes into one set. 
	causalGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
	
	svData = InputParser().getSVsFromFile(sys.argv[4], "all", codingEffectSVs)
	
	NeighborhoodDefiner(causalGenes, svData, None, 'SV') #Provide the mode to ensure that the right variant type is used (different positions used in annotation)
	
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

ruleBasedAffectedGenes = getGenesWithRuleBasedApproach()
print("rule-based affected genes: ", len(ruleBasedAffectedGenes))
ruleSvGenePairs = np.loadtxt('Output/RankedGenes/naive/BRCA//nonCoding_geneSVPairs.txt_none', dtype='object')


np.savetxt('Output/ruleSVs.txt', ruleSvGenePairs[:,0], delimiter='\t', fmt='%s')

#Save all genes in memory to prevent re-computing every time
np.savetxt('affectedGenesWindowed.txt', affectedGenesWindowed, delimiter='\t', fmt='%s')
np.savetxt('tadAffectedGenes.txt', tadAffectedGenes, delimiter='\t', fmt='%s')
np.savetxt('ruleBasedAffectedGenes.txt', ruleBasedAffectedGenes, delimiter='\t', fmt='%s')

np.savetxt('svGenePairsWindowed.txt', svGenePairsWindowed, delimiter='\t', fmt='%s')
np.savetxt('tadSVGenePairs.txt', tadSVGenePairs, delimiter='\t', fmt='%s')
np.savetxt('ruleSvGenePairs.txt', ruleSvGenePairs[:,0], delimiter='\t', fmt='%s')

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

# #Compute the chi2 p-values for these findings
# #Because we are looking at all other genes, the number of cosmic genes - genes in the true group is the negative.
# # 
# obs = np.array([[len(windowedGenesBc), len(breastCancerGenes) - len(windowedGenesBc)], [len(affectedGenesWindowed) - len(windowedGenesBc), (19286 - len(affectedGenesWindowed)- (len(breastCancerGenes) - len(windowedGenesBc)))]])
# print(obs)
# g, p, dof, expctd = chi2_contingency(obs)
# print("BC p-value windowed: ", p)
# 
# obs = np.array([[len(tadGenesBc), len(breastCancerGenes) - len(tadGenesBc)], [len(tadAffectedGenes) - len(tadGenesBc), (19286 - len(tadAffectedGenes) - (len(breastCancerGenes) - len(tadGenesBc)))]])
# g, p, dof, expctd = chi2_contingency(obs)
# print("BC p-value tad: ", p)
# 
# obs = np.array([[len(ruleGenesBc), len(breastCancerGenes) - len(ruleGenesBc)], [len(ruleBasedAffectedGenes) - len(ruleGenesBc), (19286 - len(ruleBasedAffectedGenes) - (len(breastCancerGenes) - len(ruleGenesBc)))]])
# g, p, dof, expctd = chi2_contingency(obs)
# print("BC p-value rules: ", p)



#To get the DEGs, we actually need to re-compute the DEGs based on the gene-SV pairs that we get for each method. Otherwise we are biasing towards the rule-based approach. 
#Collect the DEG genes for each SV-gene pair combination
#The DEGs  will need to be re-computed for each shuffled iteration
# windowExprCall = "python computeSVGenePairExpression_oneSet.py svGenePairsWindowed.txt " + sys.argv[8] 
# os.system(windowExprCall)
# tadExprCall = "python computeSVGenePairExpression_oneSet.py tadSVGenePairs.txt " + sys.argv[8] 
# os.system(tadExprCall)
# rulesExprCall = "python computeSVGenePairExpression_oneSet.py ruleSvGenePairs.txt " + sys.argv[8] 
# os.system(rulesExprCall)

# Read the DEG pairs and determine how many genes are DEG in total
svGenePairsWindowed = np.loadtxt("Output/windowedSVs.txt", dtype='object')
windowSVsDegPairs = np.load("svGenePairsWindowed.txt_degPairs.npy", allow_pickle=True, encoding='latin1')
tadSVsDegPairs = np.load("tadSVGenePairs.txt_degPairs.npy", allow_pickle=True, encoding='latin1')
ruleSVsDegPairs = np.load("ruleSvGenePairs.txt_degPairs.npy", allow_pickle=True, encoding='latin1')

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
		

#Get the breast cancer specific genes
breastCancerGenesFile = sys.argv[6]
breastCancerGenes = []
with open(breastCancerGenesFile, 'r') as f:
	
	for line in f:
		
		line = line.strip()
		
		breastCancerGenes.append(line)

			
windowedDegGenes = []
windowedCosmicDegGenes = []
windowedBcDegGenes = []
for pair in svGenePairsWindowed:
	if pair in windowSVsDegPairs[:,0]:
		splitPair = pair.split("_")
		if splitPair[0] not in windowedDegGenes:
			windowedDegGenes.append(splitPair[0])
			
			if splitPair[0] in cosmicGenes:
				if splitPair[0] not in windowedCosmicDegGenes:
					windowedCosmicDegGenes.append(splitPair[0])
			if splitPair[0] in breastCancerGenes:
				if splitPair[0] not in windowedBcDegGenes:
					windowedBcDegGenes.append(splitPair[0])
				

tadDegGenes = []
tadCosmicDegGenes = []
tadBcDegGenes = []
for pair in tadSVGenePairs:
	if pair in tadSVsDegPairs[:,0]:
		splitPair = pair.split("_")
		if splitPair[0] not in tadDegGenes:
			tadDegGenes.append(splitPair[0])
			
		if splitPair[0] in cosmicGenes:
			if splitPair[0] not in tadCosmicDegGenes:
				tadCosmicDegGenes.append(splitPair[0])
		if splitPair[0] in breastCancerGenes:
			if splitPair[0] not in tadBcDegGenes:
				tadBcDegGenes.append(splitPair[0])

ruleDegGenes = []
ruleCosmicDegGenes = []
ruleBcDegGenes = []
for pair in ruleSvGenePairs[:,0]:
	if pair in ruleSVsDegPairs[:,0]:
		splitPair = pair.split("_")
		if splitPair[0] not in ruleDegGenes:
			ruleDegGenes.append(splitPair[0])
			
		if splitPair[0] in cosmicGenes:
			if splitPair[0] not in ruleCosmicDegGenes:
				ruleCosmicDegGenes.append(splitPair[0])
		if splitPair[0] in breastCancerGenes:
			if splitPair[0] not in ruleBcDegGenes:
				ruleBcDegGenes.append(splitPair[0])
			
print("Number of DEG genes in the windowed approach: ", len(windowedDegGenes))
print("Number of DEG genes in the TAD approach: ", len(tadDegGenes))
print("number of DEG genes in the rule approach: ", len(ruleDegGenes))

#For each set, how many of the genes are in COSMIC
		
print("Number of cosmic genes in the windowed approach: ", len(windowedCosmicDegGenes))
print("Number of cosmic genes in the tad approach: ", len(tadCosmicDegGenes))
print("Number of cosmic genes in the rule approach: ", len(ruleCosmicDegGenes))


#Compute the chi2 p-values for these findings
#Because we are looking at all other genes, the number of cosmic genes - genes in the true group is the negative.
# 
# obs = np.array([[len(windowedGenesCosmic), len(cosmicGenes) - len(windowedGenesCosmic)], [len(affectedGenesWindowed) - len(windowedGenesCosmic), (19286 - len(affectedGenesWindowed)- (len(cosmicGenes) - len(windowedGenesCosmic)))]])
# print(obs)
# g, p, dof, expctd = chi2_contingency(obs)
# print("COSMIC p-value windowed: ", p)
# 
# obs = np.array([[len(tadGenesCosmic), len(cosmicGenes) - len(tadGenesCosmic)], [len(tadAffectedGenes) - len(tadGenesCosmic), (19286 - len(tadAffectedGenes) - (len(cosmicGenes) - len(tadGenesCosmic)))]])
# g, p, dof, expctd = chi2_contingency(obs)
# print("COSMIC p-value tad: ", p)
# 
# obs = np.array([[len(ruleGenesCosmic), len(cosmicGenes) - len(ruleGenesCosmic)], [len(ruleBasedAffectedGenes) - len(ruleGenesCosmic), (19286 - len(ruleBasedAffectedGenes) - (len(cosmicGenes) - len(ruleGenesCosmic)))]])
# g, p, dof, expctd = chi2_contingency(obs)
# print("COSMIC p-value rules: ", p)

		
print("Number of bc genes in the windowed approach: ", len(windowedBcDegGenes))
print("Number of bc genes in the tad approach: ", len(tadBcDegGenes))
print("Number of bc genes in the rule approach: ", len(ruleBcDegGenes))


#Compute enrichment for gene sets and GO terms
# print("doing enrichment: ")
# 
# import gseapy as gp
# import pandas as pd
# pd.set_option('display.max_columns', 500)
# import matplotlib.pyplot as plt
# from matplotlib import interactive
# interactive(True)
# 
# # 
# # bm = Biomart(verbose=False, host="asia.ensembl.org")
# # results = bm.query(dataset='hsapiens_gene_ensembl',
# #                    attributes=['external_gene_name','entrezgene', 'go_id'],
# #                    filters={'hgnc_symbol': ruleBasedAffectedGenes},
# #                    # save output file
# #                    filename="ruleBasedBmIDs.txt")
# # 
# # print results.head()
# 
# enr = gp.enrichr(gene_list=list(ruleDegGenes),
#                  description='ruleBasedGenes',
#                  gene_sets=['KEGG_2016','KEGG_2013'],
#                  outdir='ruleBased/enrichr_kegg'
#                 )
# 
# print(enr.results.head())
# 
# from gseapy.plot import barplot, dotplot
# a = barplot(enr.res2d,title='KEGG_2013',ofname='test.png')
# print(a)
# exit()


#Now load all the shuffled files and get the counts. Then do a t-test for the real numbers

def getAllCounts(files):
	
	#go through the files and get the number
	counts = []
	for currentFile in files:
		count = np.loadtxt(currentFile)
		counts.append(count)
	
	return counts	


import glob

shuffledPath = sys.argv[7]
# 
windowedCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'windowedCosmicDeg.txt*'))
tadCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'tadCosmicDeg.txt*'))
rulesCosmicCounts = getAllCounts(glob.glob(shuffledPath + 'rulesCosmicDeg.txt*'))

# print "no of genes in the cosmic case for rules: ", len(ruleGenesCosmic)
# plt.hist(rulesCosmicCounts)
# plt.show()
# plt.clf()

windowedBcCounts = getAllCounts(glob.glob(shuffledPath + 'windowedBcDeg.txt*'))
tadBcCounts = getAllCounts(glob.glob(shuffledPath + 'tadBcDeg.txt*'))
rulesBcCounts = getAllCounts(glob.glob(shuffledPath + 'rulesBcDeg.txt*'))

# print "no of genes in the cosmic case for rules: ", len(ruleGenesBc)
# plt.hist(rulesBcCounts)
# plt.show()
# plt.clf()

windowedDegCounts = getAllCounts(glob.glob(shuffledPath + 'windowedDeg.txt*'))
tadDegCounts = getAllCounts(glob.glob(shuffledPath + 'tadDeg.txt*'))
rulesDegCounts = getAllCounts(glob.glob(shuffledPath + 'rulesDeg.txt*'))

# print "no of genes in the cosmic case for rules: ", len(ruleDegGenes)
# plt.hist(rulesDegCounts)
# plt.show()
# plt.clf()

#Do t-tests and get the significance


z = (len(windowedCosmicDegGenes) - np.mean(windowedCosmicCounts)) / float(np.std(windowedCosmicCounts))
windowCosmicPValue = stats.norm.sf(abs(z))*2
z = (len(tadCosmicDegGenes) - np.mean(tadCosmicCounts)) / float(np.std(tadCosmicCounts))
tadCosmicPValue = stats.norm.sf(abs(z))*2
z = (len(ruleCosmicDegGenes) - np.mean(rulesCosmicCounts)) / float(np.std(rulesCosmicCounts))
rulesCosmicPValue = stats.norm.sf(abs(z))*2

print("Windowed p-value for COSMIC genes: ", windowCosmicPValue)
print("TAD p-value for COSMIC genes: ", tadCosmicPValue)
print("Rules p-value for COSMIC genes: ", rulesCosmicPValue)

#BC
z = (len(windowedBcDegGenes) - np.mean(windowedBcCounts)) / float(np.std(windowedBcCounts))
windowBcPValue = stats.norm.sf(abs(z))*2
z = (len(tadBcDegGenes) - np.mean(tadBcCounts)) / float(np.std(tadBcCounts))
tadBcPValue = stats.norm.sf(abs(z))*2
z = (len(ruleBcDegGenes) - np.mean(rulesBcCounts)) / float(np.std(rulesBcCounts))
rulesBcPValue = stats.norm.sf(abs(z))*2

print("Windowed p-value for BC genes: ", windowBcPValue)
print("TAD p-value for BC genes: ", tadBcPValue)
print("Rules p-value for BC genes: ", rulesBcPValue)

#DEGs
z = (len(windowedDegGenes) - np.mean(windowedDegCounts)) / float(np.std(windowedDegCounts))
windowDegPValue = stats.norm.sf(abs(z))*2
z = (len(tadDegGenes) - np.mean(tadDegCounts)) / float(np.std(tadDegCounts))
tadDegPValue = stats.norm.sf(abs(z))*2
z = (len(ruleDegGenes) - np.mean(rulesDegCounts)) / float(np.std(rulesDegCounts))
rulesDegPValue = stats.norm.sf(abs(z))*2

print("Windowed p-value for DEG genes: ", windowDegPValue)
print("TAD p-value for DEG genes: ", tadDegPValue)
print("Rules p-value for DEG genes: ", rulesDegPValue)

exit()
##Do enrichment for gene sets/Go terms






#Compute the chi2 p-values for these findings
# obs = np.array([[len(windowedDegGenes), len(cosmicGenes) - len(windowedDegGenes)], [len(affectedGenesWindowed) - len(windowedDegGenes), (19286 - len(affectedGenesWindowed)- (len(cosmicGenes) - len(windowedDegGenes)))]])
# g, p, dof, expctd = chi2_contingency(obs)
# print "DEG p-value windowed: ", p
# 
# obs = np.array([[len(tadDegGenes), len(cosmicGenes) - len(tadDegGenes)], [len(tadAffectedGenes) - len(tadDegGenes), (19286 - len(tadAffectedGenes) - (len(cosmicGenes) - len(tadDegGenes)))]])
# g, p, dof, expctd = chi2_contingency(obs)
# print "DEG p-value tad: ", p
# 
# obs = np.array([[len(ruleDegGenes), len(cosmicGenes) - len(ruleDegGenes)], [len(ruleBasedAffectedGenes) - len(ruleDegGenes), (19286 - len(ruleBasedAffectedGenes) - (len(cosmicGenes) - len(ruleDegGenes)))]])
# g, p, dof, expctd = chi2_contingency(obs)
# print "DEG p-value rules: ", p

#Make a venn diagram

#For the genes that are found
allCriteriaIntersect = list(set(affectedGenesWindowed) & set(tadAffectedGenes) & set(ruleBasedAffectedGenes))
windowTadIntersect = list(set(affectedGenesWindowed) & set(tadAffectedGenes))
windowRuleIntersect = list(set(affectedGenesWindowed) & set(ruleBasedAffectedGenes))
tadRuleIntersect = list(set(tadAffectedGenes) & set(ruleBasedAffectedGenes))

print("Number of genes that are in all 3 nc-based methods: ", len(allCriteriaIntersect))
print("Number of genes that are in the windowed and TAD approaches: ", len(windowTadIntersect))
print("Number of genes that are in windowed and rule approaches: ", len(windowRuleIntersect))
print("Number of genes that are in the TAD and rule approahes: ", len(tadRuleIntersect))

#Determine the venn diagram values
twr = len(allCriteriaIntersect)
tw = len(windowTadIntersect) - twr
wr = len(windowRuleIntersect) - twr
tr = len(tadRuleIntersect) - twr
t = len(tadAffectedGenes) - twr - tw - tr
r = len(ruleBasedAffectedGenes) - twr - wr - tr
w = len(affectedGenesWindowed) - twr - wr - tw

v = venn3(subsets=(w, t,tw, r,wr, tr,twr),
		  set_labels=('Windowed', 'TAD', 'Rules'))
#v.get_label_by_id('100').set_text(w)
#v.get_label_by_id('010').set_text(t)
#v.get_label_by_id('001').set_text(r)
plt.title("Overlapping genes in the nc-based approaches")
plt.savefig('genesOverlap.svg')
#plt.show()
plt.clf()

#Repeat for COSMIC

allCriteriaIntersect = list(set(windowedGenesCosmic) & set(tadGenesCosmic) & set(ruleGenesCosmic))
windowTadIntersect = list(set(windowedGenesCosmic) & set(tadGenesCosmic))
windowRuleIntersect = list(set(windowedGenesCosmic) & set(ruleGenesCosmic))
tadRuleIntersect = list(set(tadGenesCosmic) & set(ruleGenesCosmic))
print("Number of genes that are in all 3 nc-based methods: ", len(allCriteriaIntersect))
print("Number of genes that are in the windowed and TAD approaches: ", len(windowTadIntersect))
print("Number of genes that are in windowed and rule approaches: ", len(windowRuleIntersect))
print("Number of genes that are in the TAD and rule approahes: ", len(tadRuleIntersect))


twr = len(allCriteriaIntersect)
tw = len(windowTadIntersect) - twr
wr = len(windowRuleIntersect) - twr
tr = len(tadRuleIntersect) - twr
t = len(tadGenesCosmic) - twr - tw - tr
r = len(ruleGenesCosmic) - twr - wr - tr
w = len(windowedGenesCosmic) - twr - wr - tw

v = venn3(subsets=(w, t,tw, r,wr, tr,twr),
		  set_labels=('Windowed', 'TAD', 'Rules'))
plt.title("Overlapping COSMIC genes in the nc-based approaches")
plt.savefig('cosmicOverlap.svg')
#plt.show()
plt.clf()

#repeat for bc genes
allCriteriaIntersect = list(set(windowedGenesBc) & set(tadGenesBc) & set(ruleGenesBc))
windowTadIntersect = list(set(windowedGenesBc) & set(tadGenesBc))
windowRuleIntersect = list(set(windowedGenesBc) & set(ruleGenesBc))
tadRuleIntersect = list(set(tadGenesBc) & set(ruleGenesBc))
print("Number of genes that are in all 3 nc-based methods: ", len(allCriteriaIntersect))
print("Number of genes that are in the windowed and TAD approaches: ", len(windowTadIntersect))
print("Number of genes that are in windowed and rule approaches: ", len(windowRuleIntersect))
print("Number of genes that are in the TAD and rule approahes: ", len(tadRuleIntersect))


twr = len(allCriteriaIntersect)
tw = len(windowTadIntersect) - twr
wr = len(windowRuleIntersect) - twr
tr = len(tadRuleIntersect) - twr
t = len(tadGenesBc) - twr - tw - tr
r = len(ruleGenesBc) - twr - wr - tr
w = len(windowedGenesBc) - twr - wr - tw

v = venn3(subsets=(w, t,tw, r,wr, tr,twr),
		  set_labels=('Windowed', 'TAD', 'Rules'))
plt.title("Overlapping breast cancer genes in the nc-based approaches")
plt.savefig('bcOverlap.svg')
#plt.show()
plt.clf()


#Repeat for DEG genes

allCriteriaIntersect = list(set(windowedDegGenes) & set(tadDegGenes) & set(ruleDegGenes))
windowTadIntersect = list(set(windowedDegGenes) & set(tadDegGenes))
windowRuleIntersect = list(set(windowedDegGenes) & set(ruleDegGenes))
tadRuleIntersect = list(set(tadDegGenes) & set(ruleDegGenes))

print("Number of genes that are in all 3 nc-based methods: ", len(allCriteriaIntersect))
print("Number of genes that are in the windowed and TAD approaches: ", len(windowTadIntersect))
print("Number of genes that are in windowed and rule approaches: ", len(windowRuleIntersect))
print("Number of genes that are in the TAD and rule approahes: ", len(tadRuleIntersect))

twr = len(allCriteriaIntersect)
tw = len(windowTadIntersect) - twr
wr = len(windowRuleIntersect) - twr
tr = len(tadRuleIntersect) - twr
t = len(tadDegGenes) - twr - tw - tr
r = len(ruleDegGenes) - twr - wr - tr
w = len(windowedDegGenes) - twr - wr - tw

v = venn3(subsets=(w, t,tw, r,wr, tr,twr),
		  set_labels=('Windowed', 'TAD', 'Rules'))
plt.title("Overlapping DEG genes in the nc-based approaches")
#plt.show()
plt.savefig('degOverlap.svg')




