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
from scipy.stats import chi2_contingency

#1. Filter SVs based on if these cause DEGs or affect COSMIC genes in the coding way

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
					
	print "Number of genes affected in the coding way: ", len(genesAffectedByCodingSVs)
	print "Number of genes affected in the coding way that are DEG or COSMIC: ", len(genesAffectedByFilteredSVs)
	return codingEffectSVs
	
codingEffectSVs = getSVsWithCodingEffects()
print "Number of SVs filtered out with coding effects: ", len(codingEffectSVs)
np.savetxt('codingEffectSVs.txt', codingEffectSVs, delimiter='\t', fmt='%s')


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
				svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[7]
			svGenePairs.append(gene[3].name + "_" + svStr)
	
	return affectedGenes, svGenePairs
	
affectedGenesWindowed, svGenePairsWindowed = findAffectedGenesWithinWindow()
np.savetxt("Output/windowedSVs.txt", svGenePairsWindowed, delimiter='\t', fmt='%s')
print "Number of affected genes in windowed approach: ", len(affectedGenesWindowed)

#3. Find all genes within the TAD of the filtered SVs
def findAffectedGenesByTadDisruptions(codingEffectSVs):
		
	somaticSVs = InputParser().getSVsFromFile(sys.argv[4], "all", codingEffectSVs)
	tads = InputParser().getTADsFromFile(sys.argv[5])
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
	svGenePairs = []
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
				
				svGenePairs.append(gene[3].name + "_" + svStr)
	
	return affectedGenes, svGenePairs

tadAffectedGenes, tadSVGenePairs = findAffectedGenesByTadDisruptions(codingEffectSVs)
np.savetxt("Output/tadSVs.txt", tadSVGenePairs, delimiter='\t', fmt='%s')

print "TAD affected genes: ", len(tadAffectedGenes)

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
ruleSvGenePairs = np.loadtxt('Output/geneSVPairs_somatic_me_12072019_shuffled.txt_none', dtype='object')
np.savetxt('Output/ruleSVs.txt', ruleSvGenePairs[:,0], delimiter='\t', fmt='%s')


#5. Compare the resulting genes between the approaches

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
		
print "Number of cosmic genes in the windowed approach: ", len(windowedGenesCosmic)
print "Number of cosmic genes in the tad approach: ", len(tadGenesCosmic)
print "Number of cosmic genes in the rule approach: ", len(ruleGenesCosmic)

#Compute the chi2 p-values for these findings
#Because we are looking at all other genes, the number of cosmic genes - genes in the true group is the negative.
# 
obs = np.array([[len(windowedGenesCosmic), len(cosmicGenes) - len(windowedGenesCosmic)], [len(affectedGenesWindowed) - len(windowedGenesCosmic), (19286 - len(affectedGenesWindowed)- (len(cosmicGenes) - len(windowedGenesCosmic)))]])
print obs
g, p, dof, expctd = chi2_contingency(obs)
print "COSMIC p-value windowed: ", p

obs = np.array([[len(tadGenesCosmic), len(cosmicGenes) - len(tadGenesCosmic)], [len(tadAffectedGenes) - len(tadGenesCosmic), (19286 - len(tadAffectedGenes) - (len(cosmicGenes) - len(tadGenesCosmic)))]])
g, p, dof, expctd = chi2_contingency(obs)
print "COSMIC p-value tad: ", p

obs = np.array([[len(ruleGenesCosmic), len(cosmicGenes) - len(ruleGenesCosmic)], [len(ruleBasedAffectedGenes) - len(ruleGenesCosmic), (19286 - len(ruleBasedAffectedGenes) - (len(cosmicGenes) - len(ruleGenesCosmic)))]])
g, p, dof, expctd = chi2_contingency(obs)
print "COSMIC p-value rules: ", p
#To get the DEGs, we actually need to re-compute the DEGs based on the gene-SV pairs that we get for each method. Otherwise we are biasing towards the rule-based approach. 

#For now, simply load in the data
windowSVsDegPairs = np.load('Output/ShuffledCodingNonCoding/geneCodingSVPairs_somatic_me_12072019_shuffled.txt__windowedSVs.txt_degPairs.npy', allow_pickle=True)
tadSVsDegPairs = np.load('Output/ShuffledCodingNonCoding/geneCodingSVPairs_somatic_me_12072019_shuffled.txt__tadSVs.txt_degPairs.npy', allow_pickle=True)
ruleSVsDegPairs = np.load('Output/ShuffledCodingNonCoding/geneCodingSVPairs_somatic_me_12072019_shuffled.txt__ruleSVs.txt_degPairs.npy', allow_pickle=True)


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
			
#For the rule based method, we need to look at a separate output file to get the gene-SV pairs. This may be fixed in the actual tool output later.
ruleSvGenePairs = np.loadtxt('Output/geneSVPairs_somatic_me_12072019_shuffled.txt_none', dtype='object')

ruleDegGenes = []
for pair in ruleSvGenePairs[:,0]:
	if pair in ruleSVsDegPairs[:,0]:
		splitPair = pair.split("_")
		if splitPair[0] not in ruleDegGenes:
			ruleDegGenes.append(splitPair[0])
			
print "Number of DEG genes in the windowed approach: ", len(windowedDegGenes)
print "Number of DEG genes in the TAD approach: ", len(tadDegGenes)
print "number of DEG genes in the rule approach: ", len(ruleDegGenes)

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
print "Number of genes that are in all 3 nc-based methods: ", len(allCriteriaIntersect)
print "Number of genes that are in the windowed and TAD approaches: ", len(windowTadIntersect)
print "Number of genes that are in windowed and rule approaches: ", len(windowRuleIntersect)
print "Number of genes that are in the TAD and rule approahes: ", len(tadRuleIntersect)

import pylab as plt
from matplotlib_venn import venn3, venn3_circles

v = venn3(subsets=(len(affectedGenesWindowed),len(tadAffectedGenes),len(windowTadIntersect), len(ruleBasedAffectedGenes),len(windowRuleIntersect), len(tadRuleIntersect),len(allCriteriaIntersect)),
		  set_labels=('Windowed', 'TAD', 'Rules'))
v.get_label_by_id('100').set_text(len(affectedGenesWindowed))
v.get_label_by_id('010').set_text(len(tadAffectedGenes))
v.get_label_by_id('001').set_text(len(ruleBasedAffectedGenes))
plt.title("Overlapping genes in the nc-based approaches")
plt.show()

#Repeat for COSMIC

allCriteriaIntersect = list(set(windowedGenesCosmic) & set(tadGenesCosmic) & set(ruleGenesCosmic))
windowTadIntersect = list(set(windowedGenesCosmic) & set(tadGenesCosmic))
windowRuleIntersect = list(set(windowedGenesCosmic) & set(ruleGenesCosmic))
tadRuleIntersect = list(set(tadGenesCosmic) & set(ruleGenesCosmic))
print "Number of genes that are in all 3 nc-based methods: ", len(allCriteriaIntersect)
print "Number of genes that are in the windowed and TAD approaches: ", len(windowTadIntersect)
print "Number of genes that are in windowed and rule approaches: ", len(windowRuleIntersect)
print "Number of genes that are in the TAD and rule approahes: ", len(tadRuleIntersect)

import pylab as plt
from matplotlib_venn import venn3, venn3_circles

v = venn3(subsets=(len(windowedGenesCosmic),len(tadGenesCosmic),len(windowTadIntersect), len(ruleGenesCosmic),len(windowRuleIntersect), len(tadRuleIntersect),len(allCriteriaIntersect)),
		  set_labels=('Windowed', 'TAD', 'Rules'))
v.get_label_by_id('100').set_text(len(windowedGenesCosmic))
v.get_label_by_id('010').set_text(len(tadGenesCosmic))
v.get_label_by_id('001').set_text(len(ruleGenesCosmic))
plt.title("Overlapping COSMIC in the nc-based approaches")
plt.show()

#Repeat for DEG genes

allCriteriaIntersect = list(set(windowedDegGenes) & set(tadDegGenes) & set(ruleDegGenes))
windowTadIntersect = list(set(windowedDegGenes) & set(tadDegGenes))
windowRuleIntersect = list(set(windowedDegGenes) & set(ruleDegGenes))
tadRuleIntersect = list(set(tadDegGenes) & set(ruleDegGenes))
print "Number of genes that are in all 3 nc-based methods: ", len(allCriteriaIntersect)
print "Number of genes that are in the windowed and TAD approaches: ", len(windowTadIntersect)
print "Number of genes that are in windowed and rule approaches: ", len(windowRuleIntersect)
print "Number of genes that are in the TAD and rule approahes: ", len(tadRuleIntersect)

import pylab as plt
from matplotlib_venn import venn3, venn3_circles

v = venn3(subsets=(len(windowedDegGenes),len(tadDegGenes),len(windowTadIntersect), len(ruleDegGenes),len(windowRuleIntersect), len(tadRuleIntersect),len(allCriteriaIntersect)),
		  set_labels=('Windowed', 'TAD', 'Rules'))
v.get_label_by_id('100').set_text(len(windowedDegGenes))
v.get_label_by_id('010').set_text(len(tadDegGenes))
v.get_label_by_id('001').set_text(len(ruleDegGenes))
plt.title("Overlapping DEG genes in the nc-based approaches")
plt.show()




