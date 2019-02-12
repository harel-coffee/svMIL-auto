"""
	Script to go through the list of ranked genes, and determine at which threshold of number of samples we find more COSMIC genes, genes with SNVs or DEGs than when SVs are randomized. 

"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import os
from os import listdir
from os.path import isfile, join

rankedGenesFileNum = int(sys.argv[1])-1
cosmicGenesFile = sys.argv[2]
snvFile = sys.argv[3]
degFile = sys.argv[4]
permutationDataFolder = sys.argv[5]

#Get the files in the folder
geneScoreFiles = [f for f in listdir(permutationDataFolder) if isfile(join(permutationDataFolder, f))]


rankedGenesFile = geneScoreFiles[rankedGenesFileNum]

#1. Get the ranked genes 
#Read the ranking file
geneScores = dict()
maxScore = 0 #Keep a max score that we use to determine the thresholds to search through
with open(permutationDataFolder + "/" + rankedGenesFile, 'r') as rankF:
	
	for line in rankF:
		line = line.strip()
		splitLine = line.split()
		
		geneName = splitLine[0]
		geneScore = int(float(splitLine[1]))
		if geneScore > maxScore: 
			maxScore = geneScore
		
		geneScores[geneName] = geneScore

#2. Get the COSMIC genes, genes with SNVs, and DEGs

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
		
#Then, determine how many of the genes are differentially expressed
degGenes = []
with open(degFile, 'r') as inf:
	for line in inf:
		line = line.strip()
		degGenes.append(line)
		
snvGenes = []
with open(snvFile, 'r') as inf:
	for line in inf:
		line = line.strip()
		snvGenes.append(line)

#3. At every threshold, compute how many genes are in the 3 categories (also overlapping)



#Compute the intersects
def computeCategoryMatches(realGeneScores, threshold):
	#print "threshold: ", threshold
	degGenesPos = []
	degGenesNeg = []
	cosmicGenesPos = []
	cosmicGenesNeg = []
	snvGenesPos = []
	snvGenesNeg = []
	
	for gene in realGeneScores:
		
		if float(gene[1]) > threshold:
			if gene[0] in degGenes:
				degGenesPos.append(gene[0])
			if gene[0] in cosmicGenes:
				cosmicGenesPos.append(gene[0])
			if gene[0] in snvGenes:
				snvGenesPos.append(gene[0])
		else:
			if gene[0] in degGenes:
				degGenesNeg.append(gene[0])
			if gene[0] in cosmicGenes:
				cosmicGenesNeg.append(gene[0])
			if gene[0] in snvGenes:
				snvGenesNeg.append(gene[0])
				
	#Do some intersect things here
	
	allCriteriaIntersect = list(set(degGenesPos) & set(cosmicGenesPos) & set(snvGenesPos))
	cosmicSNVsIntersect = list(set(cosmicGenesPos) & set(snvGenesPos))
	cosmicDEGsIntersect = list(set(cosmicGenesPos) & set(degGenesPos))
	snvDEGsIntersect = list(set(snvGenesPos) & set(degGenesPos))
	# print "Number of genes that are in COSMIC, have SNVs and are DEG: ", len(allCriteriaIntersect)
	# print "Number of genes that are in COSMIC and have SNVs: ", len(cosmicSNVsIntersect)
	# print "Number of genes that are in COSMIC and are DEG: ", len(cosmicDEGsIntersect)
	# print "Number of genes that have SNV and are DEG: ", len(snvDEGsIntersect)

	#Return in the format of the venn diagram method
	return len(cosmicGenesPos), len(snvGenesPos), len(cosmicSNVsIntersect), len(degGenesPos), len(cosmicDEGsIntersect), len(snvDEGsIntersect), len(allCriteriaIntersect)

#Get the gene scores from the file without permutations
realGeneScores = np.loadtxt(permutationDataFolder + "/" + rankedGenesFile, dtype="object")


#4. Compute the values at each threshold and add them to the plot
print "computing intersects for real data"
# realScoreCountsCosmic = dict() #keep per threshold what the scores are
# realScoreCountsSNVs = dict()
# realScoreCountsDEGS = dict()
# realScoreCountsCosmicSNVs = dict()
# realScoreCountsCosmicDEGs = dict()
# realScoreCountsSNVDEGs = dict()
# realScoreCountsAll = dict()

#output this to a new file
if not os.path.exists('ThresholdEnrichment'):
    os.makedirs('ThresholdEnrichment')
outFile = "ThresholdEnrichment/" + rankedGenesFile
with open(outFile, 'w') as outF:
	for threshold in range(0, maxScore):
		[cosmic, snvs, cosmicSNVs, degs, cosmicDEGs, snvDEGs, allCriteria] = computeCategoryMatches(realGeneScores, threshold)
		# realScoreCountsCosmic[threshold] = cosmic
		# realScoreCountsSNVs[threshold] = snvs
		# realScoreCountsDEGS[threshold] = degs
		# realScoreCountsCosmicSNVs[threshold] = cosmicSNVs
		# realScoreCountsCosmicDEGs[threshold] = cosmicDEGs
		# realScoreCountsSNVDEGs[threshold] = snvDEGs
		# realScoreCountsAll[threshold] = allCriteria
		# 
		outF.write(str(threshold) + "\t" + str(cosmic) + "\t" + str(snvs) + "\t" + str(cosmicSNVs) + "\t" + str(degs) + "\t" + str(cosmicDEGs) + "\t" + str(snvDEGs) + "\t" + str(allCriteria) + "\n")
	
	#realScoreCounts[threshold] = [[cosmic, snvs, cosmicSNVs, degs, cosmicDEGs, snvDEGs, allCriteria]]
	# ax.plot(threshold, cosmic, 'bo')
	# ax.plot(threshold, snvs, 'ko')
	# ax.plot(threshold, cosmicSNVs, 'ro')
	# ax.plot(threshold, degs, 'yo')
	# ax.plot(threshold, cosmicDEGs, 'co')
	# ax.plot(threshold, snvDEGs, 'mo')
	# ax.plot(threshold, allCriteria, 'go')
	

	
# 
# #Read the permuted scores, and add the values to the threshold plot (and store with pkl somewhere)
# permutationDataFolder = sys.argv[5] #provide the folder which contains the output of all permutation runs.
# print "compute the intersection counts for permutation data"
# geneScoreFiles = [f for f in listdir(permutationDataFolder) if isfile(join(permutationDataFolder, f))]
# permutedScoreCountsCosmic = dict() #keep per threshold what the scores are
# permutedScoreCountsSNVs = dict()
# permutedScoreCountsDEGS = dict()
# permutedScoreCountsCosmicSNVs = dict()
# permutedScoreCountsCosmicDEGs = dict()
# permutedScoreCountsSNVDEGs = dict()
# permutedScoreCountsAll = dict()
# for threshold in range(0, maxScore): #prepare dictionaries
# 	permutedScoreCountsCosmic[threshold] = []
# 	permutedScoreCountsSNVs[threshold] = []
# 	permutedScoreCountsDEGS[threshold] = []
# 	permutedScoreCountsCosmicSNVs[threshold] = []
# 	permutedScoreCountsCosmicDEGs[threshold] = []
# 	permutedScoreCountsSNVDEGs[threshold] = []
# 	permutedScoreCountsAll[threshold] = []
# 	
# for geneScoreFile in geneScoreFiles:
# 	print "at file num: ", len(permutedScoreCountsAll[threshold])
# 	if geneScoreFile == "realSVs_geneScores.txt":
# 		continue
# 	
# 	#Read the file
# 	geneScores = np.loadtxt(permutationDataFolder + "/" + geneScoreFile, dtype="object")
# 	
# 	#Go through the same thresholds
# 	for threshold in range(0, maxScore):
# 		[cosmic, snvs, cosmicSNVs, degs, cosmicDEGs, snvDEGs, allCriteria] = computeCategoryMatches(geneScores, threshold)
# 
# 		permutedScoreCountsCosmic[threshold].append(cosmic)
# 		permutedScoreCountsSNVs[threshold].append(snvs)
# 		permutedScoreCountsDEGS[threshold].append(degs)
# 		permutedScoreCountsCosmicSNVs[threshold].append(cosmicSNVs)
# 		permutedScoreCountsCosmicDEGs[threshold].append(cosmicDEGs)
# 		permutedScoreCountsSNVDEGs[threshold].append(snvDEGs)
# 		permutedScoreCountsAll[threshold].append(allCriteria)
# 		#if threshold not in permutedScoreCounts:
# 		#	permutedScoreCounts[threshold] = []	
# 		#permutedScoreCounts[threshold].append([cosmic, snvs, cosmicSNVs, degs, cosmicDEGs, snvDEGs, allCriteria])
