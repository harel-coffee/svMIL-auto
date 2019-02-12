"""
	Script to go through the list of ranked genes, and determine at which threshold of number of samples we find more COSMIC genes, genes with SNVs or DEGs than when SVs are randomized. 

"""

import sys
import matplotlib.pyplot as plt

rankedGenesFile = sys.argv[1]
cosmicGenesFile = sys.argv[2]
snvFile = sys.argv[3]
degFile = sys.argv[4]

#1. Get the ranked genes 
#Read the ranking file
geneScores = dict()
maxScore = 0 #Keep a max score that we use to determine the thresholds to search through
with open(rankedGenesFile, 'r') as rankF:
	
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
def computeCategoryMatches(rankedGenesFile, threshold):
	print "threshold: ", threshold
	degGenesPos = []
	degGenesNeg = []
	cosmicGenesPos = []
	cosmicGenesNeg = []
	snvGenesPos = []
	snvGenesNeg = []
	
	with open(rankedGenesFile, 'rb') as f:
		
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")
			
			if float(splitLine[1]) > threshold:
				if splitLine[0] in degGenes:
					degGenesPos.append(splitLine[0])
				if splitLine[0] in cosmicGenes:
					cosmicGenesPos.append(splitLine[0])
				if splitLine[0] in snvGenes:
					snvGenesPos.append(splitLine[0])
			else:
				if splitLine[0] in degGenes:
					degGenesNeg.append(splitLine[0])
				if splitLine[0] in cosmicGenes:
					cosmicGenesNeg.append(splitLine[0])
				if splitLine[0] in snvGenes:
					snvGenesNeg.append(splitLine[0])
					
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

ax = plt.subplot(1,1,1)
#4. Compute the values at each threshold and add them to the plot
for threshold in range(0, maxScore):
	[cosmic, snvs, cosmicSNVs, degs, cosmicDEGs, snvDEGs, allCriteria] = computeCategoryMatches(rankedGenesFile, threshold)
	
	ax.plot(threshold, cosmic, 'bo')
	ax.plot(threshold, snvs, 'ko')
	ax.plot(threshold, cosmicSNVs, 'ro')
	ax.plot(threshold, degs, 'yo')
	ax.plot(threshold, cosmicDEGs, 'co')
	ax.plot(threshold, snvDEGs, 'mo')
	ax.plot(threshold, allCriteria, 'go')
	#ax.plot(0, 1, 'blue')

plt.show()

#The y axis is the number in the overlap, the x axis the threshold
#All intersection values are individual points in the plot


