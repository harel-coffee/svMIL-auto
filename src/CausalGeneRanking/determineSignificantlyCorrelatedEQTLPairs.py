"""
	Script to determine which pairs of eQTLs are observed significantly more often in the same SV than by random chance. 

"""

import sys
import numpy as np
from glob import glob


countsDir = sys.argv[1]


#1. Determine the significance by randomly shuffling the SVs
genePairScores = dict()
pairLargerThanCounts = dict() #for each pair of genes, keep an index of how often the real scores are higher than the shuffled case. 

#Read the real counts
with open(countsDir + "/realSVs_counts.txt", 'r') as inF:
	
	for line in inF:
		splitLine = line.split("\t")
		
		genePair = splitLine[0]
		genePairScore = splitLine[1]
		
		genePairScores[genePair] = genePairScore
		pairLargerThanCounts[genePair] = 0
		
#Read the counts for the shuffled SVs.
#Determine for each gene pair how often the counts in the real case are higher than in the shuffled cases.

shuffledFiles = glob(countsDir + "/shuffled*")
#Collect all the files for the shuffled cases

for shuffledFile in shuffledFiles:
	
	with open(shuffledFile, 'r') as inF:
		
		for line in inF:
			splitLine = line.split("\t")
			
			genePair = splitLine[0]
			genePairScore = splitLine[1]
			
			if genePairScore >= genePairScores[genePair]:
				pairLargerThanCounts[genePair] += 1
			

#Compute the significance
significantPairs = []
for genePair in pairLargerThanCounts:
	
	significance = pairLargerThanCounts[genePair] / float(len(shuffledFiles))
	if significance < 0.05:
		significantPairs.append([genePair, pairLargerThanCounts[genePair], significance])
	
significantPairs = np.array(significantPairs)
print significantPairs.shape
print len(genePairScores)


