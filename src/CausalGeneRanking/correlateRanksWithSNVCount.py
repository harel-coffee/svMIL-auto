"""
	Obtain the gene ranking scores and also the number of SNVs of the genes and see if there is a correlation.
	
"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

rankingFile = sys.argv[1]
snvFile = sys.argv[2]
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

geneSnvCount = dict()
with open(snvFile, 'r') as snvF:
	for line in snvF:
		line = line.strip()
		
		splitLine = line.split("\t")
		geneName = splitLine[0]
		snvCount = splitLine[1]

		geneSnvCount[geneName] = int(snvCount) 
	

#Read the ranking file
geneScores = dict()
with open(rankingFile, 'r') as rankF:
	
	for line in rankF:
		line = line.strip()
		splitLine = line.split()
		
		geneName = splitLine[0]
		geneScore = float(splitLine[1]) + float(splitLine[4])
		
		geneScores[geneName] = geneScore
	
#Plot correlation
#Convert dictionaries to same order of genes

scoresAndSnvCounts = []
for gene in geneScores:

	geneScore = geneScores[gene]
	if gene not in geneSnvCount: #many of the RNA genes won't be in the list. 
		continue

	
	snvCount = geneSnvCount[gene]
	
	if snvCount > 2000:
		print(gene)
		print(snvCount)
		print(geneScore)
		if gene in cosmicGenes:
			print("gene in COSMIC")
	
	scoresAndSnvCounts.append([geneScore, snvCount])

exit()

scoresAndSnvCounts = np.array(scoresAndSnvCounts)
print(scoresAndSnvCounts)

plt.scatter(scoresAndSnvCounts[:,0], scoresAndSnvCounts[:,1])
plt.xlabel("Gene ranking scores")
plt.ylabel("Number of SNVs")
plt.show()




