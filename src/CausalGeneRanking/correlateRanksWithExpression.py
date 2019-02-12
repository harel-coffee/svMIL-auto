"""
	Obtain the gene ranking scores and also the expression z-scores of the genes and see if there is a correlation.
	Do the genes with the highest ranking also show different expression patterns? 
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

rankingFile = sys.argv[1]
expressionFile = sys.argv[2]
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

geneExpression = dict()
diffExprCount = 0
geneCount = 0
with open(expressionFile, 'r') as exprF:
	for line in exprF:
		line = line.strip()
		
		splitLine = line.split("\t")
		geneName = splitLine[0].replace(" ", "")
		zScore = splitLine[1].replace(" ", "")

		geneExpression[geneName] = float(zScore) #does the absolute help for correlation?
	

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

scoresAndExpression = []
for gene in geneScores:

	geneScore = geneScores[gene]
	if gene not in geneExpression: #many of the RNA genes won't be in the list. 
		continue

	
	expression = geneExpression[gene]
	
	if geneScore > 300:
		if gene in cosmicGenes:
			print "gene in COSMIC:"
		print gene
		print expression
		print geneScore
	
	scoresAndExpression.append([geneScore, expression])

exit()

scoresAndExpression = np.array(scoresAndExpression)
print scoresAndExpression

plt.scatter(scoresAndExpression[:,0], scoresAndExpression[:,1])
plt.xlabel("Gene ranking scores")
plt.ylabel("Log fold change")
plt.show()




