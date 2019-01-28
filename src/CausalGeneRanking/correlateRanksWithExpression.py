"""
	Obtain the gene ranking scores and also the expression z-scores of the genes and see if there is a correlation.
	Do the genes with the highest ranking also show different expression patterns? 
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

rankingFile = sys.argv[1]
expressionFile = sys.argv[2]

geneExpression = dict()
diffExprCount = 0
geneCount = 0
with open(expressionFile, 'r') as exprF:
	lineCount = 0
	for line in exprF:
		
		if lineCount < 1:
			lineCount += 1
			continue
		
		#First get the median z-score for the gene
		line = line.strip()
		splitLine = line.split("\t")
		
		hugoGeneName = splitLine[0]
		entrezGeneName = splitLine[1]
		
		#The rest is samples.
		sampleValues = []
		geneDiffCount = 0
		for colInd in range(2, len(splitLine)):
			#print splitLine[colInd]
			if splitLine[colInd] != "NA":
				sampleValues.append(float(splitLine[colInd]))
				if abs(float(splitLine[colInd])) >= 1.96:
					geneDiffCount += 1
		if geneDiffCount > 0:
			diffExprCount += 1
		geneCount += 1
		#Take the median for now across the samples
		medianZScore = np.median(sampleValues)
		geneExpression[hugoGeneName] = np.abs(medianZScore) #does the absolute help for correlation?

print diffExprCount
print geneCount
exit()


#Read the ranking file
geneScores = dict()
with open(rankingFile, 'r') as rankF:
	
	for line in rankF:
		line = line.strip()
		splitLine = line.split()
		
		geneName = splitLine[0]
		geneScore = float(splitLine[1])
		
		geneScores[geneName] = geneScore
		
#Plot correlation
#Convert dictionaries to same order of genes

scoresAndExpression = []
for gene in geneScores:
	
	geneScore = geneScores[gene]
	if gene not in geneExpression: #many of the RNA genes won't be in the list. 
		continue
	
	
	expression = geneExpression[gene]
	
	scoresAndExpression.append([geneScore, expression])



scoresAndExpression = np.array(scoresAndExpression)

plt.scatter(scoresAndExpression[:,0], scoresAndExpression[:,1])
plt.xlabel("Gene ranking scores")
plt.ylabel("Absolute median expression z-score")
plt.show()




