"""
	Correlate the frequency of genes found in the concept with various other things

"""

import numpy as np
import sys


#1. Correlate concept frequency with the total rank

#All scores
scores = np.loadtxt(sys.argv[1], dtype="object")

#concept genes and frequency
conceptGenes = np.loadtxt(sys.argv[2], dtype="object")

ranksAndTotalScores = np.zeros([conceptGenes.shape[0], 2])

for geneInd in range(0, conceptGenes.shape[0]):
	
	ranksAndTotalScores[geneInd, 0] = float(conceptGenes[geneInd, 1])
	
	geneScore = scores[np.where(scores[:,0] == conceptGenes[geneInd,0])][0]
	
	ranksAndTotalScores[geneInd, 1] = float(geneScore[30])


print np.corrcoef(ranksAndTotalScores[:,0], ranksAndTotalScores[:,1])

#Correlation with the number of samples


ranksAndTotalScores = np.zeros([conceptGenes.shape[0], 2])

for geneInd in range(0, conceptGenes.shape[0]):
	
	ranksAndTotalScores[geneInd, 0] = float(conceptGenes[geneInd, 1])
	
	geneScore = scores[np.where(scores[:,0] == conceptGenes[geneInd,0])][0]
	
	samples = geneScore[31].split(",")
	
	ranksAndTotalScores[geneInd, 1] = len(samples)


print np.corrcoef(ranksAndTotalScores[:,0], ranksAndTotalScores[:,1])
