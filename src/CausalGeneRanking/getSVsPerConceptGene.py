"""
	Go through the concept genes and the number of SVs linked in the ranking, output how many SVs are affecting these genes in total. 

"""

import sys
import numpy as np

conceptGenes = np.loadtxt(sys.argv[1], dtype="object")
rankedGenes = np.loadtxt(sys.argv[2], dtype="object")

totalCounts = []
for gene in conceptGenes[:,0]:
	
	geneRanking = rankedGenes[np.where(rankedGenes[:,0] == gene)[0]][0]
	patients = geneRanking[len(geneRanking)-1]
	
	splitPatients = patients.split(",")
	
	svCount = len(splitPatients)
	totalCounts.append(svCount)
	
print totalCounts
print np.sum(totalCounts)
	


