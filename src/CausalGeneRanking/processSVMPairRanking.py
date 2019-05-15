

#Load the pair ranking, and process that to a top X of genes

import sys
import numpy as np

pairRanking = np.loadtxt(sys.argv[1], dtype="object")

threshold = 300 # look at the top X pairs
topGenes = dict()
count = 0
for pair in pairRanking[:,0]:
	
	splitName = pair.split("_")
	gene = splitName[len(splitName)-1]
	if gene not in topGenes:
		topGenes[gene] = 0
	topGenes[gene] += 1
	
	count += 1
	if count > threshold:
		break
	
#output in numpy accesible format
topGenesArray = np.empty([len(topGenes),2], dtype="object")
for geneInd in range(0, len(topGenes)):
	gene = topGenes.keys()[geneInd]
	topGenesArray[geneInd, 0] = gene
	topGenesArray[geneInd, 1] = topGenes[gene]


topGenesSorted = topGenesArray[np.argsort(topGenesArray[:,1])][::-1]
print topGenesSorted
print topGenesSorted.shape

np.savetxt(sys.argv[2], topGenesSorted, delimiter='\t', fmt='%s')
