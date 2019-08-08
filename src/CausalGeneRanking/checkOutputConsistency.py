"""
	Check between 2 files if the genes have the same scores

"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
from six.moves import range

fileA = np.loadtxt(sys.argv[1], dtype="object")
fileB = np.loadtxt(sys.argv[2], dtype="object")

diffGeneCount = 0
diffGenes = []
for geneA in fileA:
	
	geneB = fileB[fileB[:,0] == geneA[0]][0]
	diff = False
	for feature in range(0, len(geneA)):
		if geneA[feature] != geneB[feature]:
			diff = True
			print("Genes are different for feature: ", geneA[feature], geneB[feature])
			print(geneA)
			print(geneB)
	
	
	if diff == True:
		diffGeneCount += 1
		diffGenes.append(geneA[0])

print(diffGenes)
print(diffGeneCount)
