"""
	Compute for each gene_SV pair how often a gene is significant

"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np

pValues = np.loadtxt(sys.argv[1], dtype="object")

geneCounts = dict()

for pair in pValues[:,0]:
	
	splitPair = pair.split("_")
	gene = splitPair[0]
	if gene not in geneCounts:
		geneCounts[gene] = 0
	geneCounts[gene] += 1
	

geneCountsArray = np.empty([len(geneCounts), 2], dtype="object")
geneCountsArray[:,0] = list(geneCounts.keys())
geneCountsArray[:,1] = list(geneCounts.values())

geneCountsArray[:,1] = geneCountsArray[:,1].astype(float)

geneCountsArray = geneCountsArray[geneCountsArray[:,1].argsort()][::-1]

print(geneCountsArray)
print(geneCountsArray.shape)
print(len(np.where(geneCountsArray[:,1] > 1)[0]))

