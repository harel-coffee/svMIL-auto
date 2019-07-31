
import sys

set1 = sys.argv[1]
set2 = sys.argv[2]


set1Genes = []
with open(set1, 'r') as inF:
	for line in inF:
		line = line.strip()
		splitLine = line.split(" ")
		
		geneName = splitLine[0].replace("'", "")
		geneName = geneName.replace("[", "")
		
		set1Genes.append(geneName)

set2Genes = []
with open(set2, 'r') as inF:
	for line in inF:
		line = line.strip()
		splitLine = line.split(" ")
		
		geneName = splitLine[0].replace("'", "")
		geneName = geneName.replace("[", "")
		
		set2Genes.append(geneName)
		
#Number of genes in set1 that are also in set2
import numpy as np
intersect = np.intersect1d(set1Genes,set2Genes)
print intersect
print len(intersect)

diff = np.setdiff1d(set1Genes,set2Genes)
print diff
