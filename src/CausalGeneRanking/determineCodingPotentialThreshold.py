"""

	For every coding SV-gene pair, get the SVs uniquely.
	For each gene associated with the SV, determine how many of these are in COSMIC, and how many of these genes are DEG. 

"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

codingPairs = np.loadtxt(sys.argv[1], dtype='object')
degPairs = np.load(sys.argv[1] + "_nonCodingPairDEGs.npy", allow_pickle=True)
cosmicGenesFile = sys.argv[3]

codingSVGenes = dict()
for pair in codingPairs[:,0]:
	
	splitPair = pair.split("_")
	svEntries = splitPair[1:]
	sv = "_".join(svEntries)
	
	if sv not in codingSVGenes:
		codingSVGenes[sv] = []
	codingSVGenes[sv].append(splitPair[0])
	
print(codingSVGenes)

#get the COSMIC genes

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

#Which of these genes are DEG?
#Which of these genes are COSMIC?
#Which of these genes are both DEG + COSMIC?

svDegCount = []
svCosmicCount = []
svDegAndCosmicCount = []
svDegOrCosmicCount = []
for sv in codingSVGenes:
	
	degCount = 0
	cosmicCount = 0
	degAndCosmicCount = 0
	degOrCosmicCount = 0
	
	for gene in codingSVGenes[sv]:
		pair = gene + "_" + sv
		
		if pair in degPairs[:,0]:
			degCount += 1

		if gene in cosmicGenes:
			cosmicCount += 1
			
		if gene in cosmicGenes and pair in degPairs[:,0]:
			degAndCosmicCount += 1
			
		if gene in cosmicGenes or pair in degPairs[:,0]:
			degOrCosmicCount += 1
	
	svDegCount.append(degCount)
	svCosmicCount.append(cosmicCount)
	svDegAndCosmicCount.append(degAndCosmicCount)
	svDegOrCosmicCount.append(degOrCosmicCount)

print(svDegCount)
print(svCosmicCount)
print(svDegAndCosmicCount)


svDegCount = np.array(svDegCount)
svCosmicCount = np.array(svCosmicCount)
svDegAndCosmicCount = np.array(svDegAndCosmicCount)
svDegOrCosmicCount = np.array(svDegOrCosmicCount)

print(len(np.where(svDegCount != 0)[0]))
print(len(np.where(svCosmicCount != 0)[0]))
print(len(np.where(svDegAndCosmicCount != 0)[0]))
print(len(np.where(svDegOrCosmicCount != 0)[0]))
print(len(codingSVGenes))

exit()

ind = np.argsort(svDegCount)
sortedDegCounts = svDegCount[ind]
sortedCosmicCounts = svCosmicCount[ind]
sortedDegAndCosmicCounts = svDegAndCosmicCount[ind]




plt.plot(sortedDegCounts)
plt.plot(sortedCosmicCounts)
plt.plot(sortedDegAndCosmicCounts)
plt.show()




