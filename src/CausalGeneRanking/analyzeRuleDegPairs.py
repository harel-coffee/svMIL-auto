import numpy as np

ruleSVsDegPairs = np.load("ruleSvGenePairs.txt_degPairs.npy", allow_pickle=True, encoding='latin1')

geneCounts = dict()
for pair in ruleSVsDegPairs:
	splitPair = pair[0].split("_")
	if splitPair[0] not in geneCounts:
		geneCounts[splitPair[0]] = 0
	geneCounts[splitPair[0]] += 1
	
for gene in geneCounts:
	if geneCounts[gene] > 1:
		print("gene ", gene, " found ", geneCounts[gene], " times")