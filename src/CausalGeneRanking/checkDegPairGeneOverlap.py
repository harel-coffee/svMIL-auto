import sys
import numpy as np


svGenePairs = np.loadtxt(sys.argv[1], dtype='object')
degPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')

degData = []
nonDegData = []

for pair in svGenePairs:
	if pair[0] in degPairs[:,0]:
		degData.append(pair)
	else:
		nonDegData.append(pair)

degGenes = []
for pair in degData:
	splitPair = pair[0].split("_")
	degGenes.append(splitPair[0])

nonDegGenes = []
for pair in nonDegData:
	splitPair = pair[0].split("_")
	nonDegGenes.append(splitPair[0])

degGenes = np.unique(degGenes)
nonDegGenes = np.unique(nonDegGenes)
intersect = np.intersect1d(degGenes, nonDegGenes)
print(intersect.shape)
print(intersect)
print("no of deg genes: ", degGenes.shape)
print("no of non-deg genes: ", nonDegGenes.shape)


positivePairsFeatures = []
negativePairsFeatures = []

for pair in svGenePairs:
	features = pair[1:len(pair)]
	splitPair = pair[0].split("_")
	if pair[0] in degPairs[:,0] and splitPair[0] not in intersect:
		positivePairsFeatures.append(features)
	else:
		if splitPair[0] not in intersect:
			negativePairsFeatures.append(features)

positivePairsFeatures = np.array(positivePairsFeatures)
negativePairsFeatures = np.array(negativePairsFeatures)

print(positivePairsFeatures.shape)
print(negativePairsFeatures.shape)

np.savetxt(sys.argv[1] + '_diff_degPairsFeatures.txt', positivePairsFeatures, fmt='%s', delimiter='\t')
np.savetxt(sys.argv[1] + '_diff_nonDegPairsFeatures.txt', negativePairsFeatures, fmt='%s', delimiter='\t')