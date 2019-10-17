## Read the rule-based SV-gene pairs and split into deg & non-deg pairs

import sys
import numpy as np


svGenePairs = np.loadtxt(sys.argv[1], dtype='object')
degPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')

positivePairsFeatures = []
negativePairsFeatures = []

for pair in svGenePairs:
	features = pair
	
	#check f there are really no gains
	gains = features[24:46]
	splitSV = features[0].split("_")
	if splitSV[8] == 'del':
		print(pair[0])
		print(gains)
	
	if pair[0] in degPairs[:,0]:
		positivePairsFeatures.append(features)
		
	else:
		negativePairsFeatures.append(features)

positivePairsFeatures = np.array(positivePairsFeatures)
negativePairsFeatures = np.array(negativePairsFeatures)

print(positivePairsFeatures.shape)
print(negativePairsFeatures.shape)

np.savetxt(sys.argv[1] + '_degPairsFeatures.txt', positivePairsFeatures, fmt='%s', delimiter='\t')
np.savetxt(sys.argv[1] + '_nonDegPairsFeatures.txt', negativePairsFeatures, fmt='%s', delimiter='\t')