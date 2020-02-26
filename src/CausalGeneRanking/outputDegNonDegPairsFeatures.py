## Read the rule-based SV-gene pairs and split into deg & non-deg pairs

import sys
import numpy as np


svGenePairs = np.loadtxt(sys.argv[1], dtype='object')
#degPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')

degPairs = np.loadtxt(sys.argv[2], dtype='object')

positivePairsFeatures = []
negativePairsFeatures = []

for pair in svGenePairs:
	features = list(pair)
	
	#check for enhancer gain or loss. Otherwise, ignore.
	#f2 is enhancer loss.
	#f28 is enhancer gain
	#+1 because the first entry is the pair name. 
	
		
	if str(features[2]) != '1.0' and str(features[28]) != '1.0':
		continue 
	
	splitPair = pair[0].split('_')
	shortPair = splitPair[7] + '_' + splitPair[0]
	
	if shortPair in degPairs[:,0]:
		
		#check if this is true or not.
		degPairInfo = degPairs[degPairs[:,0] == shortPair][0]
		#features.append(np.abs(float(degPairInfo[5])))
			
		#if degPairInfo[3] == 'True':
		# 	positivePairsFeatures.append(features)
		# else:
		# 	negativePairsFeatures.append(features)

		if float(degPairInfo[5]) > 2 or float(degPairInfo[5]) < -2:
			positivePairsFeatures.append(features)
		else:
			negativePairsFeatures.append(features)


positivePairsFeatures = np.array(positivePairsFeatures)
negativePairsFeatures = np.array(negativePairsFeatures)

print(positivePairsFeatures.shape)
print(negativePairsFeatures.shape)

np.savetxt(sys.argv[1] + '_degPairsFeatures.txt', positivePairsFeatures, fmt='%s', delimiter='\t')
np.savetxt(sys.argv[1] + '_nonDegPairsFeatures.txt', negativePairsFeatures, fmt='%s', delimiter='\t')