"""
	Use machine learning approaches to find the most informative features for detecting causal non-coding SVs. 

"""

import sys
import numpy as np

#Take windowed SV-pairs and DEG pairs as input
svGenePairs = np.loadtxt(sys.argv[1], dtype='object')
svGeneDegPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')

#Split the data into positive (all DEG pairs) and negative (pairs that are not in the DEG set)
positivePairs = svGeneDegPairs[:,0]
negativePairs = np.setdiff1d(svGenePairs, svGeneDegPairs[:,0])

#Randomly subsample the SV-pairs to be balanced with the DEG pairs
negativePairsSubsampled = np.random.choice(negativePairs, positivePairs.shape[0])

print positivePairs.shape
print negativePairsSubsampled.shape

#Assign features to each pair (gains/losses of elements yes/no)
#Load the rule-based SV-gene pairs, and get the gains & losses from there.
svGenePairsRules = np.loadtxt(sys.argv[3], dtype='object')

missing = 0
present = 0
for pair in svGenePairsRules[:,0]:
	if pair not in svGenePairs:
		missing += 1
	else:
		present += 1

print missing
print present
exit()


positivePairsFeatures = []
positiveWithFeatures = 0
positiveWithoutFeatures = 0
for pair in positivePairs:
	if pair in svGenePairsRules[:,0]:
		positiveWithFeatures += 1
		#get these features
		features = svGenePairsRules[svGenePairsRules[:,0] == pair,:][0]
		#skip the pair name and also the total score
		features = features[1:len(features)-1]
		features = [float(feature) for feature in features]
		positivePairsFeatures.append(features)
	else: #if not, assign all features as 0
		positiveWithoutFeatures += 1
		positivePairsFeatures.append([0]*26)

#repeat for negative pairs
negativePairsFeatures = []
negativeWithFeatures = 0
negativeWithoutFeatures = 0
for pair in negativePairsSubsampled:
	if pair in svGenePairsRules[:,0]:
		negativeWithFeatures += 1
		#get these features
		features = svGenePairsRules[svGenePairsRules[:,0] == pair,:][0]
		#skip the pair name and also the total score
		features = features[1:len(features)-1]
		features = [float(feature) for feature in features]
		negativePairsFeatures.append(features)
	else: #if not, assign all features as 0
		negativeWithoutFeatures += 1
		negativePairsFeatures.append([0]*26)

positivePairsFeatures = np.array(positivePairsFeatures)
negativePairsFeatures = np.array(negativePairsFeatures)

print positivePairsFeatures
print negativePairsFeatures

print positiveWithFeatures
print positiveWithoutFeatures
print negativeWithFeatures
print negativeWithoutFeatures
		
#First make a PCA and tSNE
