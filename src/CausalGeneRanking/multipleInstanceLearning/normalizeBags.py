"""
	Normalize the bags prior to classification.

"""

import sys
import os
import numpy as np
import pickle as pkl

outDir = sys.argv[1]
finalOutDir = outDir + '/multipleInstanceLearning/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

#Get the bags
with open(outDir + '/linkedSVGenePairs/bags.pkl', 'rb') as handle:
	bagDict = pkl.load(handle)

instanceNum = len(bagDict[list(bagDict.keys())[0]][0])

#first check what the minimum and maximum feature values are in order to normalize
currentMax = [0]*instanceNum
currentMin = [float('inf')]*instanceNum
for pair in bagDict:
	
	for instance in bagDict[pair]:

		if instance[0] == '0' and instance[1] == '0': #skip instances where there are no gains or losses, if these slipped in somehow
			continue
		
		for featureInd in range(0, len(instance)):
			feature = instance[featureInd]
			if feature > currentMax[featureInd]:
				currentMax[featureInd] = feature
			if feature < currentMin[featureInd]:
				currentMin[featureInd] = feature

#loop through each bag, and normalize the features of each instance based on the max/min 
normalizedBagDict = dict()
for pair in bagDict:
	
	normalizedBagDict[pair] = []
	
	for instance in bagDict[pair]:
		
		normInstance = []
		
		for featureInd in range(0, len(instance)):
			
			
			feature = instance[featureInd]
			
			if currentMin[featureInd] == 0 and currentMax[featureInd] == 0: #if the min/max are 0 for this feature, the normalized value should also be 0. 
				normInstance.append(0)
				continue 
			
			#do the normalization
			normFeature = (feature-currentMin[featureInd])/(currentMax[featureInd]-currentMin[featureInd])
			normInstance.append(normFeature)
			
		normalizedBagDict[pair].append(normInstance)

bagDict = normalizedBagDict

#save to a file to prevent heavy computing load
with open(finalOutDir + '/normalizedBags.pkl', 'wb') as handle:
		pkl.dump(bagDict, handle, protocol=pkl.HIGHEST_PROTOCOL)