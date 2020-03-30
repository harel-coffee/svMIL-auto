"""
	The goal of this script is to generate the similarity matrices for MIL

	If we do feature selection, we have similarity matrices for each feature
	output these as well.

"""

import sys
import os
import numpy as np
import pickle as pkl

featureElimination = False
svTypes = ['DEL', 'DUP', 'INV', 'ITX']
svTypes = ['ITX']

outDir = sys.argv[1]
finalOutDir = outDir + '/multipleInstanceLearning/similarityMatrices/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

if featureElimination == True:
	featureEliminationOutDir = finalOutDir + '/featureSelection'
	if not os.path.exists(featureEliminationOutDir):
		os.makedirs(featureEliminationOutDir)


#input the normalized bags
with open(outDir + '/multipleInstanceLearning/normalizedBags.pkl', 'rb') as handle:
	bagDict = pkl.load(handle)

#get the information for the bag labels
degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object') #labels

print(degPairs)

print("initial number of bags: ", len(bagDict))
print('initial deg pairs: ', degPairs.shape[0])

#function to get the similarity matrix
def getSimilarityMatrix(bags, instances, reverseBagMap):

	bagIndices = np.arange(bags.shape[0])
	similarityMatrix = np.zeros([bags.shape[0], instances.shape[0]])
	print("Number of bags: ", bags.shape[0])
	for bagInd in range(0, bags.shape[0]):

		#Get the indices of the instances that are in this bag
		instanceIndices = reverseBagMap[bagInd]
		instanceSubset = instances[instanceIndices,:]

		#get the average of all instances in this bag
		instanceAvg = np.mean(instanceSubset, axis=0)

		#compute distance to all other instances from this bag average
		distance = np.abs(instanceAvg - instances)

		#sum the distances to get 1 similarity score
		summedDistance = np.sum(distance,axis=1)
		similarityMatrix[bagInd,:] = summedDistance

	print('similarity matrix: ')
	print(similarityMatrix)
	return similarityMatrix

#then, generate the similarity matrices for the SV types
for svType in svTypes:
	
	#allow for running with feature selection
	bagLabels = []
	positiveBagPairNames = []
	negativeBagPairNames = []
	positiveInstanceLabels = []
	positiveBags = []
	negativeBags = []

	#for each SV-gene pair, get the instances
	for pair in bagDict:

		#check if the SV type matches our selection
		splitPair = pair.split("_")
		shortPair = splitPair[7] + '_' + splitPair[0]

		if svType != '':
			if splitPair[12] != svType:
				continue

		#get the label of the bag by checking if it exists in degPairs, some pairs do not have a z-score because the gene is excluded due to mutations.
		if shortPair in degPairs[:,0]:

			#get the z-score of the pair.
			degPairInfo = degPairs[degPairs[:,0] == shortPair][0]

			#if the z-score matches this criterion, the SV-gene pair is positive
			if float(degPairInfo[5]) > 1.5 or float(degPairInfo[5]) < -1.5:

				#go through the instances of this SV-gene pair, and include only those that have gains and losses, and more than 1 instance. This should in principle not happen, but good to keep a check.
				instances = []
				for instance in bagDict[pair]:

					if instance[0] == 0 and instance[1] == 0:
						continue
					positiveInstanceLabels.append(pair + '_' + '_'.join([str(i) for i in instance]))
					instances.append(instance)


				if len(instances) < 1:
					continue

				positiveBagPairNames.append(pair)
				positiveBags.append(instances)

			else: #if the z-score is anything else, this bag will be labeled negative.

				#get the right number of features per instance
				instances = []
				for instance in bagDict[pair]:

					if instance[0] == 0 and instance[1] == 0:
						continue

					instances.append(instance)

				if len(instances) < 1:
					continue

				negativeBags.append(instances)
				negativeBagPairNames.append(pair)


	positiveBags = np.array(positiveBags)
	negativeBags = np.array(negativeBags)
	positiveBagPairNames = np.array(positiveBagPairNames)
	negativeBagPairNames = np.array(negativeBagPairNames)

	print('Number of positive bags: ', positiveBags.shape)
	print('Number of negative bags: ', negativeBags.shape)

	print('Number of positive instances: ', len(positiveInstanceLabels))

	#set a random seed to always subsample the same set
	np.random.seed(0)
	#subsample the negative set to the same number of positives.
	negativeBagsSubsampled = np.random.choice(negativeBags, positiveBags.shape[0])

	negativeBagsSubsampleInd = np.random.choice(np.arange(negativeBags.shape[0]), positiveBags.shape[0])
	negativeBagsSubsampled = negativeBags[negativeBagsSubsampleInd]

	negativeBagPairNamesSubsampled = negativeBagPairNames[negativeBagsSubsampleInd]
	bagPairLabels = np.concatenate((positiveBagPairNames, negativeBagPairNamesSubsampled))

	#merge the bags so that we can easily get to 1 similarity matrix and do all-to-all computations
	bags = np.concatenate((positiveBags, negativeBagsSubsampled))
	#assign bag labels
	bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBagsSubsampled.shape[0])

	#stack the instances in the bags so that we can easily compute bag-instance distances
	instances = np.vstack(bags)
	print(instances[:,2])

	#Make an index where we can lookup at which position the instances are in the concatenated bag array.
	reverseBagMap = dict() #lookup instance by bag index
	bagMap = dict() #lookup bag by instance index
	instanceInd = 0
	for bagInd in range(0, bags.shape[0]):
		reverseBagMap[bagInd] = []
		for instance in bags[bagInd]:
			reverseBagMap[bagInd].append(instanceInd)
			bagMap[instanceInd] = bagInd

			instanceInd += 1

	#if we do feature selection, randomize the features here
	featureCount = instances.shape[1]

	featureStart = featureCount
	if featureElimination == True:
		featureStart = 0 #set this to featureCount to run with all features. (make setting later)

	#if featureStart is not updated, this will run once
	#otherwise it will randomize a new feature each time
	for featureInd in range(featureStart, featureCount+1):
		print('current feature: ', featureInd+1)

		if featureElimination == True:

			#randomize one feature across the bags
			#get all values of this instance
			shuffledInstanceValues = instances[:,featureInd]
			randomInd = np.arange(0, shuffledInstanceValues.shape[0])
			np.random.shuffle(randomInd)
			print(randomInd)

			#we compute the similarity matrix based on the instances
			#but the instance values need to be reset every iteration
			shuffledInstances = np.zeros(instances.shape)
			for col in range(0, instances.shape[1]):
				if col != featureInd:
					shuffledInstances[:,col] = instances[:,col]
				else:
					shuffledInstances[:,col] = instances[randomInd,col]

			print(shuffledInstances[:,featureInd])
			print(instances[:,featureInd])

			#bags are just used for shape, so no need to shuffle those
			similarityMatrix = getSimilarityMatrix(bags, shuffledInstances, reverseBagMap)
			#output this similarity matrix to a file.
			#output to a folder specific for the feature selection data
			np.save(featureEliminationOutDir + '/similarityMatrix_' + svType + '_' + featureInd + '.npy', similarityMatrix)

		else:
			#Make similarity matrix
			print("generating similarity matrix")
			similarityMatrix = getSimilarityMatrix(bags, instances, reverseBagMap)
			np.save(finalOutDir + '/similarityMatrix_' + svType + '.npy', similarityMatrix)
			print(finalOutDir + '/similarityMatrix_' + svType + '.npy')
			