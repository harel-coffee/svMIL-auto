"""
	The goal of this script is to generate the similarity matrices for MIL

	If we do feature selection, we have similarity matrices for each feature
	output these as well.

"""

import sys
import os
import numpy as np
import pickle as pkl
import random

featureElimination = sys.argv[2]
leaveOnePatientOut = sys.argv[3] #make the similarity matrices for each left out patient
svTypes = ['DEL', 'DUP', 'INV', 'ITX']

outDir = sys.argv[1]
finalOutDir = outDir + '/multipleInstanceLearning/similarityMatrices/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

if featureElimination == "True":
	featureEliminationOutDir = finalOutDir + '/featureSelection'
	if not os.path.exists(featureEliminationOutDir):
		os.makedirs(featureEliminationOutDir)

if leaveOnePatientOut == 'True':
	leaveOnePatientOutDir = finalOutDir + '/leaveOnePatientOut'
	if not os.path.exists(leaveOnePatientOutDir):
		os.makedirs(leaveOnePatientOutDir)

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

	return similarityMatrix

def getSimilarityMatrixTest(testBags, trainInstances):
	bagIndices = np.arange(testBags.shape[0])
	similarityMatrix = np.zeros([testBags.shape[0], trainInstances.shape[0]])

	for bagInd in range(0, testBags.shape[0]):

		#get the average of all instances in this test patient bag
		testInstances = testBags[bagInd]

		instanceAvg = np.mean(testInstances, axis=0)

		#compute distance to all other instances from this bag average
		distance = np.abs(instanceAvg - trainInstances)

		#sum the distances to get 1 similarity score
		summedDistance = np.sum(distance,axis=1)
		similarityMatrix[bagInd,:] = summedDistance

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
			#if float(degPairInfo[5]) > 1 or float(degPairInfo[5]) < -1:
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

	if positiveBags.shape[0] == 0 or negativeBags.shape[0] == 0:
		continue

	#set a random seed to always subsample the same set
	if leaveOnePatientOut == 'False':
		np.random.seed(0)
		#subsample the negative set to the same number of positives.
		negativeBagsSubsampled = np.random.choice(negativeBags, positiveBags.shape[0])

		negativeBagsSubsampleInd = np.random.choice(np.arange(negativeBags.shape[0]), positiveBags.shape[0])
		negativeBagsSubsampled = negativeBags[negativeBagsSubsampleInd]

		negativeBagPairNamesSubsampled = negativeBagPairNames[negativeBagsSubsampleInd]
		bagPairLabels = np.concatenate((positiveBagPairNames, negativeBagPairNamesSubsampled))

		#save the bag pair labels for later
		np.save(finalOutDir + '/bagPairLabels_' + svType + '.npy', bagPairLabels)

		#merge the bags so that we can easily get to 1 similarity matrix and do all-to-all computations
		bags = np.concatenate((positiveBags, negativeBagsSubsampled))
		#assign bag labels
		bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBagsSubsampled.shape[0])
	else:
		#in case of leave-one-patient out, we subsample later on
		bagPairLabels = np.concatenate((positiveBagPairNames, negativeBagPairNames))

		#save the bag pair labels for later
		np.save(finalOutDir + '/bagPairLabelsNotSubsampled_' + svType + '.npy', bagPairLabels)

		#merge the bags so that we can easily get to 1 similarity matrix and do all-to-all computations
		bags = np.concatenate((positiveBags, negativeBags))
		#assign bag labels
		bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBags.shape[0])


	#output the bag labels which we can later read with the matrices
	np.save(finalOutDir + '/bagLabels_' + svType + '.npy', bagLabels)

	#stack the instances in the bags so that we can easily compute bag-instance distances
	instances = np.vstack(bags)

	#also output the instances for later
	np.save(finalOutDir + '/instances_' + svType + '.npy', instances)

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

	#save bagmap for later
	np.save(finalOutDir + '/bagMap_' + svType + '.npy', bagMap)

	#if we do feature selection, randomize the features here
	featureCount = instances.shape[1]

	featureStart = featureCount-1
	if featureElimination == "True":
		featureStart = 0 #set this to featureCount to run with all features. (make setting later)

	#if featureStart is not updated, this will run once
	#otherwise it will randomize a new feature each time
	for featureInd in range(featureStart, featureCount):
		print('current feature: ', featureInd+1)

		if featureElimination == "True":

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
			np.save(featureEliminationOutDir + '/similarityMatrix_' + svType + '_' + str(featureInd) + '.npy', similarityMatrix)
			
		elif featureElimination == 'False' and leaveOnePatientOut == 'False':
			#Make similarity matrix
			print("generating similarity matrix")
			similarityMatrix = getSimilarityMatrix(bags, instances, reverseBagMap)
			np.save(finalOutDir + '/similarityMatrix_' + svType + '.npy', similarityMatrix)
			print(finalOutDir + '/similarityMatrix_' + svType + '.npy')
		elif featureElimination == 'False' and leaveOnePatientOut == 'True':
			#go through the bags, and get all bags of all but one patient.
			#make a similarity matrix based on this.
			#so the output should be 2 matrices, one for training, one for testing
			#also make sure to output the labels separately.

			#first, get the bags and labels per patient
			perPatientPositiveBags = dict()
			for bagInd in range(0, positiveBags.shape[0]):

				#get the label of this bag
				bagPairLabel = positiveBagPairNames[bagInd]
				splitLabel = bagPairLabel.split('_')

				patientId = splitLabel[7]
				if patientId not in perPatientPositiveBags:
					perPatientPositiveBags[patientId] = dict()
					perPatientPositiveBags[patientId]['bags'] = []

				perPatientPositiveBags[patientId]['bags'].append(positiveBags[bagInd])

			perPatientNegativeBags = dict()
			for bagInd in range(0, negativeBags.shape[0]):

				#get the label of this bag
				bagPairLabel = negativeBagPairNames[bagInd]
				splitLabel = bagPairLabel.split('_')

				patientId = splitLabel[7]
				if patientId not in perPatientNegativeBags:
					perPatientNegativeBags[patientId] = dict()
					perPatientNegativeBags[patientId]['bags'] = []

				perPatientNegativeBags[patientId]['bags'].append(negativeBags[bagInd])

			#for each patient, randomly subsample as many negative bags as there are positives
			perPatientBags = dict()
			skippedPatients = 0
			for patient in perPatientPositiveBags:

				if patient not in perPatientNegativeBags:
					skippedPatients += 1
					continue

				if patient not in perPatientBags:
					perPatientBags[patient] = dict()
					perPatientBags[patient]['bags'] = []
					perPatientBags[patient]['labels'] = []

				patientNegativeBags = perPatientNegativeBags[patient]['bags']
				patientNegativeBags = np.array(patientNegativeBags)
				#add the same number of positives/negatives

				if len(perPatientPositiveBags[patient]['bags']) > patientNegativeBags.shape[0]:
					sampleCount = patientNegativeBags.shape[0]

					randomInd = random.sample(range(0, patientNegativeBags.shape[0]), sampleCount)

					randomNegativeBags = patientNegativeBags[randomInd]

					for bag in randomNegativeBags:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(0)

					for bag in perPatientPositiveBags[patient]['bags']:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(1)

				else:
					sampleCount = len(perPatientPositiveBags[patient]['bags'])

					randomInd = random.sample(range(0, patientNegativeBags.shape[0]), sampleCount)

					randomNegativeBags = patientNegativeBags[randomInd]

					for bag in randomNegativeBags:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(0)

					for bag in perPatientPositiveBags[patient]['bags']:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(1)
				
			print(skippedPatients)

			#go through each patient, and divide into train/test
			#the training set will be a merge of the bags of all other patients.
			#then for each patient, get the train/test combination
			foldSize = 10

			import math
			folds = math.ceil(len(perPatientBags) / foldSize)
			print(folds)

			testPatients = dict()
			trainPatients = dict()
			ind = 1
			foldInd = 0
			for patient in perPatientBags:
				if foldInd not in testPatients:
					testPatients[foldInd] = []
					trainPatients[foldInd] = []

				testPatients[foldInd].append(patient)
				if ind % 10 == 0:

					for patient2 in perPatientBags:

						if patient2 not in testPatients[foldInd]:
							trainPatients[foldInd].append(patient2)
					foldInd += 1

				ind += 1

			for fold in testPatients:

				testBags = []
				trainBags = []
				testLabels = []
				trainLabels = []
				for patient in perPatientBags:

					patientBags = perPatientBags[patient]['bags']
					patientLabels = perPatientBags[patient]['labels']

					patientBags = np.array(patientBags)
					patientLabels = np.array(patientLabels)



					otherPatientBags = []
					otherPatientLabels = []
					for patient2 in perPatientBags:

						if patient == patient2:
							continue
						#only if this patient is different, get those bags and labels
						otherPatientBags += perPatientBags[patient2]['bags']
						otherPatientLabels += perPatientBags[patient2]['labels']

					#mrge the other patient data together for training data
					otherPatientBags = np.array(otherPatientBags)
					otherPatientLabels = np.array(otherPatientLabels)

					#make a similarity matrix

					#get instances
					otherPatientInstances = np.vstack(otherPatientBags)

					#this needs a bag map, which is changed each time we make subsets.
					reverseBagMapOtherPatients = dict() #lookup instance by bag index
					instanceInd = 0
					for bagInd in range(0, otherPatientBags.shape[0]):
						reverseBagMapOtherPatients[bagInd] = []
						for instance in otherPatientBags[bagInd]:
							reverseBagMapOtherPatients[bagInd].append(instanceInd)

							instanceInd += 1

				#collect all this information as total bags/labels

				if ind % foldSize == 0:
					#make train/test

					similarityMatrixTrain = getSimilarityMatrix(otherPatientBags, otherPatientInstances, reverseBagMapOtherPatients)
					print(similarityMatrixTrain.shape)
					#now the curent patient bags need to be to the instances of the training set
					similarityMatrixTest = getSimilarityMatrixTest(patientBags, otherPatientInstances)
					print(similarityMatrixTest.shape)
						
					#write these data to disk so that we can access it later on
					np.save(leaveOnePatientOutDir + '/' + 'similarityMatrixTrain_' + patient + '_' + svType + '.npy', similarityMatrixTrain)
					np.save(leaveOnePatientOutDir + '/' + 'similarityMatrixTest_' + patient + '_' + svType + '.npy', similarityMatrixTest)
						
					#also save the labels
					np.save(leaveOnePatientOutDir + '/' + 'bagLabelsTrain_' + patient + '_' + svType + '.npy', otherPatientLabels)
					np.save(leaveOnePatientOutDir + '/' + 'bagLabelsTest_' + patient + '_' + svType + '.npy', patientLabels)

					ind += 1

		else:
			print('Combination of options not implemented')
			exit(1)
			