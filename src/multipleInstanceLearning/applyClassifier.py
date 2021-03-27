"""
	The goal of this script is to apply MIL trained on HMF data to another
	dataset and predict the positive SV-gene pairs.

"""

import sys
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix
import os.path

outDirTrain = sys.argv[1]
outDirTest = sys.argv[2]

def getSimilarityMatrixTest(testBags, trainInstances, labels):
	"""
		function to get the similarity matrix specific for the test case.
		The instances that we map the distance to are the provided train instances.

		testBags (numpy array): all test bags that we use for this matrix
		trainInstances (numpy array): all instances in the bags of the training data, we compute distance to these instances from the test bags
		labels (list): obsolete.
	"""

	similarityMatrix = np.zeros([testBags.shape[0], trainInstances.shape[0]])

	#print(similarityMatrix.shape)

	for bagInd in range(0, testBags.shape[0]):
		#print(labels[bagInd])
		#get the average of all instances in this test patient bag
		testInstances = testBags[bagInd]

		instanceAvg = np.mean(testInstances, axis=0)

		#compute distance to all other instances from this bag average
		distance = np.abs(instanceAvg - trainInstances)

		#sum the distances to get 1 similarity score
		summedDistance = np.sum(distance,axis=1)
		#print(summedDistance)
		similarityMatrix[bagInd,:] = summedDistance

	return similarityMatrix


#1. Get the training similarity matrix on the full training data
svTypes = ['DEL', 'DUP', 'INV', 'ITX']
for svType in svTypes:
	print('sv type: ', svType)
	#define the classifiers to use (from optimization)
	#would be nicer if these are in 1 file somewhere, since they are also used in another script

	if svType == 'DEL':
		clf = RandomForestClassifier(random_state=785, n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'deletions'
	elif svType == 'DUP':
		clf = RandomForestClassifier(random_state=785, n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'duplications'
	elif svType == 'INV':
		clf = RandomForestClassifier(random_state=785, n_estimators= 200, min_samples_split=5, min_samples_leaf=4, max_features='auto', max_depth=10, bootstrap=True)
		title = 'inversions'
	elif svType == 'ITX':
		clf = RandomForestClassifier(random_state=785, n_estimators= 1000, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'translocations'
	else:
		print('SV type not supported')
		exit(1)

	dataPath = outDirTrain + '/multipleInstanceLearning/similarityMatrices/'
	
	if not os.path.exists(dataPath + '/similarityMatrix_' + svType + '.npy'):
		continue
	
	similarityMatrix = np.load(dataPath + '/similarityMatrix_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	bagLabelsTrain = np.load(dataPath + '/bagLabelsSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	trainInstances = np.load(dataPath + '/instancesSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)

	#train the classifier on the full dataset
	clf.fit(similarityMatrix, bagLabelsTrain)
	print('Classifier performance on all data: ', clf.score(similarityMatrix, bagLabelsTrain))

	#load the test data
	dataPath = outDirTest + '/multipleInstanceLearning/similarityMatrices/'
	#load the bags
	bags = np.load(dataPath + '/bagsNotSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	bagLabelsTest = np.load(dataPath + '/bagLabelsNotSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	instances = np.load(dataPath + '/instancesNotSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	bagPairLabels = np.load(dataPath + '/bagPairLabelsNotSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)

	#generate test similarity matrix (to the training dataset)
	similarityMatrixTest = getSimilarityMatrixTest(bags, trainInstances, [])

	#test classifier.
	print('Classifier performance on test data: ', clf.score(similarityMatrixTest, bagLabelsTest))
	proba = clf.predict_proba(similarityMatrixTest)
	predictions = clf.predict(similarityMatrixTest)

	average_precision = average_precision_score(bagLabelsTest, proba[:,1])
	print('Average precision: ', average_precision)
	
	#Get the probabilities of the SV-gene pairs, this is a ranked list/prioritization!
	rankedList = []
	for pairLabelInd in range(0, len(bagPairLabels)):
		rankedList.append([bagPairLabels[pairLabelInd], proba[:,1][pairLabelInd]])

	rankedList = np.array(rankedList, dtype='object')
	sortedRankedList = rankedList[np.argsort(rankedList[:,1])][::-1]

	print('Prioritized SV-gene pairs: ')
	print(sortedRankedList)





