"""
	The goal of this script is to apply MIL trained on HMF data to another
	dataset and predict the positive SV-gene pairs.

"""

import sys
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix

sys.path.insert(1, './multipleInstanceLearning')
#import generateSimilarityMatrices

outDir = sys.argv[1]
outDirPCAWG = sys.argv[2]

##### how to properly import this from other function??
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

#get the cosmic genes to check for cancer gene overrepresentation among the positives
cosmicGenes = []
causalGeneFile = '../data/genes/CCGC.tsv'
with open(causalGeneFile, 'r') as geneFile:
	
	lineCount = 0
	header = []
	for line in geneFile:
		splitLine = line.split("\t")
		#First extract the header and store it in the dictionary to remove dependency on the order of columns in the file
		if lineCount < 1:

			header = splitLine
			lineCount += 1
			continue
			
		#Obtain the gene name and gene position
		
		geneSymbolInd = header.index('Gene Symbol')
		genePositionInd = header.index('Genome Location')
		
		geneSymbol = splitLine[geneSymbolInd]
		genePosition = splitLine[genePositionInd]
		
		#Split the gene position into chr, start and end
		
		colonSplitPosition = genePosition.split(":")
		dashSplitPosition = colonSplitPosition[1].split("-")
		
		chromosome = colonSplitPosition[0]
		start = dashSplitPosition[0].replace('"',"") #apparently there are some problems with the data, sometimes there are random quotes in there
		end = dashSplitPosition[1].replace('"', "")
		
		cancerTypeInd = header.index('Tumour Types(Somatic)')
		cancerType = splitLine[cancerTypeInd]
		
		if start == '' or end == '':
			continue
		
		cosmicGenes.append([geneSymbol, cancerType])

cosmicGenes = np.array(cosmicGenes, dtype='object')		

#1. Get the training similarity matrix on the full HMF data

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
		
	#clf = RandomForestClassifier()

	dataPath = outDir + '/multipleInstanceLearning/similarityMatrices/'
	similarityMatrix = np.load(dataPath + '/similarityMatrix_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	bagLabelsTrain = np.load(dataPath + '/bagLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	trainInstances = np.load(dataPath + '/instances_' + svType + '.npy', encoding='latin1', allow_pickle=True)

	#train the classifier on the full dataset
	clf.fit(similarityMatrix, bagLabelsTrain)
	print('Classifier performance on all data: ', clf.score(similarityMatrix, bagLabelsTrain))

	#load the PCAWG OV data
	dataPath = outDirPCAWG + '/multipleInstanceLearning/similarityMatrices/'
	#load the bags
	bags = np.load(dataPath + '/bags_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	bagLabelsTest = np.load(dataPath + '/bagLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	instances = np.load(dataPath + '/instances_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	bagPairLabels = np.load(dataPath + '/bagPairLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)

	import random
	#random.shuffle(bagLabelsTest)

	#generate test similarity matrix
	similarityMatrixTest = getSimilarityMatrixTest(bags, trainInstances, [])
	
	#test classifier.
	print('Classifier performance on PCAWG: ', clf.score(similarityMatrixTest, bagLabelsTest))
	predictions = clf.predict(similarityMatrixTest)

	average_precision = average_precision_score(bagLabelsTest, predictions)
	print(average_precision)
	
	tn, fp, fn, tp = confusion_matrix(bagLabelsTest, predictions).ravel()
	
	print(tp, fp, tn, fn)
	
	for predictionInd in range(0, len(predictions)):
		
		prediction = predictions[predictionInd]
		
		if prediction == 1:
			
			label = bagPairLabels[predictionInd]
			gene = label.split('_')[0]
			
			if gene in cosmicGenes[:,0]:
				
				print(cosmicGenes[cosmicGenes[:,0] == gene])
	

