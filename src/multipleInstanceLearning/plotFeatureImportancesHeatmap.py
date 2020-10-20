## plot a heatmap with the top 100 instance features across all SV types.

import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import random
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
import os

outDir = sys.argv[1] #fix this if we run for each cancer type


def generateHeatmap(svTypes):
	
	#this part should later be across all cancer types

	#first get the instances and their ranking across all SV types
	importances, instances = getFeatureImportances(svTypes)
	#then compute the significances of the top 100 to random 100 instances
	pValues = computeFeatureSignificances(importances, instances, 100)

	#then make a heatmap plot of the significances.
	
def computeFeatureSignificances(importances, instances, top):

	#rank the importances by score
	indices = np.argsort(importances)[::-1]

	#then get the top instances
	topInstances = instances[indices[0:top]]

	#compute the percentages in these top X instances
	avgInstances = np.sum(topInstances, axis=0)
	totalInstances = avgInstances / topInstances.shape[0]

	#then compare to 100 random instances to see if it is significant.
	#per feature, have a distribution
	nullDistributions = dict()
	for i in range(0,top):

		if i == 0:
			for featureInd in range(0, len(totalInstances)):
				nullDistributions[featureInd] = []

		#sample as much random instances as in our filtered instances
		randomIndices = random.sample(range(0,instances.shape[0]), topInstances.shape[0])

		randomTopInstances = instances[randomIndices]

		#compute the percentages in these top X instances
		avgRandomInstances = np.sum(randomTopInstances, axis=0)

		totalRandomInstances = avgRandomInstances / randomTopInstances.shape[0]

		for featureInd in range(0, len(totalRandomInstances)):
			nullDistributions[featureInd].append(totalRandomInstances[featureInd])

	#for each feature, compute a z-score
	featurePValues = []
	featureZScores = []
	for featureInd in range(0, len(nullDistributions)):

		if np.std(nullDistributions[featureInd]) == 0:
			z = 0
			pValue = 1
			featureZScores.append(z)
			featurePValues.append(pValue)
			continue

		z = (totalInstances[featureInd] - np.mean(nullDistributions[featureInd])) / float(np.std(nullDistributions[featureInd]))
		pValue = stats.norm.sf(abs(z))*2

		featureZScores.append(z)
		featurePValues.append(pValue)

	#do MTC on the p-values

	reject, pAdjusted, _, _ = multipletests(featurePValues, method='bonferroni')

	#for cosmic, remove the cosmic feature
	#del featureZScores[31]
	#pAdjusted = np.delete(pAdjusted, 31)

	print(pAdjusted)
	return pAdjusted

def getFeatureImportances(svTypes):
	
	#gather the top 100 instances across all SV types
	#also return the instances themselves to get the features
	allInstances = []
	allImportances = []

	for svType in svTypes:
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

		#load the similarity matrix of this SV type
		dataPath = outDir + '/multipleInstanceLearning/similarityMatrices/'
		similarityMatrix = np.load(dataPath + '/similarityMatrix_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagLabels = np.load(dataPath + '/bagLabelsSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		instances = np.load(dataPath + '/instancesSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagPairLabels = np.load(dataPath + '/bagPairLabelsSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagMap = np.load(dataPath + '/bagMap_' + svType + '.npy', encoding='latin1', allow_pickle=True).item()
		filteredFeatures = np.loadtxt(dataPath + '/lowVarianceIdx_' + svType + '.txt')

		#train the classifier on the full dataset
		clf.fit(similarityMatrix, bagLabels)
		print('Classifier performance on all data: ', clf.score(similarityMatrix, bagLabels))
		#get the feature importances
		importances = list(clf.feature_importances_)
		
		allImportances += importances
		
		#because low index features are removed, add them back here if necessary
		#to retain an overview of all used features
		fixedInstances = []
		for instance in instances:
			finalLength = len(instance) + len(filteredFeatures)
			instanceMal = np.zeros(finalLength)
			addedFeatures = 0
			for featureInd in range(0, finalLength):

				if featureInd in filteredFeatures:
					instanceMal[featureInd] = 0
					addedFeatures += 1
				else:
					instanceMal[featureInd] = instance[featureInd-addedFeatures]

			fixedInstances.append(instanceMal)

		allInstances += fixedInstances

	allImportances = np.array(allImportances)
	allInstances = np.array(allInstances)
	return allImportances, allInstances

generateHeatmap(['DEL', 'DUP', 'INV', 'ITX'])

		