### test alternative to MIL

import sys
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn import model_selection
from sklearn.metrics import plot_roc_curve, auc, average_precision_score
import matplotlib.pyplot as plt
from scipy import interp
import random
from sklearn.ensemble import RandomForestClassifier

outDir = sys.argv[1]

positivePairs = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt__degPairsFeatures.txt', dtype='object')
negativePairs = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt__nonDegPairsFeatures.txt', dtype='object')

##Some features were not used in MIL, so we skip them here too to make it comparable.
#keep the features in the file to be able to plot them for figure 2.
allowedFeatures = list(np.arange(1,3)) #skip name (0) and promoter losses (3)
allowedFeatures += list(np.arange(4,29)) #skip promoter gains (29)
allowedFeatures += list(np.arange(30,69))

#normalize the data
#get the min and max from all pairs
currentMax = [0]*positivePairs.shape[1]
currentMin = [float('inf')]*positivePairs.shape[1]
for pair in positivePairs:

	for featureInd in range(1, len(pair)): #skip the name

		if float(pair[featureInd]) < currentMin[featureInd]:
			currentMin[featureInd] = float(pair[featureInd])
		if float(pair[featureInd]) > currentMax[featureInd]:
			currentMax[featureInd] = float(pair[featureInd])
for pair in negativePairs:

	for featureInd in range(1, len(pair)): #skip name

		if float(pair[featureInd]) < currentMin[featureInd]:
			currentMin[featureInd] = float(pair[featureInd])
		if float(pair[featureInd]) > currentMax[featureInd]:
			currentMax[featureInd] = float(pair[featureInd])

#the normalize the data
normalizedPositivePairs = []
for pair in positivePairs:

	normalizedValues = [pair[0]]
	for featureInd in range(1, len(pair)):
		if currentMin[featureInd] == 0 and currentMax[featureInd] == 0: #if the min/max are 0 for this feature, the normalized value should also be 0.
			normalizedValues.append(0)
			continue

		if currentMax[featureInd] - currentMin[featureInd] == 0:
			normalizedValues.append(0) #min and max should be the same, so no norm necessary
			continue

		normFeature = (float(pair[featureInd])-currentMin[featureInd])/(currentMax[featureInd]-currentMin[featureInd])
		normalizedValues.append(normFeature)

	normalizedPositivePairs.append(normalizedValues)

normalizedNegativePairs = []
for pair in negativePairs:

	normalizedValues = [pair[0]]
	for featureInd in range(1, len(pair)):
		if currentMin[featureInd] == 0 and currentMax[featureInd] == 0: #if the min/max are 0 for this feature, the normalized value should also be 0.
			normalizedValues.append(0)
			continue

		if currentMax[featureInd] - currentMin[featureInd] == 0:
			normalizedValues.append(0) #min and max should be the same, so no norm necessary
			continue

		normFeature = (float(pair[featureInd])-currentMin[featureInd])/(currentMax[featureInd]-currentMin[featureInd])
		normalizedValues.append(normFeature)

	normalizedNegativePairs.append(normalizedValues)

normalizedPositivePairs = np.array(normalizedPositivePairs, dtype='object')
normalizedNegativePairs = np.array(normalizedNegativePairs, dtype='object')

#try simple classification on this.
svTypes  = ['DEL', 'DUP', 'INV', 'ITX']
for svType in svTypes:

	if svType == 'DEL':
		classifier = RandomForestClassifier(random_state=785, n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
	elif svType == 'DUP':
		classifier = RandomForestClassifier(random_state=785, n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
	elif svType == 'INV':
		classifier = RandomForestClassifier(random_state=785, n_estimators= 200, min_samples_split=5, min_samples_leaf=4, max_features='auto', max_depth=10, bootstrap=True)
	elif svType == 'ITX':
		classifier = RandomForestClassifier(random_state=785, n_estimators= 1000, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
	else:
		print('Please provide a valid SV type')
		exit(1)


	#divide into chromosomes.

	positivePairsPerChromosome = dict()
	for pair in normalizedPositivePairs:

		splitPair = pair[0].split('_')

		if splitPair[12] != svType:
			continue

		if splitPair[1] not in positivePairsPerChromosome:
			positivePairsPerChromosome[splitPair[1]] = []

		allFeatures = []
		for featureInd in range(1, len(pair)):

			if featureInd in allowedFeatures:
				allFeatures.append(pair[featureInd])

		positivePairsPerChromosome[splitPair[1]].append(allFeatures)


	negativePairsPerChromosome = dict()
	matchPairsNeg = 0
	for pair in normalizedNegativePairs:

		splitPair = pair[0].split('_')

		if splitPair[12] != svType:
			continue

		if splitPair[1] not in negativePairsPerChromosome:
			negativePairsPerChromosome[splitPair[1]] = []

		allFeatures = []
		for featureInd in range(1, len(pair)):

			#if featureInd-1 in badIdx: #skip 0 variance
			#	continue

			if featureInd in allowedFeatures:
				allFeatures.append(pair[featureInd])

		negativePairsPerChromosome[splitPair[1]].append(allFeatures)

	chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
				   'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
				   'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
				   'chr20', 'chr21', 'chr22']



	#go through the chromosomes and make the positive/negative sets
	#subsample the negative set to the same size.
	aucs = []
	#for the test set, keep statistics for TPR and FPR to compare between methods
	totalTP = 0
	totalFP = 0
	totalTN = 0
	totalFN = 0
	totalP = 0
	totalN = 0
	for chromosome in chromosomes:

		if chromosome not in positivePairsPerChromosome or chromosome not in negativePairsPerChromosome:
			continue

		testPositivePairs = np.array(positivePairsPerChromosome[chromosome])
		testNegativePairs = np.array(negativePairsPerChromosome[chromosome])

		#randomly subsample
		randInd = random.sample(range(0, testNegativePairs.shape[0]), testPositivePairs.shape[0])
		testNegativeSubset = testNegativePairs[randInd]

		testSet = np.concatenate((testPositivePairs, testNegativeSubset))
		testLabels = [1]*testPositivePairs.shape[0] + [0]*testNegativeSubset.shape[0]

		totalP += testPositivePairs.shape[0]
		totalN += testNegativeSubset.shape[0]

		trainingSet = []
		trainingLabels = []
		for chromosome2 in positivePairsPerChromosome:

			if chromosome2 == chromosome:
				continue

			if chromosome2 not in positivePairsPerChromosome or chromosome2 not in negativePairsPerChromosome:
				continue

			chrPositivePairs = np.array(positivePairsPerChromosome[chromosome2])
			chrNegativePairs = np.array(negativePairsPerChromosome[chromosome2])

			#randomly subsample
			randInd = random.sample(range(0, chrNegativePairs.shape[0]), chrPositivePairs.shape[0])
			chrNegativeSubset = chrNegativePairs[randInd]

			#chrTrainingSet = np.concatenate((chrPositivePairs, chrNegativeSubset))

			#trainingSet.append(list(chrTrainingSet))

			for pair in chrPositivePairs:
				trainingSet.append(pair)
			for pair in chrNegativeSubset:
				trainingSet.append(pair)

			chrLabels = [1]*chrPositivePairs.shape[0] + [0]*chrNegativeSubset.shape[0]
			trainingLabels += chrLabels

		trainingSet = np.array(trainingSet)

		allPairs = np.concatenate((trainingSet, testSet))
		from sklearn.feature_selection import VarianceThreshold
		t = 0
		vt = VarianceThreshold(threshold=t)
		vt.fit(allPairs)
		idx = np.where(vt.variances_ > t)[0]
		badIdx = np.where(vt.variances_ <= t)[0]

		trainingSet = trainingSet[:,idx]
		testSet = testSet[:,idx]

		classifier.fit(trainingSet, trainingLabels)

		#check true/false positives/negatives
		predictions = classifier.predict(testSet)
		for labelInd in range(0, len(testLabels)):

			if testLabels[labelInd] == 1 and predictions[labelInd] == 1:
				totalTP += 1
			elif testLabels[labelInd] == 0 and predictions[labelInd] == 1:
				totalFP += 1
			elif testLabels[labelInd] == 1 and predictions[labelInd] == 0:
				totalFN += 1
			else:
				totalTN += 1

		fig, ax = plt.subplots()
		viz = plot_roc_curve(classifier, testSet, testLabels,
							 name='roc',
							 alpha=0.3, lw=1, ax=ax)

		#print('auc: ', np.mean(viz.roc_auc))
		aucs.append(np.mean(viz.roc_auc))
		plt.close()

	print(np.mean(aucs))

	#report on true positive and false positie rate
	print(totalTP, totalFP, totalFN, totalTN)
	tpr = totalTP / (totalTP + totalFN)
	fpr = totalFP / (totalTN + totalFP)
	print(tpr, fpr)

