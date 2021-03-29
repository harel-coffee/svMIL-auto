import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import random
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import os.path
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap

import matplotlib
matplotlib.use('Agg')

def computeFeatureSignificances(importances, instances, top):
	"""
		Compute the significance of the total occurrence of features in the instances
		within the provided top X instances with highest feature importance compared to
		X random instances.

		Parameters:
		- importances: allImportances output from self.getFeatureImportances()
		- instances: allInstances output from self.getFeatureImportances()
		- top: integer value of top X instances with highest importance to select.

		Return:
		- pAdjusted: bonferroni corrected p-values of each feature in the true instances
		compared to the random instances.
		- featureZScores: z-scores used to compute the p-values.
	"""

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
	random.seed(785)
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

	from math import sqrt


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
		if pValue == 0:
			pValue = 0.00000000000000001
		featurePValues.append(pValue)

	#do MTC on the p-values

	reject, pAdjusted, _, _ = multipletests(featurePValues, method='bonferroni')

	return pAdjusted, featureZScores


def getFeatureImportances(svTypes, cancerType, gainLossType):
	"""
		For the given cancer type, compute the random forest feature importances for
		each model of each SV type. Obtain the full similarity matrix (e.g. not subsampled)
		and train the RF classifier. Merge the importances of each SV type
		with the rest, and disregard SV type information as the importances point to
		similar instances for SV types.

		Parameters:
		- svTypes: list with SV types to get the importances for
		- cancerType: name of cancer type output folder to get data from

		Return:
		- allInstances: numpy array with all instances across SV types concatenated
		- allImportances: feature importance score of instances in allInstances

	"""

	#set the directory to look in for this cancer type
	outDir = 'output/' + cancerType

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


		#check for sv types for which we have no SVs
		if not os.path.isfile(dataPath + '/similarityMatrix_' + svType + '.npy'):
			continue

		similarityMatrix = np.load(dataPath + '/similarityMatrix_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagLabels = np.load(dataPath + '/bagLabelsSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		instances = np.load(dataPath + '/instancesSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagPairLabels = np.load(dataPath + '/bagPairLabelsSubsampled_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagMap = np.load(dataPath + '/bagMap_' + svType + '.npy', encoding='latin1', allow_pickle=True).item()
		filteredFeatures = np.loadtxt(dataPath + '/lowVarianceIdx_' + svType + '.txt')

		#train the classifier on the full dataset
		clf.fit(similarityMatrix, bagLabels)
		#get the feature importances
		importances = list(clf.feature_importances_)

		allImportances += importances

		#because low index features are removed, add them back here if necessary
		#to retain an overview of all used features
		fixedInstances = []
		for instance in instances:

			finalLength = len(instance) + filteredFeatures.size
			instanceMal = np.zeros(finalLength) #mal including missing features
			addedFeatures = 0
			for featureInd in range(0, finalLength):

				if featureInd in filteredFeatures:
					#if the filtered feature is number 1 and this is a DEL or DUP, it is always a
					#gain because these sv types cannot cause losses. Therefore the value is 1.
					if svType in ['DEL', 'DUP'] and featureInd == 1:
						instanceMal[featureInd] = 1
					else:
						instanceMal[featureInd] = 0
					addedFeatures += 1
				else:
					instanceMal[featureInd] = instance[featureInd-addedFeatures]

			fixedInstances.append(instanceMal)

		allInstances += fixedInstances

	allImportances = np.array(allImportances)
	allInstances = np.array(allInstances)


	###For separating based on gains/losses
	filteredInstances = []
	filteredImportances = []
	for instanceInd in range(0, allInstances.shape[0]):
		instance = allInstances[instanceInd]

		if gainLossType == 'loss':
			if instance[0] != 1:
				continue
		elif gainLossType == 'gain':
			if instance[1] != 1:
				continue
		filteredInstances.append(instance)
		filteredImportances.append(allImportances[instanceInd])

	allImportances = np.array(filteredImportances)
	allInstances = np.array(filteredInstances)

	return allImportances, allInstances

def computeOverallMin(pValuesPerCancerType):

	overallMin = float('inf')
	overallMax = 0

	for cancerType in pValuesPerCancerType:

		pValues = pValuesPerCancerType[cancerType][0]
		zScores = pValuesPerCancerType[cancerType][1]
		directionalAdjustedP = -np.log(pValues) * np.sign(zScores)

		if np.min(directionalAdjustedP) < overallMin:
			overallMin = np.min(directionalAdjustedP)
		if np.max(directionalAdjustedP) > overallMax:
			overallMax = np.max(directionalAdjustedP)

	print(overallMin, overallMax)

	return overallMin, overallMax

#compute the zero-crossing based on all cancer types.
def determineZeroCrossing(pValuesPerCancerType, overallMin, overallMax):

	zeroOffset = 0
	for cancerType in pValuesPerCancerType:

		pValues = pValuesPerCancerType[cancerType][0]
		zScores = pValuesPerCancerType[cancerType][1]

		directionalAdjustedP = -np.log(pValues) * np.sign(zScores)
		#normalize between 0 and 1
		directionalAdjustedP = (directionalAdjustedP - overallMin) / (overallMax - overallMin)

		if len(np.where(pValues == 1)) > 0:
			zeroOffsetInd = np.where(pValues == 1)[0][0]

			zeroOffset = directionalAdjustedP[zeroOffsetInd]
			break #we only need this for one cancer type, it should be the same for all.

	return zeroOffset

def plotFeatureWheel(pValuesPerCancerType, overallMin, overallMax, zeroOffset, svType, gainLossType):

	#Scale the p-values to match across cancer types
	scaledP = dict()
	for cancerType in pValuesPerCancerType:

		pValues = pValuesPerCancerType[cancerType][0]
		zScores = pValuesPerCancerType[cancerType][1]

		directionalAdjustedP = -np.log(pValues) * np.sign(zScores)
		print('p-values: ', pValues)
		print('before: ', directionalAdjustedP)
		#normalize between 0 and 1
		directionalAdjustedP = (directionalAdjustedP - overallMin) / (overallMax - overallMin)
		print('after:', directionalAdjustedP)

		#remove gains and losses
		directionalAdjustedP = np.delete(directionalAdjustedP, 1)
		directionalAdjustedP = np.delete(directionalAdjustedP, 0)

		scaledP[cancerType] = directionalAdjustedP

	border = zeroOffset

	#for the borders, also normalize between 0 and 1
	signBorderTop = (((-np.log(0.05) - overallMin) / (overallMax-overallMin)))
	signBorderBottom = border - (signBorderTop - border)

	print(signBorderTop)
	print(signBorderBottom)
	print(border)


	blockSize = 360 / len(directionalAdjustedP)
	xRange = np.arange(0, 360, blockSize)

	#add the last value, which is missed in the range
	xRange = np.append(xRange, xRange[xRange.shape[0]-1]+blockSize)
	centers = np.deg2rad(np.ediff1d(xRange)//2 + xRange[:-1])

	fig = plt.figure(figsize=(15,13))
	ax = fig.add_subplot(111, projection='polar')

	# allColors = np.array(['blue', 'red', 'magenta', 'black'])
	# allOffset = np.array([-0.02, -0.01, 0.01, 0.02])
	#
	# colors = []
	# offset = []
	# for svType in usedSVTypes:
	# 	if svType == 'DEL':
	# 		colors.append(allColors[0])
	# 		offset.append(allOffset[0])
	# 	if svType == 'DUP':
	# 		colors.append(allColors[1])
	# 		offset.append(allOffset[1])
	# 	if svType == 'INV':
	# 		colors.append(allColors[2])
	# 		offset.append(allOffset[2])
	# 	if svType == 'ITX':
	# 		colors.append(allColors[3])
	# 		offset.append(allOffset[3])

	palette=sns.color_palette("hls", len(pValuesPerCancerType))
	colors = palette.as_hex()
	filledMarkers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')

	#offset = [0]*len(pValuesPerCancerType)
	offset = [-0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06]

	import math
	ind = 0
	#for svType in usedSVTypes:
	#	print('area')
		# if cosmicName == 'cosmic':
		# 	area = 10 * (1 + (-np.log(adjustedPValues[svType])))
		# else:
	for cancerType in pValuesPerCancerType:
		pValues = pValuesPerCancerType[cancerType][0]
		#remove gains and losses
		pValues = np.delete(pValues, 1)
		pValues = np.delete(pValues, 0)
		area = 2 * (1 + (-np.log(pValues)))

		ax.scatter(centers+offset[ind], scaledP[cancerType], color=colors[ind], alpha=0.8, s=area, marker=filledMarkers[ind])
		ind += 1

	xlabels = ['CpG', 'TF', 'CTCF', 'DNAseI', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1',
				'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repressed', 'Transcribed', 'RNA pol II',
				'CTCF strength', 'RNA pol II strength', 'h3k4me3 strength', 'h3k27ac strength', 'h3k27me3 strength', 'h3k4me1 strength', 'Enhancer type', 'eQTL type', 'Super enhancer type',
				'Instance count']

	ax.set_xticks(centers)
	ax.set_xticklabels(xlabels, fontsize=5)

	ax.set_ylim(-0.5, 1.2)
	ax.set_yticklabels(['p < 0.05'])
	ax.set_yticks([signBorderBottom, border, signBorderTop])
	gridlines = ax.yaxis.get_gridlines()
	gridlines[0].set_color("red")
	gridlines[0].set_linewidth(0.5)
	gridlines[0].set_linestyle('--')

	gridlines[2].set_color("red")
	gridlines[2].set_linewidth(0.5)
	gridlines[2].set_linestyle('--')

	#ax.set_rorigin(-200)
	#ax.set_theta_zero_location('N', offset=10)
	ax.set_theta_zero_location('N')
	plt.savefig('output/figures/' + svType + '_' + gainLossType + '.svg')


def plotFigure(svTypes, cancerTypes, gainLossType):

	pValuesPerCancerType = dict()
	for cancerType in cancerTypes:
		print('Processing cancer type: ', cancerType)
		#first get the instances and their ranking across all SV types
		importances, instances = getFeatureImportances(svTypes, cancerType, gainLossType)
		#then compute the significances of the top 100 to random 100 instances

		pValues, zScores = computeFeatureSignificances(importances, instances, 100)
		pValuesPerCancerType[cancerType] = [pValues, zScores]

	overallMin, overallMax = computeOverallMin(pValuesPerCancerType)
	zeroOffset = determineZeroCrossing(pValuesPerCancerType, overallMin, overallMax)

	plotFeatureWheel(pValuesPerCancerType, overallMin, overallMax, zeroOffset, "allSVTypes", gainLossType)


#1. Make the figures for each combination
svTypes = ['DEL', 'DUP', 'INV', 'ITX']
cancerTypes = ['HMF_Breast_hmec', 'HMF_Ovary_ov', 'HMF_Lung_luad', 'HMF_Colorectal_coad',
 				'HMF_UrinaryTract_urinaryTract', 'HMF_Prostate_prostate',
 				'HMF_Esophagus_esophagus', 'HMF_Skin_skin',
 				'HMF_Pancreas_pancreas', 'HMF_Uterus_uterus',
 				'HMF_Kidney_kidney', 'HMF_NervousSystem_nervousSystem']

plotFigure(svTypes, cancerTypes, 'gain')
plotFigure(svTypes, cancerTypes, 'loss')

#repeat for CTCF

cancerTypes = ['HMF_Breast_CTCF']

plotFigure(svTypes, cancerTypes, 'gain')
plotFigure(svTypes, cancerTypes, 'loss')




