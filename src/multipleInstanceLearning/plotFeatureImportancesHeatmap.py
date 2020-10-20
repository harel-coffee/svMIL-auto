## plot a heatmap with the top 100 instance features across all SV types.

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

def generateHeatmap(cancerTypes, svTypes = ['DEL', 'DUP', 'INV', 'ITX']):

	#get the significances for each cancer type
	pValuesPerCancerType = dict()
	for cancerType in cancerTypes:
		print('cancer type: ', cancerType)
		#first get the instances and their ranking across all SV types
		importances, instances = getFeatureImportances(svTypes, cancerType)
		#then compute the significances of the top 100 to random 100 instances

		pValues, zScores = computeFeatureSignificances(importances, instances, 100)
		pValuesPerCancerType[cancerType] = [pValues, zScores]

	#then make a heatmap plot of the significances.
	plotHeatmap(pValuesPerCancerType)

def plotHeatmap(pValuesPerCancerType):

	#re-format the p-values to a binary style for plotting
	significanceMatrix = []
	for cancerType in pValuesPerCancerType:
		significances = []
		pValues = pValuesPerCancerType[cancerType][0]
		zScores = pValuesPerCancerType[cancerType][1]

		#below this we call it 'very' significant.
		signCutoff = 1e-5
		for pValueInd in range(0, len(pValues)):
			pValue = pValues[pValueInd]
			if pValue < 0.05 and zScores[pValueInd] > 0:

				if pValue < signCutoff:
					significances.append(2)
				else:
					significances.append(1)
			elif pValue < 0.05 and zScores[pValueInd] < 0:

				if pValue < signCutoff:
					significances.append(-2)
				else:
					significances.append(-1)

			else:
				significances.append(0)

		significanceMatrix.append(significances)

	significanceMatrix = np.array(significanceMatrix)

	fig =plt.figure(figsize=(15,10))
	
	data = pd.DataFrame(significanceMatrix) #exclude translocations, these are not there for germline.
	g=sns.heatmap(data,annot=False,square=True, linewidths=0.5,
				  cmap=ListedColormap(['#0055d4ff', '#0055d47d', '#f7f6f6ff', '#c8373780', '#c83737ff']),
				  yticklabels=list(pValuesPerCancerType.keys()))

	g.set_yticklabels(g.get_yticklabels(), horizontalalignment='right',fontsize='small')
	plt.xticks(np.arange(0, len(pValues))+0.5, ['Gains', 'Losses', 'cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
			'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
			'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k26me3_s', 'enhancerType', 'eQTLType', 'superEnhancerType',
			'instCount'], rotation=45, horizontalalignment='right')

	plt.tight_layout()


	plt.savefig('featureImportanceHeatmap.svg')
	plt.show()

	
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


	print(pAdjusted)
	return pAdjusted, featureZScores

def getFeatureImportances(svTypes, cancerType):

	#set the directory to look in for this cancer type
	outDir = sys.argv[1] + '/' + cancerType
	
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

cancerTypes = ['HMF_Breast', 'HMF_Ovary', 'HMF_Liver', 'HMF_Lung', 'HMF_Colorectal',
				   'HMF_UrinaryTract', 'HMF_Prostate', 'HMF_Esophagus', 'HMF_Skin',
				   'HMF_Pancreas', 'HMF_Uterus', 'HMF_Kidney', 'HMF_NervousSystem']

generateHeatmap(cancerTypes)

		