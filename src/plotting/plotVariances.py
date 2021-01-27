"""
	Script to make a plot of the feature variances across the bags.
	Use the normalized bags as input, and then show which features have low
	variances, and are therefore not useful for the model.
"""

import sys
import os
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.feature_selection import VarianceThreshold

import matplotlib
matplotlib.use('Agg')

outDir = sys.argv[1]

#1. First gather the normalized bags

svTypes = ['DEL', 'DUP', 'INV', 'ITX']

#input the normalized bags
with open(outDir + '/linkedSVGenePairs/normalizedBags.pkl', 'rb') as handle:
	bagDict = pkl.load(handle)

#get the information for the bag labels
degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object') #labels

mutDir = outDir + '/patientGeneMutationPairs/'
cnvPatientsAmp = np.load(mutDir + 'cnvPatientsAmp.npy', allow_pickle=True, encoding='latin1').item()
svPatientsDup = np.load(mutDir + 'svPatientsDup.npy', allow_pickle=True, encoding='latin1').item()
svGenePairs = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_', dtype='object')
splitSVGenePairs = []
for pair in svGenePairs:

	splitPair = pair[0].split('_')

	splitSVGenePairs.append(splitPair[0] + '_' + splitPair[7] + '_' + splitPair[12])

#2. Divide into positive and negative instances
for svType in svTypes:

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

		if svType != '' and svType != 'ALL':
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


					instances.append(instance)

				if len(instances) < 1:
					continue

				###Here do an extra check:
				#we only look at TADs with SVs across the boundary when computing z-scores, so those z-scores are in the set.
				#BUT some of these genes are not actually affected by the SV, since this doesn't lead to
				#regulatory elements gained/lost. SO, we need to remove those here to get the actual pairs.
				#This only goes wrong for duplications, because we also keep CNV amps that could be the same event,
				#but then the duplication does not lead to gains/losses, while the CNV amp does because it is slightly
				#longer. So if there is evidence of a cnv AMP, but no non-coding duplication linked, we can remove
				#this as a positive pair.
				if splitPair[7] not in cnvPatientsAmp:
					positiveBagPairNames.append(pair)
					positiveBags.append(instances)
				else:

					dupMatch = splitPair[0] + '_' + splitPair[7] + '_DUP'

					if splitPair[0] in cnvPatientsAmp[splitPair[7]] and dupMatch not in splitSVGenePairs:
						negativeBags.append(instances)
						negativeBagPairNames.append(pair)

					else:
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

	#add the number of instances per bag as feature to the instances, which
	#is not part of the normalized bags. 
	#normalize by the number of bags equal to how other features are normalized. 
	for bag in positiveBags:
		instCount = len(bag)

		for instance in bag:
			instance.append(instCount / positiveBags.shape[0])

	for bag in negativeBags:
		instCount = len(bag)

		for instance in bag:
			instance.append(instCount / negativeBags.shape[0])

	#remove instances with no variance
	posInstances = np.vstack(positiveBags)
	negInstances = np.vstack(negativeBags)

	allInstances = np.concatenate((posInstances, negInstances))

	#3. Calculate variances
	t = 0
	vt = VarianceThreshold(threshold=t)
	vt.fit(allInstances)

	#normalize variances for visualization purposes
	normalizedVariances = np.log(vt.variances_+0.0000001)

	#Then make a plot of the variances for each SV type model. 
	plotData = []
	for ind in range(0, len(vt.variances_)):
		plotData.append([ind, normalizedVariances[ind]])

	data = pd.DataFrame(plotData)
	data.columns = ['Features', 'Variance']

	sns.scatterplot(data=data, x='Features',
					y='Variance', s = 30, edgecolor = 'k', facecolor='k')
	plt.xticks(np.arange(0, len(normalizedVariances)), ['Gains', 'Losses', 'CpG', 'TF', 'HiC', 'CTCF', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
				'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'RNA pol II',
				'Enhancer strength', 'CTCF strength', 'RNA pol II strength', 'h3k9me3 strength', 'h3k4me3 strength', 'h3k27ac strength', 'h3k27me3 strength', 'h3k4me1 strength', 'h3k36me3 strength', 'Enhancer type', 'eQTL type', 'Super enhancer type',
				'Instance count'], rotation=45, horizontalalignment='right')
	plt.ylabel('Log(variance)')
	
	plt.tight_layout()
	plt.savefig('figureS1_' + svType + '.svg')
	plt.clf()




