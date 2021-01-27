"""
	Gather all bags per SV type. Group into positive and negative. Then take the averages per bag across the instances. Plot in a boxplot

"""

import sys
import numpy as np
import pickle as pkl
import random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

outDir = sys.argv[1]

#1. load the normalized bags
with open(outDir + '/linkedSVGenePairs/normalizedBags.pkl', 'rb') as handle:
	bagDict = pkl.load(handle)

#get the information for the bag labels
degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object') #labels

print(degPairs)

print("initial number of bags: ", len(bagDict))
print('initial deg pairs: ', degPairs.shape[0])

mutDir = outDir + '/patientGeneMutationPairs/'
cnvPatientsAmp = np.load(mutDir + 'cnvPatientsAmp.npy', allow_pickle=True, encoding='latin1').item()
svPatientsDup = np.load(mutDir + 'svPatientsDup.npy', allow_pickle=True, encoding='latin1').item()
svGenePairs = np.loadtxt(outDir + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_', dtype='object')
splitSVGenePairs = []
for pair in svGenePairs:

	splitPair = pair[0].split('_')

	splitSVGenePairs.append(splitPair[0] + '_' + splitPair[7] + '_' + splitPair[12])

finalOutDir = outDir + '/multipleInstanceLearning/similarityMatrices/'
svTypes = ['DEL', 'DUP', 'INV', 'ITX']
#Generate the similarity matrices for the SV types
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
			#if float(degPairInfo[5]) > 2 or float(degPairInfo[5]) < -2:
				#go through the instances of this SV-gene pair, and include only those that have gains and losses, and more than 1 instance. This should in principle not happen, but good to keep a check.
				instances = []
				for instance in bagDict[pair]:

					if instance[0] == 0 and instance[1] == 0:
						continue


					instances.append(instance)

				if len(instances) < 1:
					continue

				###Here do an extra check:
				#to make fig2A, we only look at TADs with SVs across the boundary, so those z-scores are in the set.
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
	positiveBagPairNames = np.array(positiveBagPairNames)
	negativeBagPairNames = np.array(negativeBagPairNames)

	posInstances = np.vstack(positiveBags)
	negInstances = np.vstack(negativeBags)
	allInstances = np.concatenate((posInstances, negInstances))


	#fail-safe in case there are not enough SVs of this type
	if positiveBags.shape[0] < 2 or negativeBags.shape[0] < 2:
		continue
	
	print(positiveBags.shape, negativeBags.shape)
	
	random.seed(785) #somehow the global setting doesn't work in the second loop? so set it here.

	#subsample the negative set to the same number of positives.
	negativeBagsSubsampled = np.random.choice(negativeBags, positiveBags.shape[0])

	negativeBagsSubsampleInd = np.random.choice(np.arange(negativeBags.shape[0]), positiveBags.shape[0])
	negativeBagsSubsampled = negativeBags[negativeBagsSubsampleInd]
	
	
	negativeBagPairNamesSubsampled = negativeBagPairNames[negativeBagsSubsampleInd]
	
	#add the number of instances per bag as feature to the instances
	for bag in positiveBags:
		instCount = len(bag)
	
		for instance in bag:
			
			instance.append(instCount / positiveBags.shape[0])
			#print(instance)
			#exit()
	
	
	for bag in negativeBagsSubsampled:
		instCount = len(bag)
	
		for instance in bag:
			instance.append(instCount / negativeBagsSubsampled.shape[0])
	
	bagPairLabelsSubsampled = np.concatenate((positiveBagPairNames, negativeBagPairNamesSubsampled))
	
	#merge the bags so that we can easily get to 1 similarity matrix and do all-to-all computations
	bagsSubsampled = np.concatenate((positiveBags, negativeBagsSubsampled))
	
	#instances = np.vstack(bagsSubsampled)
	#print(instances.shape)
	features = ['Losses', 'Gains', 'cpg', 'tf', 'ctcf', 'dnaseI', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1',
				'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
				'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'enhancerType', 'eQTLType', 'superEnhancerType',
				'instCount']
	print(len(features))

	plotData = []
	labels = []
	for bag in positiveBags:
		avg = np.mean(bag, axis = 0)
		labels.append(1)
		plotData.append(avg)
	for bag in negativeBagsSubsampled:
		avg = np.mean(bag, axis = 0)
		labels.append(0)
		plotData.append(avg)

	#plotData = np.array(plotData)
	#print(plotData[0])
	
	data = pd.DataFrame(plotData)
	print(data)
	data.columns = features
	data['Label'] = labels

	df_long = pd.melt(data, "Label", var_name="a", value_name="c")
	print(df_long)
	sns.boxplot(x="a", hue="Label", y="c", data=df_long)
	plt.xticks(np.arange(0, len(features)), features, rotation=90)
	plt.show()

	
	
