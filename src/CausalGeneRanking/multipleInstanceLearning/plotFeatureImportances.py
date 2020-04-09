"""

	Make the polar plots showing the feature importances.
	This script can output the plots for the COSMIC, non-COSMIC, and all instance case.

	Based on a parameter, we can select to only output the plot from pre-made plotting data,
	or to also generate the plotting data.

"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
import random
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

path = sys.argv[5]
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVsGenes/')

import settings
from inputParser import InputParser

import matplotlib
matplotlib.use('Agg')

outDir = sys.argv[1]
generatePlottingData = sys.argv[2]
gainLoss = sys.argv[3] #ALL, GAIN, LOSS
cosmic = sys.argv[4] #ALL, COSMIC, NONCOSMIC
outputTop1000FeatureImportance = sys.argv[6]
degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object')

#set the right output files based on these parameters
if gainLoss == 'ALL':
	gainLossName = 'all'
elif gainLoss == 'GAIN':
	gainLossName = 'gain'
elif gainLoss == 'LOSS':
	gainLossName = 'loss'
else:
	print('invalid option for gains/losses')
	exit(1)

if cosmic == 'ALL':
	cosmicName = 'all'
elif cosmic == 'COSMIC':
	cosmicName = 'cosmic'
elif cosmic == 'NONCOSMIC':
	cosmicName = 'noncosmic'
else:
	print('invalid option for COSMIC')
	exit(1)

#maybe output the .npy files to a tmp dir to make them easier to find
outFilePrefix = outDir + '/' + gainLossName + '_' + cosmicName + '_'

if gainLossName != 'loss':
	svTypes = ['DEL', 'DUP', 'INV', 'ITX']
	
else:
	svTypes = ['INV', 'ITX']
svType = np.array(svTypes)
usedSVTypes = [] #use this to later determine which colors need to be used in the plot in case we skip an sv type, e.g.
#for cosmic when translocations are not linked to any cosmic gene.

if cosmicName != 'all':
	#read the cosmic files to split instances into cosmic/non-cosmic.
	cosmicGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
	cosmicGeneNames = []
	for gene in cosmicGenes:
		cosmicGeneNames.append(gene[3].name)


if generatePlottingData == "True":
	adjustedPValues = dict()
	allFeatureZScores = dict()

	svTypeInd = 0
	for svType in svTypes:
		#define the classifiers to use (from optimization)
		#would be nicer if these are in 1 file somewhere, since they are also used in another script

		if svType == 'DEL':
			clf = RandomForestClassifier(random_state=785, n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
			title = 'deletions'
		elif svType == 'DUP':
			#classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=2, min_samples_leaf=2, max_features='sqrt', max_depth=110, bootstrap=False)
			clf = RandomForestClassifier(random_state=785, n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
			#classifier = RandomForestClassifier(n_estimators= 1200, min_samples_split=2, min_samples_leaf=4, max_features='sqrt', max_depth=70, bootstrap=False)
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
		bagLabels = np.load(dataPath + '/bagLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		instances = np.load(dataPath + '/instances_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagPairLabels = np.load(dataPath + '/bagPairLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagMap = np.load(dataPath + '/bagMap_' + svType + '.npy', encoding='latin1', allow_pickle=True).item()
		filteredFeatures = np.loadtxt(dataPath + '/lowVarianceIdx_' + svType + '.txt')

		print(filteredFeatures)

		#train the classifier on the full dataset
		clf.fit(similarityMatrix, bagLabels)
		print('Classifier performance on all data: ', clf.score(similarityMatrix, bagLabels))
		#get the feature importances
		importances = clf.feature_importances_

		#rank these importances by score
		indices = np.argsort(importances)[::-1]

		###Figure S4
		if outputTop1000FeatureImportance == 'True':
			plt.figure()
			plt.title("Feature importances")
			plt.bar(range(similarityMatrix.shape[1])[0:1000], importances[indices[0:1000]],
				   color="r", align="center")
			plt.xticks(range(similarityMatrix.shape[1])[0:1000][::25], range(0,1000, 25), rotation=90)
			plt.tight_layout()
			plt.savefig(outDir + '/feature_importances_top1000_' + svType + '.svg')
	
			continue

		## because the low variance features were removed, we need to add them back here to show all features.
		#so these can be set to none, or masked for now.

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

		fixedInstances = np.array(fixedInstances)


		#### here we determine if we look at gains/losses, cosmic/non-cosmic
		# all these need to be parameters to automatically output the right plots
		top = 100
		topInstances = fixedInstances[indices[0:top]]

		filteredInstances = []
		uniqueGenes = []
		uniqueCosmicGenes = []
		cosmicGenes = []
		topPairLabels = []
		hallmarkCount = 0
		for instanceInd in range(0, topInstances.shape[0]):
			#get the ranked index of this instance
			rankedInstanceInd = indices[instanceInd]

			#get the label of the sv-gene pair this instance comes from
			bagLabel = bagPairLabels[bagMap[rankedInstanceInd]]
			splitPair = bagLabel.split('_')

			topPairLabels.append(bagLabel)

			shortPair = splitPair[7] + '_' + splitPair[0]

			if splitPair[0] not in uniqueGenes:
				uniqueGenes.append(splitPair[0])

			if gainLossName == 'gain':
				#in case of deletions and duplications, there are always noly gains.
				#so we only need to chekc for INV and ITX
				if svType not in ['DEL', 'DUP']:
					if topInstances[instanceInd][1] == 0:
						continue
			if gainLossName == 'loss':
				if topInstances[instanceInd][0] == 0:
					continue

			if cosmicName == 'cosmic':

				#check in the pair labels if this is a cosmic gene.
				if splitPair[0] in cosmicGeneNames:
					print(bagLabel)
					filteredInstances.append(topInstances[instanceInd])

			elif cosmicName == 'noncosmic':
				if splitPair[0] not in cosmicGeneNames:
					filteredInstances.append(topInstances[instanceInd])
			else:
				filteredInstances.append(topInstances[instanceInd])

		np.savetxt(outDir + '/pairLabels_top100_' + svType + '.txt', topPairLabels, fmt='%s', delimiter='\t')
		filteredInstances = np.array(filteredInstances)
		print('number of instances used: ', filteredInstances.shape)

		if len(filteredInstances) < 1:
			continue
		usedSVTypes.append(svType)

		#remove the gains and losses features, we do not need those anymore.
		filteredInstances = np.delete(filteredInstances, 1, 1)
		filteredInstances = np.delete(filteredInstances, 0, 1)


		#compute the percentages in these top X instances
		avgInstances = np.sum(filteredInstances, axis=0)

		totalInstances = avgInstances / filteredInstances.shape[0]
		print(totalInstances)

		#100 times, randomly sample
		#per feature, have a distribution
		nullDistributions = dict()
		for i in range(0,100):

			randomHallmarkCount = 0

			if i == 0:
				for featureInd in range(0, len(totalInstances)):
					nullDistributions[featureInd] = []

			#sample as much random instances as in our filtered instances
			randomIndices = random.sample(range(0,fixedInstances.shape[0]), filteredInstances.shape[0])

			randomTopInstances = fixedInstances[randomIndices]

			#here also skip gains/losses
			randomTopInstances = np.delete(randomTopInstances, 1, 1)
			randomTopInstances = np.delete(randomTopInstances, 0, 1)

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

		allFeatureZScores[svType] = featureZScores
		adjustedPValues[svType] = pAdjusted
		
		svTypeInd += 1

	#set the right names based on the parameters
	np.save(outFilePrefix + 'featureZScores.npy', allFeatureZScores)
	np.save(outFilePrefix + 'adjustedPValues.npy', adjustedPValues)

#load the data and make the polar plots
#use the parameters to determine which plot to make

### ALWAYS load the right ALL set here of gains/losses to normalize properly
allSet = outDir + '/all_' + cosmicName + '_'

allFeatureZScores = np.load(allSet + 'featureZScores.npy', allow_pickle=True, encoding='latin1').item()
adjustedPValues = np.load(allSet + 'adjustedPValues.npy', allow_pickle=True, encoding='latin1').item()

#if we are looking at cosmic/non-cosmic, do not show the cosmic feature.
#these labels are never added to the correct place in the plot for some reason, so I
#post-edit them and give them clearer names there. Too long names will make the plot unreadable.

xlabels = ['cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
			'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
			'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'enhancerType', 'eQTLType', 'superEnhancerType',
			'instCount']

#to plot, show the p-values, direction based on the z-score.
overallMin = float('inf')
overallMax = 0
zeroOffset = 0
for svType in usedSVTypes:

	pAdjusted = adjustedPValues[svType]

	directionalAdjustedP = -np.log(pAdjusted) * np.sign(allFeatureZScores[svType])

	if np.min(directionalAdjustedP) < overallMin:
		overallMin = np.min(directionalAdjustedP)
	if np.max(directionalAdjustedP) > overallMax:
		overallMax = np.max(directionalAdjustedP)


allFeatureZScores = np.load(outFilePrefix + 'featureZScores.npy', allow_pickle=True, encoding='latin1').item()
adjustedPValues = np.load(outFilePrefix + 'adjustedPValues.npy', allow_pickle=True, encoding='latin1').item()

scaledP = dict()

for svType in usedSVTypes:

	pAdjusted = adjustedPValues[svType]

	directionalAdjustedP = -np.log(pAdjusted) * np.sign(allFeatureZScores[svType])

	directionalAdjustedP += np.abs(overallMin) + 50

	if len(np.where(pAdjusted == 1)) > 0:
		zeroOffsetInd = np.where(pAdjusted == 1)[0][0]

		zeroOffset = directionalAdjustedP[zeroOffsetInd]

	scaledP[svType] = directionalAdjustedP

border = zeroOffset

signBorderTop = -np.log(0.05) + zeroOffset
signBorderBottom = border - (signBorderTop - border)

blockSize = 360 / len(directionalAdjustedP)
xRange = np.arange(0, 360, blockSize)

#add the last value, which is missed in the range
xRange = np.append(xRange, xRange[xRange.shape[0]-1]+blockSize)
centers = np.deg2rad(np.ediff1d(xRange)//2 + xRange[:-1])

fig = plt.figure(figsize=(15,13))
ax = fig.add_subplot(111, projection='polar')

allColors = np.array(['blue', 'red', 'magenta', 'black'])
allOffset = np.array([-0.02, -0.01, 0.01, 0.02])

colors = []
offset = []
for svType in usedSVTypes:
	if svType == 'DEL':
		colors.append(allColors[0])
		offset.append(allOffset[0])
	if svType == 'DUP':
		colors.append(allColors[1])
		offset.append(allOffset[1])
	if svType == 'INV':
		colors.append(allColors[2])
		offset.append(allOffset[2])
	if svType == 'ITX':
		colors.append(allColors[3])
		offset.append(allOffset[3])

ind = 0
for svType in usedSVTypes:
	print('area')
	area = 2 * (1 + (-np.log(adjustedPValues[svType])))
	print(area)
	ax.scatter(centers+offset[ind], scaledP[svType], color=colors[ind], alpha=0.4, s=area)
	ind += 1

ax.set_xticks(centers)
ax.set_xticklabels(xlabels, fontsize=5)

ax.set_yticklabels(['p < 0.05'])
ax.set_yticks([signBorderBottom, border, signBorderTop])
gridlines = ax.yaxis.get_gridlines()
gridlines[0].set_color("red")
gridlines[0].set_linewidth(0.5)
gridlines[0].set_linestyle('--')

gridlines[2].set_color("red")
gridlines[2].set_linewidth(0.5)
gridlines[2].set_linestyle('--')

ax.set_rorigin(-200)
ax.set_theta_zero_location('N', offset=10)
plt.savefig(outFilePrefix + 'featureImportances.svg')