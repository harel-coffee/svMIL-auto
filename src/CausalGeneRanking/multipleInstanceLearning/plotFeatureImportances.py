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

outDir = sys.argv[1]
generatePlottingData = sys.argv[2]
gainLoss = sys.argv[3] #ALL, GAIN, LOSS
cosmic = sys.argv[4] #ALL, COSMIC, NONCOSMIC
degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object')

#set the right output files based on these parameters
if gainLoss == 'ALL':
	gainLossName = 'all'
elif gainLoss == 'GAIN':
	gainLossName = 'gain'
elif gainLossName == 'LOSS':
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

svTypes = ['DEL', 'DUP', 'INV', 'ITX']
svTypes = ['ITX']

if generatePlottingData == "True":
	adjustedPValues = dict()
	allFeatureZScores = dict()

	for svType in svTypes:
		#define the classifiers to use (from optimization)
		#would be nicer if these are in 1 file somewhere, since they are also used in another script

		if svType == 'DEL':
			clf = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
			title = 'deletions'
		elif svType == 'DUP':
			#classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=2, min_samples_leaf=2, max_features='sqrt', max_depth=110, bootstrap=False)
			clf = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
			#classifier = RandomForestClassifier(n_estimators= 1200, min_samples_split=2, min_samples_leaf=4, max_features='sqrt', max_depth=70, bootstrap=False)
			title = 'duplications'
		elif svType == 'INV':
			clf = RandomForestClassifier(n_estimators= 200, min_samples_split=5, min_samples_leaf=4, max_features='auto', max_depth=10, bootstrap=True)
			title = 'inversions'
		elif svType == 'ITX':
			clf = RandomForestClassifier(n_estimators= 1000, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
			title = 'translocations'
		else:
			clf = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
			title = 'All SV types'

		#load the similarity matrix of this SV type
		dataPath = outDir + '/multipleInstanceLearning/similarityMatrices/'
		similarityMatrix = np.load(dataPath + '/similarityMatrix_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagLabels = np.load(dataPath + '/bagLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		instances = np.load(dataPath + '/instances_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagPairLabels = np.load(dataPath + '/bagPairLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)
		bagMap = np.load(dataPath + '/bagMap_' + svType + '.npy', encoding='latin1', allow_pickle=True).item()

		#train the classifier on the full dataset
		clf.fit(similarityMatrix, bagLabels)
		print('Classifier performance on all data: ', clf.score(similarityMatrix, bagLabels))
		#get the feature importances
		importances = clf.feature_importances_

		#rank these importances by score
		indices = np.argsort(importances)[::-1]

		########### we eventually need this part for a supplementary figure (use parameter switch)

		# plt.figure()
		# plt.title("Feature importances")
		# plt.bar(range(similarityMatrix.shape[1])[0:1000], importances[indices[0:1000]],
		# 	   color="r", align="center")
		# plt.xticks(range(similarityMatrix.shape[1])[0:1000][::25], range(0,1000, 25), rotation=90)
		# plt.tight_layout()
		# plt.savefig('feature_importances_top1000_' + svType + '.svg')
		# plt.show()
		#
		# exit()


		#split the type back into 4 features
		# enhancerTypes = []
		# eQTLTypes = []
		# promoterTypes = []
		# superEnhancerTypes = []
		# uniqueGenes = []
		# uniqueCosmicGenes = []
		# instanceInd = 0
		# for instance in instances:
		#
		# 	if instance[33] == 0:
		# 		enhancerTypes.append(1)
		# 		eQTLTypes.append(0)
		# 		promoterTypes.append(0)
		# 		superEnhancerTypes.append(0)
		# 	elif instance[33] > 0 and instance[33] < 0.34:
		# 		enhancerTypes.append(0)
		# 		eQTLTypes.append(0)
		# 		promoterTypes.append(1)
		# 		superEnhancerTypes.append(0)
		# 	elif instance[33] > 0.33 and instance[33] < 0.68:
		# 		enhancerTypes.append(0)
		# 		eQTLTypes.append(1)
		# 		promoterTypes.append(0)
		# 		superEnhancerTypes.append(0)
		# 	else:
		# 		enhancerTypes.append(0)
		# 		eQTLTypes.append(0)
		# 		promoterTypes.append(0)
		# 		superEnhancerTypes.append(1)
		#
		# newInstances = np.zeros([instances.shape[0], instances.shape[1]+4])
		#
		# newInstances[:,0:instances.shape[1]] = instances
		#
		# newInstances[:,instances.shape[1]] = enhancerTypes
		# newInstances[:,instances.shape[1]+1] = promoterTypes
		# newInstances[:,instances.shape[1]+2] = eQTLTypes
		# newInstances[:,instances.shape[1]+3] = superEnhancerTypes
		#
		# #reomove the old type column
		# newInstances = np.delete(newInstances, 33, 1)

		#top X instances to look into
		instanceCount = 100

		#get the instances by their instance sorting
		topInstances = instances[indices[0:instanceCount]]

		#### here we determine if we look at gains/losses, cosmic/non-cosmic
		# all these need to be parameters to automatically output the right plots

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

			if topInstances[instanceInd][topInstances.shape[1]-1] > 0:
				cosmicGenes.append(splitPair[0])
				if splitPair[0] not in uniqueCosmicGenes:
					uniqueCosmicGenes.append(splitPair[0])

			#get z-score
			degPairInfo = degPairs[degPairs[:,0] == shortPair][0]

			#cosmic yes/no
			#if topInstances[instanceInd,33] > 0:
				#if topInstances[instanceInd,1] > 0:
			filteredInstances.append(topInstances[instanceInd])

		filteredInstances = np.array(filteredInstances)
		print('number of instances used: ', filteredInstances.shape)

		#remove the gains and losses features, we do not need those anymore.
		filteredInstances = np.delete(filteredInstances, 1, 1)
		filteredInstances = np.delete(filteredInstances, 0, 1)

		print(svType, ':')
		print('number of cosmic genes (unique): ', len(uniqueCosmicGenes))
		print('number of cosmic genes per instance: ', len(cosmicGenes))

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
			randomIndices = random.sample(range(0,instances.shape[0]), filteredInstances.shape[0])

			randomTopInstances = instances[randomIndices]

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
if cosmicName == 'all':
	xlabels = ['cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
			   'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
			   'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'enhancerType', 'promoterType', 'eQTLType', 'superEnhancerType', 'cosmic']
else:
	xlabels = ['cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
			   'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
			   'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'enhancerType', 'promoterType', 'eQTLType', 'superEnhancerType']

print(xlabels)
print(len(xlabels))
print(allFeatureZScores['ITX'])
print(len(allFeatureZScores['ITX']))
exit()


#to plot, show the p-values, direction based on the z-score.
overallMin = float('inf')
overallMax = 0
zeroOffset = 0
for svType in svTypes:

	pAdjusted = adjustedPValues[svType]

	directionalAdjustedP = -np.log(pAdjusted) * np.sign(allFeatureZScores[svType])

	if np.min(directionalAdjustedP) < overallMin:
		overallMin = np.min(directionalAdjustedP)
	if np.max(directionalAdjustedP) > overallMax:
		overallMax = np.max(directionalAdjustedP)


allFeatureZScores = np.load(outFilePrefix + 'featureZScores.npy', allow_pickle=True, encoding='latin1').item()
adjustedPValues = np.load(outFilePrefix + 'adjustedPValues.npy', allow_pickle=True, encoding='latin1').item()

scaledP = dict()
for svType in svTypes:

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

#colors = ['blue', 'red', 'magenta', 'black']
colors = ['magenta', 'black']
#offset = [-0.02, -0.01, 0.01, 0.02]
offset = [0.01, 0.02]
ind = 0
for svType in svTypes:
	print('area')
	area = 4 * (1 + (-np.log(adjustedPValues[svType])))
	print(area)
	ax.scatter(centers+offset[ind], scaledP[svType], color=colors[ind], alpha=0.3, s=area)
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
plt.show()