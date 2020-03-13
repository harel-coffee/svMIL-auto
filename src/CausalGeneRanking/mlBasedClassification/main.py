"""
	Interface for Multiple Instance Learning on SV-gene pairs to classify these as having transcriptional consequences or not. 
	
	To do:
	- move parts to separate functions
	- move things around to settings
	- feature elimination
	- feature importance to learn more about biology

"""

import numpy as np
import sys
from inputParser import InputParser
from featureMatrixDefiner import FeatureMatrixDefiner
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import math
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import Lasso
from sklearn.metrics import auc, precision_recall_curve
from sklearn.svm import LinearSVC
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from scipy import interp
import pickle as pkl
from random import shuffle

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from inspect import signature
from sklearn import model_selection
from sklearn.metrics import average_precision_score
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import train_test_split

#settings for running in different scenarios
svType = sys.argv[3]
normalize = False #re-normalize, or use the saved file for speed? 
optimize = False #optimize classifier? 
test = False #test classifier performance with CV? 
featureImportance = True

degPairs = np.loadtxt(sys.argv[2], dtype='object') #labels

if normalize == True:
	
	#Get the bags
	with open(sys.argv[1], 'rb') as handle:
		bagDict = pkl.load(handle)
	
	#first check what the minimum and maximum feature values are in order to normalize
	currentMax = [0]*35 #would be better if we read this from the input file, this is how many features we have for each instance
	currentMin = [float('inf')]*35
	for pair in bagDict:
		
		for instance in bagDict[pair]:
	
			if instance[0] == '0' and instance[1] == '0': #skip instances where there are no gains or losses, if these slipped in somehow
				continue
			
			for featureInd in range(0, len(instance)):
				feature = instance[featureInd]
				if feature > currentMax[featureInd]:
					currentMax[featureInd] = feature
				if feature < currentMin[featureInd]:
					currentMin[featureInd] = feature
	
	#loop through each bag, and normalize the features of each instance based on the max/min 
	normalizedBagDict = dict()
	for pair in bagDict:
		
		normalizedBagDict[pair] = []
		
		for instance in bagDict[pair]:
			
			normInstance = []
			
			for featureInd in range(0, len(instance)):
				
				
				feature = instance[featureInd]
				
				if currentMin[featureInd] == 0 and currentMax[featureInd] == 0: #if the min/max are 0 for this feature, the normalized value should also be 0. 
					normInstance.append(0)
					continue 
				
				#do the normalization
				normFeature = (feature-currentMin[featureInd])/(currentMax[featureInd]-currentMin[featureInd])
				normInstance.append(normFeature)
				
			normalizedBagDict[pair].append(normInstance)

	bagDict = normalizedBagDict
	
	#save to a file to prevent heavy computing load
	with open('normalizedBags.pkl', 'wb') as handle:
			pkl.dump(bagDict, handle, protocol=pkl.HIGHEST_PROTOCOL)
else:
	with open('normalizedBags.pkl', 'rb') as handle:
		bagDict = pkl.load(handle)

print("initial number of bags: ", len(bagDict))
print('initial deg pairs: ', degPairs.shape[0])

#allow for running with feature selection
featureCount = len(bagDict[list(bagDict.keys())[0]][0])
featureStart = featureCount #set this to featureCount to run with all features. (make setting later)
similarityMatrices = dict() #store the similarity matrices for each feature selection run
bagLabels = []
positiveBagPairNames = []
negativeBagPairNames = []
positiveInstanceLabels = []
for featureInd in range(featureStart, featureCount+1):
	
	positiveBags = []
	negativeBags = []
	
	#for each SV-gene pair, get the instances
	for pair in bagDict:
		
		#check if the SV type matches our selection
		splitPair = pair.split("_")
		shortPair = splitPair[7] + '_' + splitPair[0]
		
		if svType != '':	
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
					positiveInstanceLabels.append(pair + '_' + '_'.join([str(i) for i in instance]))
					instances.append(instance[0:featureInd])

					
				if len(instances) < 1:
					continue
				
				positiveBagPairNames.append(pair)
				positiveBags.append(instances)
				
			else: #if the z-score is anything else, this bag will be labeled negative. 
				
				#get the right number of features per instance
				instances = []
				for instance in bagDict[pair]:
					
					if instance[0] == 0 and instance[1] == 0:
						continue

					instances.append(instance[0:featureInd])
				
				if len(instances) < 1:
					continue

				negativeBags.append(instances)
				negativeBagPairNames.append(pair)
				
	
	positiveBags = np.array(positiveBags)
	negativeBags = np.array(negativeBags)
	positiveBagPairNames = np.array(positiveBagPairNames)
	negativeBagPairNames = np.array(negativeBagPairNames)
	
	print('Number of positive bags: ', positiveBags.shape)
	print('Number of negative bags: ', negativeBags.shape)
	
	#set a random seed to always subsample the same set
	np.random.seed(0)
	#subsample the negative set to the same number of positives. 
	negativeBagsSubsampled = np.random.choice(negativeBags, positiveBags.shape[0])
	
	negativeBagsSubsampleInd = np.random.choice(np.arange(negativeBags.shape[0]), positiveBags.shape[0])
	negativeBagsSubsampled = negativeBags[negativeBagsSubsampleInd]
	
	negativeBagPairNamesSubsampled = negativeBagPairNames[negativeBagsSubsampleInd]
	bagPairLabels = np.concatenate((positiveBagPairNames, negativeBagPairNamesSubsampled))

	#merge the bags so that we can easily get to 1 similarity matrix and do all-to-all computations
	bags = np.concatenate((positiveBags, negativeBagsSubsampled))
	#assign bag labels
	bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBagsSubsampled.shape[0])
	
	#stack the instances in the bags so that we can easily compute bag-instance distances
	instances = np.vstack(bags)

	#Make similarity matrix	
	print("generating similarity matrix")
	
	#Make an index where we can lookup at which position the instances are in the concatenated bag array. 
	reverseBagMap = dict() #lookup instance by bag index
	bagMap = dict() #lookup bag by instance index
	instanceInd = 0
	for bagInd in range(0, bags.shape[0]):
		reverseBagMap[bagInd] = []
		for instance in bags[bagInd]:
			reverseBagMap[bagInd].append(instanceInd)
			bagMap[instanceInd] = bagInd
			
			instanceInd += 1
			
	
	bagIndices = np.arange(bags.shape[0])
	similarityMatrix = np.zeros([bags.shape[0], instances.shape[0]])
	print("Number of bags: ", bags.shape[0])
	for bagInd in range(0, bags.shape[0]):
		
		#Get the indices of the instances that are in this bag
		instanceIndices = reverseBagMap[bagInd]		
		instanceSubset = instances[instanceIndices,:]

		#get the average of all instances in this bag
		instanceAvg = np.mean(instanceSubset, axis=0)
		
		#compute distance to all other instances from this bag average
		distance = np.abs(instanceAvg - instances)
		
		#sum the distances to get 1 similarity score
		summedDistance = np.sum(distance,axis=1)
		similarityMatrix[bagInd,:] = summedDistance
		
	print('similarity matrix: ')
	print(similarityMatrix)
	similarityMatrices[featureInd] = similarityMatrix



def cvClassification(similarityMatrix, bagLabels, clf, svType, title):
	
	#get the kfold model
	kfold = model_selection.StratifiedKFold(n_splits=10, shuffle=True, random_state=10)
	
	#make the ROC curve and compute the AUC
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)
	
	fig, ax = plt.subplots()
	importances = []
	for i, (train, test) in enumerate(kfold.split(similarityMatrix, bagLabels)):
		clf.fit(similarityMatrix[train], bagLabels[train])
		viz = plot_roc_curve(clf, similarityMatrix[test], bagLabels[test],
							 name='ROC fold {}'.format(i),
							 alpha=0.3, lw=1, ax=ax)
		interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(viz.roc_auc)
		importances.append(clf.feature_importances_)
	print('aucs: ')
	print(aucs)
	print('mean auc: ', np.mean(aucs))
	print('std of auc: ', np.std(aucs))

	ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
			label='Chance', alpha=.8)
	
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	
	ax.plot(mean_fpr, mean_tpr, color='b',
			label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (np.mean(aucs), np.std(aucs)),
			lw=2, alpha=.8)
	
	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
					label=r'$\pm$ 1 std. dev.')
	
	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
		   title="Receiver operating characteristic: " + title)
	ax.legend(loc="lower right")
	plt.tight_layout()
	plt.savefig('miles_' + svType + '.svg')
	plt.show()


#do RF optimization with random parameter search
if optimize == True:
	#Number of trees in random forest
	n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
	# Number of features to consider at every split
	max_features = ['auto', 'sqrt']
	# Maximum number of levels in tree
	max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
	max_depth.append(None)
	# Minimum number of samples required to split a node
	min_samples_split = [2, 5, 10]
	# Minimum number of samples required at each leaf node
	min_samples_leaf = [1, 2, 4]
	# Method of selecting samples for training each tree
	bootstrap = [True, False]# Create the random grid
	random_grid = {'n_estimators': n_estimators,
				   'max_features': max_features,
				   'max_depth': max_depth,
				   'min_samples_split': min_samples_split,
				   'min_samples_leaf': min_samples_leaf,
				   'bootstrap': bootstrap}
	
	#initial classifier to optimize from
	rfClassifier = RandomForestClassifier(n_estimators= 500)
	rf_random = RandomizedSearchCV(estimator = rfClassifier, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)
	
	X_train, X_test, y_train, y_test = train_test_split(similarityMatrices[featureStart], bagLabels, test_size=0.33, random_state=42)
	rf_random.fit(X_train, y_train)
	print('best params; ')
	print(rf_random.best_params_)
	print('new score: ')
	print(rf_random.score(X_test, y_test))
	
	print('base score :')
	rfClassifier.fit(X_train, y_train)
	print(rfClassifier.score(X_test, y_test))
	exit()
	
	
#some results of optimization: 
#best DEL: n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True
#best ITX: n_estimators= 1000, min_samples_split=2, min_samples_leaf=1, max_features='sqrt', max_depth=110, bootstrap=True
#best inv: n_estimators= 400, min_samples_split=10, min_samples_leaf=2, max_features='auto', max_depth=40, bootstrap=False
#best DUP: n_estimators= 1200, min_samples_split=2, min_samples_leaf=4, max_features='sqrt', max_depth=70, bootstrap=False

#run2

#best itx: n_estimators= 1000, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True
#best inv: {'n_estimators': 200, 'min_samples_split': 5, 'min_samples_leaf': 4, 'max_features': 'auto', 'max_depth': 10, 'bootstrap': True}
#best dup: 'n_estimators': 600, 'min_samples_split': 2, 'min_samples_leaf': 2, 'max_features': 'sqrt', 'max_depth': 110, 'bootstrap': False
#best del:

#test classifier performance
if test == True:
	
	if svType == 'DEL':
		classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'deletions'
	elif svType == 'DUP':
		#classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=2, min_samples_leaf=2, max_features='sqrt', max_depth=110, bootstrap=False)
		classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		#classifier = RandomForestClassifier(n_estimators= 1200, min_samples_split=2, min_samples_leaf=4, max_features='sqrt', max_depth=70, bootstrap=False)
		title = 'duplications'
	elif svType == 'INV':
		classifier = RandomForestClassifier(n_estimators= 200, min_samples_split=5, min_samples_leaf=4, max_features='auto', max_depth=10, bootstrap=True)
		title = 'inversions'
	elif svType == 'ITX':
		classifier = RandomForestClassifier(n_estimators= 1000, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'translocations'
	else:
		classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'All SV types'

	#for each feature selection round, show the performance
	for featureInd in range(featureStart, featureCount+1):
		
		cvClassification(similarityMatrices[featureInd], bagLabels, classifier, svType, title)
		
	#repeat, but then with random labels.
	shuffle(bagLabels)
	for featureInd in range(featureStart, featureCount+1):
		cvClassification(similarityMatrices[featureInd], bagLabels, classifier, svType, title)

if featureImportance == True:
	
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
	
	
	clf.fit(similarityMatrix, bagLabels)
	print(clf.score(similarityMatrix, bagLabels))
	importances = clf.feature_importances_
	std = np.std([tree.feature_importances_ for tree in clf.estimators_],axis=0)
	
	#rank these importances, get the first 100, see which instances these are. 
	indices = np.argsort(importances)[::-1]
	
	# plt.figure()
	# plt.title("Feature importances")
	# plt.bar(range(similarityMatrix.shape[1]), importances[indices],
	# 	   color="r", align="center")
	# plt.xticks(range(similarityMatrix.shape[1]), indices)
	# plt.xlim([-1, similarityMatrix.shape[1]])
	# plt.show()
	
	#get the features of the top 100 ranked instances.
	#for each feature, how many % of the top 100 has this feature?
	#topInstances = instances[indices[0:100]]
	#topInstances = instances[indices]
	
	#split the type back into 4 features
	enhancerTypes = []
	eQTLTypes = []
	promoterTypes = []
	superEnhancerTypes = []
	for instance in instances:
		
		if instance[33] == 0:
			enhancerTypes.append(1)
			eQTLTypes.append(0)
			promoterTypes.append(0)
			superEnhancerTypes.append(0)
		elif instance[33] > 0 and instance[33] < 0.34:
			enhancerTypes.append(0)
			eQTLTypes.append(0)
			promoterTypes.append(1)
			superEnhancerTypes.append(0)
		elif instance[33] > 0.33 and instance[33] < 0.68:
			enhancerTypes.append(0)
			eQTLTypes.append(1)
			promoterTypes.append(0)
			superEnhancerTypes.append(0)
		else:
			enhancerTypes.append(0)
			eQTLTypes.append(0)
			promoterTypes.append(0)
			superEnhancerTypes.append(1)
	
	newInstances = np.zeros([instances.shape[0], instances.shape[1]+4])
	
	newInstances[:,0:instances.shape[1]] = instances
	
	newInstances[:,instances.shape[1]] = enhancerTypes
	newInstances[:,instances.shape[1]+1] = promoterTypes
	newInstances[:,instances.shape[1]+2] = eQTLTypes
	newInstances[:,instances.shape[1]+3] = superEnhancerTypes
	
	###only here do the subsetting
	
	avgInstances = np.sum(newInstances, axis=0)

	totalInstances = avgInstances / newInstances.shape[0]
	
	print(totalInstances)
	
	xlabels = ['loss', 'gain', 'cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
			   'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
			   'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'type', 'cosmic', 'enhancerType', 'promoterType', 'eQTLType', 'superEnhancerType']
	
	plt.bar(range(len(totalInstances)), totalInstances)
	plt.xticks(range(len(totalInstances)), xlabels, rotation=90)
	plt.tight_layout()
	plt.show()
	
	#output:
	#1. file with genesets, which are all the positive instances
	
	#2. file with the ranked values, .rnk file
	
	#get the names of the instances
	
	#first, we need to know in which bag the instance is to get the label of the 
	#so based on the index of the instance, get the bag label. 
	#use the bag map for this.
	
	#e.g. get the bag index of the second instance
	
	#merge this with the features of this instance to get a label for the instance.
	instancesWithValues = []
	
	for instanceInd in range(0, newInstances.shape[0]):
		instance = newInstances[instanceInd]
		
		#label ths instance
		bagLabel = bagPairLabels[bagMap[instanceInd]]
		instanceLabel = bagLabel + '_' + '_'.join([str(i) for i in instance])
		
		#get the feature importance for this instance
		#first test with 1 feature only, say gains
		instancesWithValues.append([instanceLabel, instance[1]])
		
	#rank the labeled instances by feature importances
	instancesWithValues = np.array(instancesWithValues, dtype='object')
	
	rankedLabeledInstances = instancesWithValues[indices]
	
	print(rankedLabeledInstances)
	
	np.savetxt('rankedInstances_gains.rnk', rankedLabeledInstances, fmt='%s', delimiter='\t')
	np.savetxt('positiveInstances.grp', positiveInstanceLabels, fmt='%s')
	