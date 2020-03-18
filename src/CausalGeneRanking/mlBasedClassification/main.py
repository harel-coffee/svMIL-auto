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
import random
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

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

import gseapy
import pandas as pd

#testing with polar plots

# #e.g. we have our bar charts
# barData = np.array([10, 2, 5, 6, 8, -20, 30, 15, 14, 20, 1,2, 3,4, 5,6, 4,1, 10, 2, 5, 6, 8, -20, 30, 15, 14, 20, 1,2, 3,4, 5,6, 4,1])
# blockSize = 360 / len(barData)
# print(blockSize)
# xRange = np.arange(0, 360, blockSize)
# print(xRange)
# #normal plot
# #plt.bar(xRange, barData)
# #plt.show()
# 
# degrees = np.random.randint(0, 360, size=200)
# 
# bin_size = 20
# 
# a , b=np.histogram(degrees, bins=np.arange(0, 360+bin_size, bin_size))
# centers = np.deg2rad(np.ediff1d(b)//2 + b[:-1])
# 
# print(a.shape)
# print(b.shape)
# print(centers.shape)
# 
# a = barData
# b = xRange
# #centers = np.deg2rad(np.ediff1d(b)//2 + b[:-1])
# 
# b = np.append(xRange, xRange[xRange.shape[0]-1]+blockSize)
# centers = np.deg2rad(np.ediff1d(b)//2 + b[:-1])
# print(np.ediff1d(b)//2 + b[:-1])
# 
# print('features: ', a.shape)
# print('range: ', b.shape)
# print('centers: ', centers.shape)
# print(centers)
# 
# 
# # # fig = plt.figure(figsize=(10,8))
# # # ax = fig.add_subplot(111, projection='polar')
# # # ax.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, color='.8', edgecolor='k')
# # # ax.set_xticks(centers)
# # # ax.set_xticklabels(barData)
# # # ax.set_yticklabels([])
# # # ax.set_theta_zero_location("N")
# # # ax.set_theta_direction(-1)
# # # plt.show()
# # 
# #instead of bars, plot a line
# #the y position is now the bar data
# fig = plt.figure(figsize=(10,8))
# ax = fig.add_subplot(111, projection='polar')
# #ax.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, color='.8', edgecolor='k')
# area = 0.25 * barData**2
# #area2 = 0.25 * barData2**2
# ax.scatter(centers-0.05, barData, color='blue', alpha=0.3, s=area)
# #ax.scatter(centers, barData2, color='red', alpha=0.3, s=area2)
# ax.set_xticks(centers)
# ax.set_xticklabels(barData)
# ax.set_yticks([-40,0,40])
# ax.set_yticklabels([])
# ax.set_theta_zero_location("N")
# ax.set_theta_direction(-1)
# plt.show()



#settings for running in different scenarios
svTypes = [sys.argv[3]]
svTypes = ['DEL', 'DUP', 'INV', 'ITX']
normalize = False #re-normalize, or use the saved file for speed? 
optimize = False #optimize classifier? 
test = False #test classifier performance with CV? 
featureImportance = True
featureLoad = True #re-create feature importances or just load from file? 

adjustedPValues = dict()
allFeatureZScores = dict()

degPairs = np.loadtxt(sys.argv[2], dtype='object') #labels

if featureLoad == False:
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


	
	for svType in svTypes:
		
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
			
			print('Number of positive instances: ', len(positiveInstanceLabels))
			
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
		
		if featureImportance == True and featureLoad == False:
			
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
			# plt.bar(range(similarityMatrix.shape[1])[0:1000], importances[indices[0:1000]],
			# 	   color="r", align="center")
			# plt.xticks(range(similarityMatrix.shape[1])[0:1000][::25], range(0,1000, 25), rotation=90)
			# plt.tight_layout()
			# plt.savefig('feature_importances_top1000_' + svType + '.svg')
			# plt.show()
			# 
			# exit()
			
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
		
			#reomove the old type column
			newInstances = np.delete(newInstances, 33, 1)
		
			avgInstances = np.sum(newInstances, axis=0)
		
			totalInstances = avgInstances / newInstances.shape[0]
			
			print(totalInstances)
			
			xlabels = ['loss', 'gain', 'cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
					   'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
					   'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'cosmic', 'enhancerType', 'promoterType', 'eQTLType', 'superEnhancerType']
			
			# plt.bar(range(len(totalInstances)), totalInstances)
			# plt.xticks(range(len(totalInstances)), xlabels, rotation=90)
			# plt.tight_layout()
			# #plt.show()
			# plt.clf()
			
			
			#get the percentage in the top 100.
			#then, repeat this 100 times, randomly sampling 100 instances.
			
			#do a t-test, and compute a p-value.
			instanceCount = 100 #top X to check
			
			#get the original top X features
			topInstances = newInstances[indices[0:instanceCount]]
			
			#split into up/down regulation
			#filter the instances by up/downregulation if needed.
			#first get the bag pair name for the instance
			filteredInstances = []
			for instanceInd in range(0, topInstances.shape[0]):
				
				bagLabel = bagPairLabels[bagMap[instanceInd]]
				splitPair = bagLabel.split('_')
				
				shortPair = splitPair[7] + '_' + splitPair[0]
				
				#get z-score
				degPairInfo = degPairs[degPairs[:,0] == shortPair][0]
		
				#if the z-score matches this criterion, the SV-gene pair is positive
				if float(degPairInfo[5]) > 1.5:
					filteredInstances.append(topInstances[instanceInd])
			
			filteredInstances = np.array(filteredInstances)
			print(filteredInstances)
			#compute the percentages in these top X instances
			avgInstances = np.sum(filteredInstances, axis=0)
		
			totalInstances = avgInstances / filteredInstances.shape[0]
			
			print(totalInstances)
			
			#100 times, randomly sample
			#per feature, have a distribution
			nullDistributions = dict()
			for i in range(0,100):
				
				if i == 0:
					for featureInd in range(0, len(totalInstances)):
						nullDistributions[featureInd] = []
				
				#sample as much random instances as in our filtered instances
				randomIndices = random.sample(range(0,newInstances.shape[0]), filteredInstances.shape[0])
			
				randomTopInstances = newInstances[randomIndices]
				
				#compute the percentages in these top X instances
				avgRandomInstances = np.sum(randomTopInstances, axis=0)
			
				totalRandomInstances = avgRandomInstances / randomTopInstances.shape[0]
				
				for featureInd in range(0, len(totalRandomInstances)):
					nullDistributions[featureInd].append(totalRandomInstances[featureInd])
			
		
			#for each feature, compute a z-score
			featurePValues = []
			featureZScores = []
			for featureInd in range(0, len(nullDistributions)):
				
				if np.mean(nullDistributions[featureInd]) == 0 or np.std(nullDistributions[featureInd]) == 0:
					z = 0
					pValue = 1
					featureZScores.append(z)
					featurePValues.append(pValue)
					continue
				
				z = (totalInstances[featureInd] - np.mean(nullDistributions[featureInd])) / float(np.std(nullDistributions[featureInd]))
				pValue = stats.norm.sf(abs(z))*2
				
				featureZScores.append(z)
				featurePValues.append(pValue)
			
			#then get a p-value
			print(featureZScores)
			print(featurePValues)
			#do MTC on the p-values
			
			reject, pAdjusted, _, _ = multipletests(featurePValues, method='bonferroni')
		
			print(reject)
			print(pAdjusted)
			
			allFeatureZScores[svType] = featureZScores
			adjustedPValues[svType] = pAdjusted
	
		np.save('featureZScores.npy', allFeatureZScores)
		np.save('adjustedPValues.npy', adjustedPValues)
			
if featureImportance == True:
	
	allFeatureZScores = np.load('featureZScores.npy', allow_pickle=True, encoding='latin1').item()
	adjustedPValues = np.load('adjustedPValues.npy', allow_pickle=True, encoding='latin1').item()
	
	print(allFeatureZScores)
	print(adjustedPValues)
	
	xlabels = ['loss', 'gain', 'cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
			   'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
			   'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'cosmic', 'enhancerType', 'promoterType', 'eQTLType', 'superEnhancerType']
	
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

	print(overallMin)
	print(overallMax)

	scaledP = dict()
	for svType in svTypes:
		
		pAdjusted = adjustedPValues[svType]
		
		directionalAdjustedP = -np.log(pAdjusted) * np.sign(allFeatureZScores[svType])
		print(directionalAdjustedP)
		directionalAdjustedP += np.abs(overallMin) + 50
		
		print(directionalAdjustedP)
		
		if len(np.where(pAdjusted == 1)) > 0:
			zeroOffsetInd = np.where(pAdjusted == 1)[0][0]
			print(zeroOffsetInd)
			zeroOffset = directionalAdjustedP[zeroOffsetInd]
			print(zeroOffset)
		
		scaledP[svType] = directionalAdjustedP
	
	print(zeroOffset)
	
	border = zeroOffset

	signBorderTop = -np.log(0.05) + zeroOffset
	signBorderBottom = border - (signBorderTop - border)
	
	print(border)
	print(signBorderTop)
	print(signBorderBottom)
	
	print(scaledP)
	

	###try making polar scatter plots for the features
	blockSize = 360 / len(directionalAdjustedP)
	xRange = np.arange(0, 360, blockSize)
	
	#add the last value, which is missed in the range
	xRange = np.append(xRange, xRange[xRange.shape[0]-1]+blockSize)
	centers = np.deg2rad(np.ediff1d(xRange)//2 + xRange[:-1])
	
	
	#instead of bars, plot a line
	#the y position is now the bar data
	fig = plt.figure(figsize=(15,13))
	ax = fig.add_subplot(111, projection='polar')
	#ax.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, color='.8', edgecolor='k')
	#area = 0.005 * directionalAdjustedP**2
	
	#exit()
	colors = ['blue', 'red', 'magenta', 'black']
	offset = [-0.02, -0.01, 0.01, 0.02]
	ind = 0
	for svType in svTypes:
		print('area')
		area = 4 * (1 + (-np.log(adjustedPValues[svType])))
		print(area)
		ax.scatter(centers+offset[ind], scaledP[svType], color=colors[ind], alpha=0.3, s=area)
		ind += 1

	ax.set_xticks(centers)
	ax.set_xticklabels(xlabels, fontsize=5)
	#ax.set_yticks([-1.5,np.log(0.05),0,-np.log(0.05), 1.5])
	#ax.set_yticks([-1.5,-0.05,0,0.05, 1.5])
	#ax.set_yticklabels(['', 'P < 0.05', '', 'P < 0.05', ''])
	
	ax.set_yticklabels(['p < 0.05'])
	ax.set_yticks([signBorderBottom, border, signBorderTop])
	gridlines = ax.yaxis.get_gridlines()
	gridlines[0].set_color("red")
	gridlines[0].set_linewidth(0.5)
	gridlines[0].set_linestyle('--')
	
	gridlines[2].set_color("red")
	gridlines[2].set_linewidth(0.5)
	gridlines[2].set_linestyle('--')
	
	# plt.gcf().canvas.draw()
	# angles = np.linspace(0,2*np.pi,len(ax.get_xticklabels())+1)
	# angles[np.cos(angles) < 0] = angles[np.cos(angles) < 0] + np.pi
	# angles = np.rad2deg(angles)
	# labels = []
	# for label, angle in zip(ax.get_xticklabels(), angles):
	# 	x,y = label.get_position()
	# 	lab = ax.text(x,y, label.get_text(), transform=label.get_transform(),
	# 				  ha=label.get_ha(), va=label.get_va())
	# 	#lab.set_rotation(90)
	# 	labels.append(lab)
	# ax.set_xticklabels([])
	
	
	
	#ax.set_theta_zero_location("N")
	#ax.set_theta_direction(-1)
	ax.set_rorigin(-300)
	ax.set_theta_zero_location('N', offset=10)
	plt.savefig('featureImportances_allTypes.svg')
	plt.show()

		
		
		