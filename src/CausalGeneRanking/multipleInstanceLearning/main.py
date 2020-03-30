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
from copy import deepcopy

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

import pandas as pd

"""

	TO DO: update this code into functions!!

"""

#settings for running in different scenarios
svTypes = [sys.argv[3]]
#svTypes = ['DEL', 'DUP', 'INV', 'ITX']
#svTypes = ['DEL', 'INV', 'ITX']
svTypes = ['ITX']
normalize = False #re-normalize, or use the saved file for speed? 
optimize = False #optimize classifier? 
test = True #test classifier performance with CV?
featureSelection = True #randomize features 1 by 1? 
featureImportance = False
featureLoad = False #re-create feature importances or just load from file? 

adjustedPValues = dict()
allFeatureZScores = dict()


if featureLoad == False:

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
			
				plot = True
				if featureSelection == True:
					plot = False
				
				cvClassification(similarityMatrix, bagLabels, classifier, svType, title, plot)
					
				#repeat, but then with random labels.
				#mind here, if multiple iterations, the baglabels are permanently shuffled!
				if featureSelection == False:
					shuffle(bagLabels)
				
					cvClassification(similarityMatrix, bagLabels, classifier, svType, title, plot)
			
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
			
			#read the hallmarks file
			hallmarksFile = 'h.all.v7.0.symbols.gmt'
			hallmarkGenes = []
			with open(hallmarksFile, 'r') as inF:
				
				for line in inF:
					line = line.strip()
					splitLine = line.split('\t')
					
					for gene in splitLine[2:len(splitLine)]:
						
						hallmarkGenes.append(gene)
			
			#split the type back into 4 features
			enhancerTypes = []
			eQTLTypes = []
			promoterTypes = []
			superEnhancerTypes = []
			uniqueGenes = []
			uniqueCosmicGenes = []
			instanceInd = 0
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
			
			#get the percentage in the top 100.
			#then, repeat this 100 times, randomly sampling 100 instances.
			
			#do a t-test, and compute a p-value.
			instanceCount = 100 #top X to check
			
			#get the original top X features
			topInstances = newInstances[indices[0:instanceCount]]
			
			avgInstances = np.sum(topInstances, axis=0)
			
			#split into up/down regulation
			#filter the instances by up/downregulation if needed.
			#first get the bag pair name for the instance
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
					
				if topInstances[instanceInd][33] > 0:
					cosmicGenes.append(splitPair[0])
					if splitPair[0] not in uniqueCosmicGenes:
						uniqueCosmicGenes.append(splitPair[0])
				
				#get z-score
				degPairInfo = degPairs[degPairs[:,0] == shortPair][0]
		
				#cosmic yes/no
				if topInstances[instanceInd,33] > 0:
					if topInstances[instanceInd,1] > 0:
						filteredInstances.append(topInstances[instanceInd])
						
						print(bagLabel)
						print(topInstances[instanceInd])
						
					#check for hallmark yes/no
					#if splitPair[0] in hallmarkGenes:
					#	hallmarkCount += 1
				
			exit()
			filteredInstances = np.array(filteredInstances)
			print('number of instances used: ', filteredInstances.shape)
			
			#remove the gains and losses features, we do not need those anymore.
			filteredInstances = np.delete(filteredInstances, 1, 1)
			filteredInstances = np.delete(filteredInstances, 0, 1)
			
			print(svType, ':')
			print('number of cosmic genes (unique): ', len(uniqueCosmicGenes))
			print('number of cosmic genes per instance: ', len(cosmicGenes))
			
			np.savetxt('pairLabels_top100_' + svType + '.txt', topPairLabels, fmt='%s')
			
			#np.savetxt('allGenesTop100_' + svType + '.txt', uniqueGenes, fmt='%s')
			
			#compute the percentages in these top X instances
			avgInstances = np.sum(filteredInstances, axis=0)
		
			totalInstances = avgInstances / filteredInstances.shape[0]
			print(totalInstances)
			
			#100 times, randomly sample
			#per feature, have a distribution
			nullDistributions = dict()
			randomHallmarkCounts = []
			for i in range(0,100):
				
				randomHallmarkCount = 0
				
				if i == 0:
					for featureInd in range(0, len(totalInstances)):
						nullDistributions[featureInd] = []
				
				#sample as much random instances as in our filtered instances
				randomIndices = random.sample(range(0,newInstances.shape[0]), filteredInstances.shape[0])
			
				randomTopInstances = newInstances[randomIndices]
				
				#here also skip gains/losses
				randomTopInstances = np.delete(randomTopInstances, 1, 1)
				randomTopInstances = np.delete(randomTopInstances, 0, 1)
				
				#compute the percentages in these top X instances
				avgRandomInstances = np.sum(randomTopInstances, axis=0)
			
				totalRandomInstances = avgRandomInstances / randomTopInstances.shape[0]
				
				for featureInd in range(0, len(totalRandomInstances)):
					nullDistributions[featureInd].append(totalRandomInstances[featureInd])
			
				#check for hallmark yes/no
				for instanceInd in randomIndices:
					bagLabel = bagPairLabels[bagMap[instanceInd]]
					splitPair = bagLabel.split('_')
					
					shortPair = splitPair[7] + '_' + splitPair[0]
					
					if splitPair[0] in hallmarkGenes:
						randomHallmarkCount += 1
		
				randomHallmarkCounts.append(randomHallmarkCount)
		
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

			#output the top 100 instances to a file as well.
			#obtain the instance of the bag 
			
			#check hallmark enrichment in true set vs randomly sampled
			z = (hallmarkCount - np.mean(randomHallmarkCounts) / float(np.std(randomHallmarkCounts)))
			pValue = stats.norm.sf(abs(z))*2
			
			print(svType, ':')
			print('hallmark stats: ', hallmarkCount, np.mean(randomHallmarkCounts), np.std(randomHallmarkCounts), z, pValue)
			
	
		np.save('featureZScores_loss.npy', allFeatureZScores)
		np.save('adjustedPValues_loss.npy', adjustedPValues)
			
if featureImportance == True:
	#load the full set so that we normalize properly
	allFeatureZScores = np.load('featureZScores_all.npy', allow_pickle=True, encoding='latin1').item()
	adjustedPValues = np.load('adjustedPValues_all.npy', allow_pickle=True, encoding='latin1').item()

	xlabels = ['cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
	 		   'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
	 		   'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'cosmic', 'enhancerType', 'promoterType', 'eQTLType', 'superEnhancerType']
	 
	#xlabels = ['cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
	#		   'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
	#		   'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'enhancerType', 'promoterType', 'eQTLType', 'superEnhancerType']
	
	
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

	
	allFeatureZScores = np.load('featureZScores_loss.npy', allow_pickle=True, encoding='latin1').item()
	adjustedPValues = np.load('adjustedPValues_loss.npy', allow_pickle=True, encoding='latin1').item()
	
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
	ax.set_rorigin(-200)
	ax.set_theta_zero_location('N', offset=10)
	plt.savefig('featureImportances_allTypes_loss.svg')
	plt.show()

		
		
		