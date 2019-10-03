"""

	Make features matrices for all features in a 1 MB window.
	Try training a deep learning model

"""

import numpy as np
import sys
from inputParser import InputParser
from featureMatrixDefiner import FeatureMatrixDefiner
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import math
from tsne import bh_sne
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import Lasso
from sklearn.metrics import auc, precision_recall_curve
from sklearn.svm import LinearSVC
from cleanlab.classification import LearningWithNoisyLabels
from cleanlab.noise_generation import generate_noise_matrix_from_trace
from cleanlab.noise_generation import generate_noisy_labels
from cleanlab.util import print_noise_matrix
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from scipy import interp

#1. Process the data
#Use the input parser to parse the SV-gene pairs and return labels
pairs, labels = InputParser().processSVGenePairs(sys.argv[1], sys.argv[2])

#2. Define a feature matrix
featureMatrixDefiner = FeatureMatrixDefiner()
featureMatrixDefiner.setFeatureData()
svGenePairsRules = np.loadtxt(sys.argv[3], dtype='object')
bags, instances, bagLabels = featureMatrixDefiner.defineFeatureMatrix(pairs, labels, svGenePairsRules)

print(bags.shape)
print(instances.shape)
print(bagLabels.shape)

#Make similarity matrix

print("generating similarity matrix")

#Unfold the training bags so that we can compute the distance matrix at once to all genes
bagMap = dict()
reverseBagMap = dict()
geneInd = 0
for bagInd in range(0, bags.shape[0]):
	reverseBagMap[bagInd] = []
	for gene in bags[bagInd]:
		bagMap[geneInd] = bagInd
		reverseBagMap[bagInd].append(geneInd)
		
		geneInd += 1

bagIndices = np.arange(bags.shape[0])
similarityMatrix = np.zeros([bags.shape[0], instances.shape[0]])
print("Number of bags: ", bags.shape[0])
for bagInd in range(0, bags.shape[0]):
	
	#Get the indices of the instances that are in this bag
	instanceIndices = reverseBagMap[bagInd]
	
	instanceSubset = instances[instanceIndices,:]
	otherInstances = np.vstack(bags[bagIndices != bagInd])
	
	instanceAvg = [np.mean(instanceSubset[:,0]), np.mean(instanceSubset[:,1])]
	
	#compute distance to all other instances
	distance = np.abs(instanceAvg - instances)

	summedDistance = np.sum(distance,axis=1)
	similarityMatrix[bagInd,:] = summedDistance
	continue
	
	
	
	#Compute the pairwise distance matrix here
	minDistance = float("inf")
	minDistanceInd = 0
	for instanceInd in range(0, instanceSubset.shape[0]):
		instance = instanceSubset[instanceInd]
		distance = np.abs(instance - instances) #compute the distances to the train instances, otherwise we are not in the same similarity space. 

		#distance = np.abs(instance - otherInstances)

		summedDistance = np.sum(distance,axis=1)

		currentMinDistance = np.mean(summedDistance)
		if currentMinDistance < np.mean(minDistance):
			minDistance = summedDistance
			minDistanceInd = instanceInd

	#This instance will be used as representative for this bag. We use this value as the similarity to all other instances.  
	similarityMatrix[bagInd] = minDistance

print(similarityMatrix)
print(bagLabels)
# 
# plt.scatter(similarityMatrix[0,:], similarityMatrix[1,:])
# plt.show()
# exit()


# pca = PCA(n_components=2)
# projected = pca.fit_transform(similarityMatrix)
# 
# colorLabels = []
# 
# for label in bagLabels:
# 	
# 	if label == 1:
# 		colorLabels.append('r')
# 	else:
# 		colorLabels.append('b')
# 
# fig,ax=plt.subplots(figsize=(7,5))
# plt.scatter(projected[:, 0], projected[:, 1], c=colorLabels)
# plt.show()

from random import shuffle
shuffle(bagLabels)

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits=10)
np.random.seed(500)

accs = []
aucs = []
coeffs = []
predDiffs = []
for train, test in cv.split(similarityMatrix, bagLabels):
	
	rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
	rfClassifier.fit(similarityMatrix[train], bagLabels[train]) #Use the bag labels, not the instance labels

	predictions = rfClassifier.predict(similarityMatrix[test])
	precision, recall, thresholds = precision_recall_curve(bagLabels[test], predictions)
	aucScore = auc(recall, precision)
	predsDiff = np.average(bagLabels[test] == np.sign(predictions))
	#Now select the most important features with random forest
	importances = rfClassifier.feature_importances_
	std = np.std([tree.feature_importances_ for tree in rfClassifier.estimators_],
				 axis=0)
	indices = np.argsort(importances)[::-1]
	
	nonZeroIndices = []
	for index in indices:
		if importances[index] > 0:
			nonZeroIndices.append(index)
	
	aucs.append(aucScore)
	predDiffs.append(predsDiff)

print("Actual acc: ", np.mean(predDiffs))
print("Mean AUC: ", np.mean(aucs))
exit()
# 
# 
# 
# vis_data = bh_sne(matrix)
# 
# # plot the result
# vis_x = vis_data[:, 0]
# vis_y = vis_data[:, 1]
# 
# plt.scatter(vis_x, vis_y, c=colorLabels)
# plt.show()
# 
# exit()


#Implement leave-one-patient out CV for all classifiers, make a function
def loopoCV(sortedPairs, sortedLabels, allFeatures, clf, aucBool= True):
	
	#gather all unique patients
	patients = np.unique(sortedPairs[:,7])
	patientScores = []
	patientAucs = []
	patientAuprcs = []
	posZero = 0
	for patient in patients:
		#gather a training set for all but this patient
		X_train = allFeatures[sortedPairs[:,7] != patient]
		y_train = sortedLabels[sortedPairs[:,7] != patient]
		
		X_test = allFeatures[sortedPairs[:,7] == patient]
		y_test = sortedLabels[sortedPairs[:,7] == patient]
		
		trainFeatPos = X_train[y_train == 1]
		testFeatPos = X_test[y_test == 1]
		
		if len(testFeatPos) == 0 or len(trainFeatPos) == 0:
			continue #not enough positive samples for classification
		
		
		#Train classifier on each fold.
		clf.fit(X_train, y_train)

		predictions = clf.predict(X_test)
		predsDiff = np.average(y_test == np.sign(predictions))
		patientScores.append(predsDiff)
		
		if aucBool == True: #not available for all classifiers
			preds = clf.predict_proba(X_test)[:,1]
			fpr, tpr, thresholds = metrics.roc_curve(y_test, preds, pos_label=1)
			aucS = metrics.auc(fpr, tpr)
			patientAucs.append(aucS)
		
		precision, recall, thresholds = precision_recall_curve(y_test, predictions)
		aucScore = auc(recall, precision)
		patientAuprcs.append(aucScore)
	
	#report averages and std
	print("Average CV score: ", np.mean(patientScores))
	print("std CV score: ", np.std(patientScores))
	if aucBool == True:
		print("Average CV AUC: ", np.mean(patientAucs))
		print("std CV AUC: ", np.std(patientAucs))
		
	print("Average CV AUPRC: ", np.mean(patientAuprcs))
	print("std CV AUPRC: ", np.std(patientAuprcs))


def performCV(featureMatrix, labels, clf):
	
	cv = StratifiedKFold(n_splits=10)
	
	#dirty
	X = featureMatrix
	y = np.array(labels)
	
	tprs = []
	aucs = []
	scores = []
	mean_fpr = np.linspace(0, 1, 100)
	
	i = 0
	for train, test in cv.split(X, y):
		clf.fit(X[train], y[train])
		
		predictions = clf.predict(X[test])
		score = np.average(y[test] == np.sign(predictions))
		scores.append(score)
		
		probas_ = clf.predict_proba(X[test])
		
		# Compute ROC curve and area the curve
		fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
		tprs.append(interp(mean_fpr, fpr, tpr))
		tprs[-1][0] = 0.0
		roc_auc = auc(fpr, tpr)
		aucs.append(roc_auc)
		
		i += 1
	
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	mean_score = np.mean(scores)
	
	print("Score: ", mean_score)
	print("AUC: ", mean_auc)
	
	
	

#Try each classifier
print("Random forest")
from sklearn.ensemble import RandomForestClassifier
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
#loopoCV(pairs, labels, matrix, rfClassifier, True)
performCV(matrix, labels, rfClassifier)

print("lasso")

currentAlpha = 1e-2
lasso = Lasso(alpha=currentAlpha)
#loopoCV(pairs, labels, matrix, lasso, False)
#performCV(matrix, labels, lasso)
print("linear SVC")

#clf = LinearSVC()
#loopoCV(pairs, labels, matrix, clf, False)

print("cleaned labels RF")
rp = LearningWithNoisyLabels(rfClassifier)
#loopoCV(pairs, labels, matrix, rp, False)
performCV(matrix, labels, rp)

#print("cleaned labels lasso")
#rp = LearningWithNoisyLabels(lasso)
#loopoCV(sortedPairs, sortedLabels, allFeatures, rp, False)
# 
# print("cleaned labels SVC")
# rp = LearningWithNoisyLabels(clf)
# loopoCV(pairs, labels, matrix, rp, False)

exit()

# 
# 
# from sklearn.model_selection import train_test_split
# from sklearn import metrics
# from sklearn.metrics import roc_curve, auc
# from sklearn.linear_model import Lasso
# from sklearn.metrics import auc, precision_recall_curve
# from sklearn.svm import LinearSVC
# 
# #Randomize labels
# from random import shuffle
# #shuffle(labels)
# 
# print labels
# 
# X_train, X_test, y_train, y_test = train_test_split(matrix, labels, test_size=0.4, random_state=42)
# 
# #1. Random forest
# from sklearn.ensemble import RandomForestClassifier
# rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
# rfClassifier.fit(X_train, y_train) #Use the bag labels, not the instance labels
# 
# predictions = rfClassifier.predict(X_test)
# predsDiff = np.average(y_test == np.sign(predictions))
# print("RF score: ", predsDiff)
# 
# precision, recall, thresholds = precision_recall_curve(y_test, predictions)
# aucScore = auc(recall, precision)
# print("RF AUC: ", aucScore)
# 
# #2. Lasso
# 
# currentAlpha = 1e-2
# lasso = Lasso(alpha=currentAlpha)
# lasso.fit(X_train,y_train)
# 
# test_score=lasso.score(X_test,y_test)
# coeff_used = np.sum(lasso.coef_!=0)
# preds = lasso.predict(X_test)
# predsDiff = np.average(y_test == np.sign(preds))
# print("lasso score: ", predsDiff)
# 
# precision, recall, thresholds = precision_recall_curve(y_test, preds)
# aucScore = auc(recall, precision)
# print("lasso AUC: ", aucScore)
# 
# #3. SVM
# 
# 
# clf = LinearSVC()
# clf.fit(X_train, y_train)
# score = clf.score(X_test, y_test)
# preds = clf.predict(X_test)
# predsDiff = np.average(y_test == np.sign(preds))
# print("SVM score: ", score)
# 
# precision, recall, thresholds = precision_recall_curve(y_test, preds)
# aucScore = auc(recall, precision)
# 
# print("SVM AUC: ", aucScore)
# 
# 
# 
# 
# 
# exit()
#3. Run deep learning model

from sklearn.model_selection import train_test_split
from mcfly import modelgen, find_architecture, storage
from keras.models import load_model
import os

#use just 1 channel for now, later split into 2 and see if it makes a difference
X_train, X_test, y_train_list, y_test_list = train_test_split(matrix, labels, test_size=.4, random_state=42)

#the training labels need to be a vector as well. For each gene we have a 1 or 0 for each class. We have 2 classes, so this will be genes * 2

y_train = np.zeros([len(y_train_list), 2])
for labelInd in range(0, len(y_train_list)):
	
	label = y_train_list[labelInd]
	
	if label == 1:
		y_train[labelInd, 0] = 0
		y_train[labelInd, 1] = 1
	if label == 0:
		y_train[labelInd, 0] = 1
		y_train[labelInd, 1] = 0

y_test = np.zeros([len(y_test_list), 2])
for labelInd in range(0, len(y_test_list)):
	
	label = y_test_list[labelInd]
	
	if label == 1:
		y_test[labelInd, 0] = 0
		y_test[labelInd, 1] = 1
	if label == 0:
		y_test[labelInd, 0] = 1
		y_test[labelInd, 1] = 0
	


num_classes = y_train.shape[1]
X_train = np.array(X_train)

X_test = np.array(X_test)


models = modelgen.generate_models(X_train.shape,
								  number_of_classes=num_classes,
								  number_of_models = 2)


models_to_print = list(range(len(models)))
for i, item in enumerate(models):
	if i in models_to_print:
		model, params, model_types = item
		print("-------------------------------------------------------------------------------------------------------")
		print(("Model " + str(i)))
		print(" ")
		print("Hyperparameters:")
		print(params)
		print(" ")
		print("Model description:")
		model.summary()
		print(" ")
		print("Model type:")
		print(model_types)
		print(" ")

# Define directory where the results, e.g. json file, will be stored
resultpath = os.path.join('.', 'models')
if not os.path.exists(resultpath):
		os.makedirs(resultpath)
		
outputfile = os.path.join(resultpath, 'modelcomparison.json')
histories, val_accuracies, val_losses = find_architecture.train_models_on_samples(X_train, y_train,
																		   X_test, y_test,
																		   models,nr_epochs=5,
																		   subset_size=300,
																		   verbose=True,
																		   outputfile=outputfile)
print(('Details of the training process were stored in ',outputfile))

