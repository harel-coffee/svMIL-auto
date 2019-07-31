"""
	1. Read the gene scores
	2. Read the DEGs based on patients with and patients without SVs
	3. Try a random forest classifier
"""

import sys
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
# 
# scoreFile = sys.argv[1]
# degFile = sys.argv[2]
# 
# geneScores = np.loadtxt(scoreFile, dtype="object")
# 
# #Filter the gene scores to a feature matrix
# featureMatrix = []
# for gene in geneScores:
# 	
# 	samples = gene[len(gene)-1]
# 	splitSamples = samples.split(",")
# 	sampleNum = len(splitSamples)
# 	
# 	geneFeatures = list(gene[4:len(gene)-3]) #everything minus the samples, and the total score
# 	#geneFeatures.append(sampleNum)
# 	featureMatrix.append(geneFeatures)
# 	
# featureMatrix = np.array(featureMatrix, dtype="float")		
# print featureMatrix
# 
# #Get the labels based on the differential expression
# #degs = np.loadtxt(degFile, dtype="object")
# 
# #Get the cosmic genes
# cosmicGenes = []
# with open(degFile, 'rb') as f:
# 	lineCount = 0
# 	
# 	for line in f:
# 		line = line.strip()
# 		if lineCount == 0:
# 			lineCount += 1
# 			continue
# 		
# 		splitLine = line.split("\t")
# 		
# 		geneName = splitLine[0]
# 		cosmicGenes.append(geneName)
# print cosmicGenes
# # 
# # geneLabels = []
# # for gene in geneScores:
# # 	if gene[0] in degs[:,0]:
# # 		geneLabels.append(1)
# # 	else:
# # 		geneLabels.append(0)
# 
# #For COSMIC
# geneLabels = []
# for gene in geneScores:
# 	if gene[0] in cosmicGenes:
# 		geneLabels.append(1)
# 	else:
# 		geneLabels.append(0)
# 
# 		
# #print geneLabels		
# 
# 
# #Try a random forest
# X_train, X_test, y_train, y_test = train_test_split(featureMatrix, geneLabels, test_size=.4, random_state=42) #60%/40%
# print X_train.shape
# print X_test.shape
# 
# print np.where(np.array(y_train) == 1)[0].shape
# print np.where(np.array(y_test) == 1)[0].shape
# 
# 
# 
# clf = RandomForestClassifier()
# 
# clf.fit(X_train, y_train)
# importances = clf.feature_importances_
# std = np.std([tree.feature_importances_ for tree in clf.estimators_],
#              axis=0)
# indices = np.argsort(importances)[::-1]
# 
# # Print the feature ranking
# import matplotlib.pyplot as plt
# print("Feature ranking:")
# 
# for f in range(X_train.shape[1]):
#     print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
# 
# # Plot the feature importances of the forest
# plt.figure()
# plt.title("Feature importances")
# plt.bar(range(X_train.shape[1]), importances[indices],
#        color="r", yerr=std[indices], align="center")
# plt.xticks(range(X_train.shape[1]), indices)
# plt.xlim([-1, X_train.shape[1]])
# plt.show()
# 
# score = clf.score(X_test, y_test)
# print("Classification score: ", score)
# 
# #With cross validation 10 fold
# clf = RandomForestClassifier()
# folds = 10 #change this value to use a different number of folds
# scores = cross_val_score(clf, featureMatrix, geneLabels, cv=folds)
# print('classification score per fold: ', scores) #show accuracy of each fold
# print('total score: ', sum(scores) / len(scores))
# 
# import sklearn.metrics as metrics
# clf.fit(X_train, y_train)
# probs = clf.predict_proba(X_test)
# preds = probs[:,1]
# fpr, tpr, threshold = metrics.roc_curve(y_test, preds)
# roc_auc = metrics.auc(fpr, tpr)
# 
# import matplotlib.pyplot as plt
# plt.title('Receiver Operating Characteristic')
# plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
# plt.legend(loc = 'lower right')
# plt.plot([0, 1], [0, 1],'r--')
# plt.xlim([0, 1])
# plt.ylim([0, 1])
# plt.ylabel('True Positive Rate')
# plt.xlabel('False Positive Rate')
# plt.show()
# 
# exit()

#Test the gene patient pair ml setup

scoreFile = sys.argv[1]
degFile = sys.argv[2]

geneScores = np.loadtxt(scoreFile, dtype="object")
labels = []
pValues = dict()
with open(degFile, 'r') as degF:
	
	for line in degF:
		line = line.strip()
		splitLine = line.split("\t")
		
		if float(splitLine[1]) < 0.05:
			
			pValues[splitLine[0]] = 1
		else:
			pValues[splitLine[0]] = 0
print len(pValues)
print geneScores.shape
count = 0
for gene in geneScores:
	if count % 10000 == 0:
		print "count: ", count
	if gene[0] in pValues:
		labels.append(pValues[gene[0]])
	else:
		labels.append(0)
	count += 1

print labels

excludedSampleData = np.loadtxt(sys.argv[3], dtype="object")

geneIndDict = dict()

#Exclude samples that are not in the excluded sample list, sum across genes
geneInd = -1
prevGene = ""
labels = []
for genePair in geneScores:
	geneName = genePair[0].split("_")[0]
	sampleName = genePair[0].split("_")[1]
	
	if geneName != prevGene:
		geneInd += 1
		prevGene = geneName
	geneIndDict[geneName] = geneInd
	
trainingFeatureMatrix = np.empty([len(geneIndDict), geneScores.shape[1]-1], dtype="object")
trainingFeatureMatrix[:,0:30] = 0
testFeatureMatrix = np.empty([len(geneIndDict), geneScores.shape[1]-1], dtype="object")
testFeatureMatrix[:,0:30] = 0
labels = np.zeros(len(geneIndDict))
for genePair in geneScores:
	geneName = genePair[0].split("_")[0]
	sampleName = genePair[0].split("_")[1]
	
	geneInd = geneIndDict[geneName]
	
	if genePair[0] in pValues:
		labels[geneInd] = pValues[genePair[0]]
		#print pValues[geneName]
	else:
		labels[geneInd] = 0
	
	if len(excludedSampleData[excludedSampleData[:,0] == geneName,1]) < 1:
		#testFeatureMatrix[geneInd,0] = geneName
		#trainingFeatureMatrix[geneInd,0] = geneName
		continue #no samples, so the feature scores are all 0.
	excludedSamples = excludedSampleData[excludedSampleData[:,0] == geneName,1][0]
	
	excludedSampleList = excludedSamples.split(",")
	
	if sampleName in excludedSampleList:
		testFeatureMatrix[geneInd,0:30] += np.array(genePair[1:30],dtype=float)
		#testFeatureMatrix[geneInd,0] = geneName
	else:
		trainingFeatureMatrix[geneInd,0:30] += np.array(genePair[1:30],dtype=float)
		#trainingFeatureMatrix[geneInd,0] = geneName
	
#print trainingFeatureMatrix
#print testFeatureMatrix
#print labels

print len(np.where(np.array(labels) == 1)[0])
print len(np.where(np.array(labels) == 0)[0])
exit()
clf = RandomForestClassifier()

clf.fit(trainingFeatureMatrix, list(labels))
score = clf.score(testFeatureMatrix, list(labels))
print("Classification score: ", score)
import sklearn.metrics as metrics
probs = clf.predict_proba(testFeatureMatrix)
preds = probs[:,1]
fpr, tpr, threshold = metrics.roc_curve(list(labels), preds)
roc_auc = metrics.auc(fpr, tpr)

import matplotlib.pyplot as plt
plt.title('Receiver Operating Characteristic')
plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
plt.legend(loc = 'lower right')
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()

# featureMatrix = geneScores[:,1:geneScores.shape[1]-1]
# 
# 
# # print geneScores.shape
# # print labels.shape
# # print list(labels)		
# 
# clf = RandomForestClassifier()
# folds = 10 #change this value to use a different number of folds
# scores = cross_val_score(clf, featureMatrix,labels, cv=folds)
# print('classification score per fold: ', scores) #show accuracy of each fold
# print('total score: ', sum(scores) / len(scores))		
# 
# import sklearn.metrics as metrics
# X_train, X_test, y_train, y_test = train_test_split(featureMatrix, labels, test_size=.4, random_state=42) #60%/40%
# print X_train.shape
# print X_test.shape
# 
# print np.where(np.array(y_train) == 1)[0].shape
# print np.where(np.array(y_test) == 1)[0].shape
# 
# clf = RandomForestClassifier()
# 
# clf.fit(X_train, y_train)
# score = clf.score(X_test, y_test)
# print("Classification score: ", score)
# 
# 
# probs = clf.predict_proba(X_test)
# preds = probs[:,1]
# fpr, tpr, threshold = metrics.roc_curve(y_test, preds)
# roc_auc = metrics.auc(fpr, tpr)
# 
# import matplotlib.pyplot as plt
# plt.title('Receiver Operating Characteristic')
# plt.plot(fpr, tpr, 'b', label = 'AUC = %0.2f' % roc_auc)
# plt.legend(loc = 'lower right')
# plt.plot([0, 1], [0, 1],'r--')
# plt.xlim([0, 1])
# plt.ylim([0, 1])
# plt.ylabel('True Positive Rate')
# plt.xlabel('False Positive Rate')
# plt.show()



