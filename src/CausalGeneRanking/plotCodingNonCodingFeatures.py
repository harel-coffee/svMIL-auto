"""
	Get the feature values of the coding & non-coding SVs, and see if these can be separated using PCA. 

"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
from six.moves import range
import re
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_score
from scipy import interp
import random

degData = np.loadtxt(sys.argv[1], dtype='object')
nonDegData = np.loadtxt(sys.argv[2], dtype='object')
snvAnnotation = np.loadtxt(sys.argv[3], dtype='object')
pathwayAnnotation = np.loadtxt(sys.argv[4], dtype='object')

#1. Split the DEG dataset into 4: 1 for SNVs, 1 for pathways, 1 for both, and 1 for neither. 
#if the sets are empty because we are not using it or there are no SNVs for example, this should not affect things such as shuffling.

snvDEGs = []
pathwayDEGs = []
snvAndPathwayDEGs = []
unAnnotatedDEGs = []
for degPair in degData:
	
	#Here include per-SV type if necessary
	features = degPair[1:]
	
	if degPair[0] in snvAnnotation[:,0] and degPair[0] in pathwayAnnotation[:,0]:
		snvAndPathwayDEGs.append(features)
	elif degPair[0] in snvAnnotation[:,0]:
		snvDEGs.append(features)
	elif degPair[0] in pathwayAnnotation[:,0]:
		pathwayDEGs.append(features)
	else:
		unAnnotatedDEGs.append(features)
	
snvDEGs = np.array(snvDEGs)
pathwayDEGs = np.array(pathwayDEGs)
snvAndPathwayDEGs = np.array(snvAndPathwayDEGs)
unAnnotatedDEGs = np.array(unAnnotatedDEGs)

# 
# 
# degData = np.array(filteredDegData, dtype='object')
# 
# #Split into features for deg pairs and non-deg pairs.
# #allow option to split this per SV type
# 
# svTypeDeg = ""
# svTypeNonDeg = ''
# 
# #Shuffle labels, randomly assign each entry to DEG or non-DEG. Will remain unbalanced!
# shuffle = False
# if svTypeDeg != "":
# 	
# 	degFeatures = []
# 	nonDegFeatures = []
# 	
# 	if shuffle == False:
# 		for sv in degData:
# 			
# 			splitSV = sv[0].split("_")
# 			features = []
# 			if len(splitSV) == 9: #for del, dup, inv
# 				if re.search(splitSV[8], svTypeDeg, re.IGNORECASE):
# 					features = sv[1:]
# 			else:
# 				if re.search("_".join([splitSV[8], splitSV[9]]), svTypeDeg, re.IGNORECASE):
# 					features = sv[1:]
# 			
# 			if len(features) < 1:
# 				continue
# 			
# 			if shuffle == False:
# 				
# 				degFeatures.append(features)
# 			if shuffle == True:
# 				flip = random.randint(0,1)
# 				if flip == 0:
# 					degFeatures.append(features)
# 				else:
# 					nonDegFeatures.append(features)
# 		
# 		for sv in nonDegData:
# 			
# 			splitSV = sv[0].split("_")
# 			features = []
# 			if len(splitSV) == 9: #for del, dup, inv
# 				if re.search(splitSV[8], svTypeNonDeg, re.IGNORECASE):
# 					features = sv[1:]
# 			else:
# 				if re.search("_".join([splitSV[8], splitSV[9]]), svTypeNonDeg, re.IGNORECASE):
# 					features = sv[1:]
# 			
# 			if len(features) < 1:
# 				continue
# 			
# 			if shuffle == False:
# 				nonDegFeatures.append(features)
# 			if shuffle == True:
# 				flip = random.randint(0,1)
# 				if flip == 0:
# 					degFeatures.append(features)
# 				else:
# 					nonDegFeatures.append(features)		
# 	
# 	degFeatures = np.array(degFeatures)
# 	nonDegFeatures = np.array(nonDegFeatures)
# else:
# 	
# 	
# 	if shuffle == False:
# 		degFeatures = degData[:,1:]
# 		nonDegFeatures = nonDegData[:,1:]
# 	else:
# 		degFeatures = []
# 		nonDegFeatures = []
# 		for sv in degData:
# 			flip = random.randint(0,1)
# 			if flip == 0:
# 				degFeatures.append(sv[1:])
# 			else:
# 				nonDegFeatures.append(sv[1:])
# 				
# 		for sv in nonDegData:
# 			flip = random.randint(0,1)
# 			if flip == 0:
# 				degFeatures.append(sv[1:])
# 			else:
# 				nonDegFeatures.append(sv[1:])
# 				
# 		degFeatures = np.array(degFeatures)
# 		nonDegFeatures = np.array(nonDegFeatures)
# 
# print(degFeatures.shape)
# print(nonDegFeatures.shape)

#To make the bar plots, add SNV/pathway annotation as stacked bars.
#First add the SNVs, then the pathwyas, a set with both pathways and SNVs, and the remaining are the pairs that are in neither. 

def getLossData(features, totalFeatures): #provide the features for the set of pairs that we are interested in. E.g. for SNVs only, or for pathways only
	
	if len(features) < 1:
		return np.array([])
		
	leftLosses = []
	for i in range(0, 23): #+1 because this goes until 22
		leftLosses.append(features[:,i].astype(float))
	
	for i in range(46, 51):
		leftLosses.append(features[:,i].astype(float))
	
	leftLosses = np.array(leftLosses)
	
	lossData = []
	for loss in leftLosses:
		lossData.append(np.sum(loss))
	
	lossData = np.array(lossData)
	lossData = (lossData / float(totalFeatures.shape[0]))
	
	return np.array(lossData)

#1. Set with SNVs only
snvLosses = getLossData(snvDEGs, degData)
pathwayLosses = getLossData(pathwayDEGs, degData)
snvAndPathwayLosses = getLossData(snvAndPathwayDEGs, degData)
unAnnotatedLosses = getLossData(unAnnotatedDEGs, degData)


#Stack the bars for the DEG data
width = 0.35
plt.barh(np.arange(len(unAnnotatedLosses)), unAnnotatedLosses, width, label='DEG pairs without SNV or pathway effects', color='red')
plt.barh(np.arange(len(snvLosses)), snvLosses, width, unAnnotatedLosses, label='DEG pairs with SNVs', color='cyan')
plt.barh(np.arange(len(pathwayLosses)), pathwayLosses, width, unAnnotatedLosses+snvLosses, label='DEG pairs with pathway effects', color='orange')
if len(snvAndPathwayLosses) > 0:
	plt.barh(np.arange(len(snvAndPathwayLosses)), snvAndPathwayLosses, width, unAnnotatedLosses+snvLosses+pathwayLosses, label='DEG pairs with both SNV and pathway effects', color='green')



#non-DEGs
nonDegLosses = getLossData(nonDegData[:,1:], nonDegData)

plt.barh(np.arange(len(nonDegLosses)) + width, nonDegLosses, width, label='non-DEG pairs', color='blue')
plt.yticks(np.arange(len(nonDegLosses) + width / 2), ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
												  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'CTCF', 'CTCF+Enhancer',
												  'CTCF+Promoter', 'chromHMM Enhancer', 'heterochromatin', 'poised promoter',
												  'chromHMM Promoter', 'repeat', 'repressed', 'transcribed',
												  'Deletions', 'Duplications', 'Inversions', 'Translocations',
												  'COSMIC'])
plt.xlim([0,1])
plt.legend(loc='best')
plt.tight_layout()
plt.show()
#plt.savefig('Output/degNonDeg_losses.svg')
plt.clf()

#Gains

def getGainData(features, totalFeatures): #provide the features for the set of pairs that we are interested in. E.g. for SNVs only, or for pathways only
	
	if len(features) < 1:
		return np.array([])
			
		
	leftGains = []
	for i in range(23, 46): #+1 because this goes until 45
		leftGains.append(features[:,i].astype(float))
	
	for i in range(46, 51):
		leftGains.append(features[:,i].astype(float))
	
	leftGains = np.array(leftGains)
	
	gainData = []
	for loss in leftGains:
		gainData.append(np.sum(loss))
	
	gainData = np.array(gainData)
	gainData = (gainData / float(totalFeatures.shape[0]))

	return np.array(gainData)

#1. Set with SNVs only
snvGains = getGainData(snvDEGs, degData)
pathwayGains = getGainData(pathwayDEGs, degData)
snvAndPathwayGains = getGainData(snvAndPathwayDEGs, degData)
unAnnotatedGains = getGainData(unAnnotatedDEGs, degData)


#Stack the bars for the DEG data
width = 0.35
plt.barh(np.arange(len(unAnnotatedGains)), unAnnotatedGains, width, label='DEG pairs without SNV or pathway effects', color='red')
plt.barh(np.arange(len(snvGains)), snvGains, width, unAnnotatedGains, label='DEG pairs with SNVs', color='cyan')
plt.barh(np.arange(len(pathwayGains)), pathwayGains, width, unAnnotatedGains+snvGains, label='DEG pairs with pathway effects', color='orange')
if len(snvAndPathwayGains) > 0:
	plt.barh(np.arange(len(snvAndPathwayGains)), snvAndPathwayGains, width, unAnnotatedGains+snvGains+pathwayGains, label='DEG pairs with both SNV and pathway effects', color='green')

#non-DEGs
nonDegGains = getGainData(nonDegData[:,1:], nonDegData)

plt.barh(np.arange(len(nonDegGains)) + width, nonDegGains, width, label='non-DEG pairs', color='blue')
plt.yticks(np.arange(len(nonDegGains) + width / 2), ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
												  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'CTCF', 'CTCF+Enhancer',
												  'CTCF+Promoter', 'chromHMM Enhancer', 'heterochromatin', 'poised promoter',
												  'chromHMM Promoter', 'repeat', 'repressed', 'transcribed',
												  'Deletions', 'Duplications', 'Inversions', 'Translocations',
												  'COSMIC'])
plt.xlim([0,1])
plt.legend(loc='best')
plt.tight_layout()
plt.show()
#plt.savefig('Output/degNonDeg_losses.svg')
plt.clf()




allFeatures = np.concatenate((nonDegData[:,1:], degData[:,1:]), axis=0)

#Make a PCA plot for the left/right set and see if these are really different

from sklearn.decomposition import PCA

#Get subset of PCA

pca = PCA(n_components=2)

projected = pca.fit_transform(allFeatures)
# projectedWithOffset = projected
# 
# jitter = [0.01, -0.01]
# for row in range(0, projected.shape[0]):
# 	for col in range(0, projected.shape[1]):
# 		projectedWithOffset[row][col] += np.random.normal(-1, 1) * 0.1
# 		
# projected = projectedWithOffset

colorLabels = []
labels = []

for i in range(0, allFeatures.shape[0]):
	
	if i < nonDegData.shape[0]:
		colorLabels.append('b')
		labels.append(0)
	elif i >= nonDegData.shape[0] and i < (nonDegData.shape[0] + degData.shape[0]):
		colorLabels.append('r')
		labels.append(1)

fig,ax=plt.subplots(figsize=(7,5))
plt.scatter(projected[:, 0], projected[:, 1], edgecolors=colorLabels, facecolors='none')
plt.show()

#rasterize the PCA plot and make a density heatmap
import math

#
colorLabels = np.array(colorLabels)

#Get the minimum and maximum to determine the bounds of the plot.
xmin = np.min(projected[:,0])
xmax = np.max(projected[:,0])
ymin = np.min(projected[:,1])
ymax = np.max(projected[:,1])

#Define the box size and how many boxes we should make
print(xmin, xmax, ymin, ymax)

#round the values to get covering boxes
xmin = round(xmin)
xmax = round(xmax)
ymin = round(ymin)
ymax = round(ymax)

boxWidth = 0.2
#Take the ceil to get the maximum possible without leaving out points
xBoxNum = int(math.ceil((xmax - xmin) / boxWidth))
yBoxNum = int(math.ceil((ymax - ymin) / boxWidth))

#Placeholder for smoothed data
plotGrid = np.zeros([xBoxNum, yBoxNum])

#Loop through the data and show the data in the boxes
yBoxStart = ymin
yBoxEnd = ymin + boxWidth
xBoxStart = xmin
xBoxEnd = xmin + boxWidth
for yInd in range(0, yBoxNum):
	for xInd in range(0, xBoxNum):
		
		#Find all data points that are within the current box
		xStartMatches = projected[:,0] >= xBoxStart
		xEndMatches = projected[:,0] <= xBoxEnd
		
		xMatches = xStartMatches * xEndMatches
		
		yStartMatches = projected[:,1] >= yBoxStart
		yEndMatches = projected[:,1] <= yBoxEnd
		
		yMatches = yStartMatches * yEndMatches
		
		dataInBox = projected[xMatches * yMatches]
		boxLabels = colorLabels[xMatches * yMatches]
		
		if len(dataInBox) > 0:
			#print dataInBox
			
			posCount = len(np.where(boxLabels == 'r')[0]) + 0.01
			negCount = len(np.where(boxLabels == 'b')[0]) + 0.01
			
			#Normalize for the total count of that label
			posCount = posCount / len(np.where(colorLabels == 'r')[0])
			negCount = negCount / len(np.where(colorLabels == 'b')[0])
			
			if negCount > 0:
				plotGrid[xInd,yInd] = np.log(posCount / float(negCount))
			

		#Move the box along x
		xBoxStart += boxWidth
		xBoxEnd += boxWidth
	
	yBoxStart += boxWidth
	yBoxEnd += boxWidth
	#Reset the box on x
	xBoxStart = xmin
	xBoxEnd = xmin + boxWidth

plotGrid = np.ma.masked_where(plotGrid == 0, plotGrid)
cmap = plt.cm.seismic
cmap.set_bad(color='white')
print(plotGrid)
plt.imshow(plotGrid, cmap=cmap, interpolation='nearest')		
plt.show()

#very simple ml test
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import auc, precision_recall_curve

from random import shuffle


def performCV(featureMatrix, labels, clf):
	
	cv = StratifiedKFold(n_splits=10)
	
	#dirty
	X = featureMatrix
	y = np.array(labels)
	
	tprs = []
	aucs = []
	auprcs = []
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

		precision, recall, thresholds = precision_recall_curve(y[test], predictions)
		aucScore = auc(recall, precision)
		auprcs.append(aucScore)
		
		i += 1
	
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	mean_score = np.mean(scores)
	
	
	
	print("Score: ", mean_score)
	print("AUC: ", mean_auc)
	print("AUPRC: ", np.mean(auprcs))

print("Random forest")
from sklearn.ensemble import RandomForestClassifier
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
performCV(allFeatures, labels, rfClassifier)

print("Shuffled:")
shuffle(labels)
performCV(allFeatures, labels, rfClassifier)
