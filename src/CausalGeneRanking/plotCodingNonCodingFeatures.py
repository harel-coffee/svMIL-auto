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

svType = ''
shuffle = False

#to shuffle, make sure that we keep the same set size. 
allDegPairs = np.concatenate((degData[:,0], nonDegData[:,0]))
allDegData = np.concatenate((degData, nonDegData))

shuffledDegPairs = np.random.choice(allDegPairs, degData.shape[0], replace=False)
shuffledNonDegData = []
shuffledDegData = []
for pair in allDegData:
	if pair[0] not in shuffledDegPairs:
		shuffledNonDegData.append(pair)
	else:
		shuffledDegData.append(pair)
		
shuffledNonDegData = np.array(shuffledNonDegData, dtype='object')
shuffledDegData = np.array(shuffledDegData, dtype='object')

nonDEGs = []

for degPair in nonDegData:
	
	sv = degPair[0].split("_")
	if svType != '':
		if sv[8] != svType:
			continue
	nonDEGs.append(degPair)

degs = []
for degPair in degData:
	
	sv = degPair[0].split("_")
	if svType != '':
		if sv[8] != svType:
			continue

	degs.append(degPair)
	
degs = np.array(degs)
nonDEGs = np.array(nonDEGs)

snvDEGs = []
pathwayDEGs = []
snvAndPathwayDEGs = []
unAnnotatedDEGs = []
for degPair in degs:

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

print(snvDEGs.shape[0] + pathwayDEGs.shape[0] + snvAndPathwayDEGs.shape[0] + unAnnotatedDEGs.shape[0])
print(nonDEGs.shape[0])

#get for shuffled
unAnnotatedDEGsPairsShuffled = np.random.choice(shuffledDegData[:,0], unAnnotatedDEGs.shape[0], replace=False)
unAnnotatedDEGsShuffled = []
for pair in unAnnotatedDEGsPairsShuffled:
	unAnnotatedDEGsShuffled.append(shuffledDegData[shuffledDegData[:,0] == pair][0][1:])

unAnnotatedDEGsShuffled = np.array(unAnnotatedDEGsShuffled, dtype='object')


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

#non-DEGs
nonDegLosses = getLossData(nonDEGs[:,1:], nonDegData)

#normalize
#1. get random losses of the same size as the unannotated DEGs
#degLossesShuffled = getLossData(unAnnotatedDEGsShuffled, degData)
degLossesShuffled = getLossData(unAnnotatedDEGsShuffled, shuffledNonDegData)
unAnnotatedLossesNormShuffled = np.log(unAnnotatedLosses / degLossesShuffled)

unAnnotatedLossesNorm = np.log(unAnnotatedLosses / nonDegLosses)
snvAndPathwayLossesNorm = snvAndPathwayLosses / nonDegLosses
pathwayLossesNorm = pathwayLosses / nonDegLosses
snvLossesNorm = snvLosses / nonDegLosses

print(unAnnotatedLosses)
print(degLossesShuffled)
print(unAnnotatedLosses / degLossesShuffled)

#Stack the bars for the DEG data
width = 0.35
# if len(unAnnotatedLossesNorm) > 0:
# 	plt.barh(np.arange(len(unAnnotatedLosses)), unAnnotatedLossesNorm, width, label='DEG pairs', color='red')

if shuffle == True:
	if len(unAnnotatedLossesNorm) > 0:
		plt.barh(np.arange(len(unAnnotatedLosses)), unAnnotatedLossesNorm, width, label='DEG pairs / non-DEG pairs', color='red')
	if len(unAnnotatedLossesNormShuffled) > 0:
		plt.barh(np.arange(len(unAnnotatedLossesNormShuffled))+ width, unAnnotatedLossesNormShuffled, width, label='DEG pairs / shuffled DEG labels', color='blue')
else:
	if len(unAnnotatedLossesNorm) > 0:
		plt.barh(np.arange(len(unAnnotatedLosses)), unAnnotatedLossesNorm, width, label='DEG pairs / non-DEG pairs', color='red')
# if len(snvLossesNorm) > 0:
# 	plt.barh(np.arange(len(snvLosses)), snvLossesNorm, width, unAnnotatedLossesNorm, label='DEG pairs + SNVs', color='purple')
# if len(pathwayLossesNorm) > 0:
# 	plt.barh(np.arange(len(pathwayLosses)), pathwayLossesNorm, width, (unAnnotatedLossesNorm+snvLossesNorm), label='DEG pairs + pathways', color='cyan')
# if len(snvAndPathwayLossesNorm) > 0:
# 	plt.barh(np.arange(len(snvAndPathwayLosses)), snvAndPathwayLossesNorm, width, (unAnnotatedLossesNorm+snvLossesNorm+pathwayLossesNorm), label='DEG pairs SNVs + pathways', color='green')

# 
# plt.barh(np.arange(len(nonDegLosses)) + width, nonDegLosses, width, label='non-DEG pairs', color='blue')
plt.yticks(np.arange(len(nonDegLosses) + width / 2), ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
												  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'CTCF', 'CTCF+Enhancer',
												  'CTCF+Promoter', 'chromHMM Enhancer', 'heterochromatin', 'poised promoter',
												  'chromHMM Promoter', 'repeat', 'repressed', 'transcribed',
												  'Deletions', 'Duplications', 'Inversions', 'Translocations',
												  'COSMIC'])
plt.xlim([-2,1])
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

#non-DEGs
nonDegGains = getGainData(nonDEGs[:,1:], nonDegData)

#normalize

#degGainsShuffled = getGainData(unAnnotatedDEGsShuffled, degData)
degGainsShuffled = getGainData(unAnnotatedDEGsShuffled, shuffledNonDegData)
unAnnotatedGainsNormShuffled = np.log(unAnnotatedGains / degGainsShuffled)

unAnnotatedGainsNorm = np.log(unAnnotatedGains / nonDegGains)
snvAndPathwayGainsNorm = snvAndPathwayGains / nonDegGains
pathwayGainsNorm = pathwayGains / nonDegGains
snvGainsNorm = snvGains / nonDegGains


#Stack the bars for the DEG data
width = 0.35
if shuffle == True:
	if len(unAnnotatedGains) > 0:
		plt.barh(np.arange(len(unAnnotatedGainsNorm)), unAnnotatedGainsNorm, width, label='DEG pairs / non-DEG pairs', color='red')
	if len(unAnnotatedGains) > 0:
		plt.barh(np.arange(len(unAnnotatedGainsNormShuffled)) + width, unAnnotatedGainsNormShuffled, width, label='DEG pairs / shuffled DEG labels', color='blue')
else:
	if len(unAnnotatedGains) > 0:
		plt.barh(np.arange(len(unAnnotatedGainsNorm)), unAnnotatedGainsNorm, width, label='DEG pairs / non-DEG pairs', color='red')
# if len(snvGains) > 0:
# 	plt.barh(np.arange(len(snvGainsNorm)), snvGainsNorm, width, unAnnotatedGainsNorm, label='DEG pairs + SNVs', color='purple')
# if len(pathwayGains) > 0:
# 	plt.barh(np.arange(len(pathwayGainsNorm)), pathwayGainsNorm, width, unAnnotatedGainsNorm+snvGainsNorm, label='DEG pairs + pathways', color='cyan')
# if len(snvAndPathwayGains) > 0:
# 	plt.barh(np.arange(len(snvAndPathwayGainsNorm)), snvAndPathwayGainsNorm, width, unAnnotatedGainsNorm+snvGainsNorm+pathwayGainsNorm, label='DEG pairs SNVs + pathways', color='green')

#plt.barh(np.arange(len(nonDegGains)) + width, nonDegGains, width, label='non-DEG pairs', color='blue')
plt.yticks(np.arange(len(nonDegGains) + width / 2), ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
												  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'CTCF', 'CTCF+Enhancer',
												  'CTCF+Promoter', 'chromHMM Enhancer', 'heterochromatin', 'poised promoter',
												  'chromHMM Promoter', 'repeat', 'repressed', 'transcribed',
												  'Deletions', 'Duplications', 'Inversions', 'Translocations',
												  'COSMIC'])
plt.xlim([-2,1])
plt.legend(loc='best')
plt.tight_layout()
plt.show()
#plt.savefig('Output/degNonDeg_losses.svg')
plt.clf()

allFeatures = np.concatenate((nonDEGs[:,1:], unAnnotatedDEGs), axis=0)

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
	
	if i < nonDEGs.shape[0]:
		colorLabels.append('b')
		labels.append(0)
	elif i >= nonDEGs.shape[0] and i < (nonDEGs.shape[0] + unAnnotatedDEGs.shape[0]):
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
from sklearn import model_selection
from sklearn.metrics import average_precision_score
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score

from random import shuffle

def cvClassification(similarityMatrix, bagLabels, clf):
	
	
	scoring = {'accuracy' : make_scorer(accuracy_score), 
			   'precision' : make_scorer(precision_score),
			   'recall' : make_scorer(recall_score), 
			   'f1_score' : make_scorer(f1_score),
			   'average_precision' : make_scorer(average_precision_score)}
	
	kfold = model_selection.StratifiedKFold(n_splits=10, random_state=42)
	
	results = model_selection.cross_validate(estimator=clf,
											  X=similarityMatrix,
											  y=bagLabels,
											  cv=kfold,
											  scoring=scoring)

	print('accuracy: ', np.mean(results['test_accuracy']), np.std(results['test_accuracy']))
	print('precision: ', np.mean(results['test_precision']), np.std(results['test_precision']))
	print('recall: ', np.mean(results['test_recall']), np.std(results['test_recall']))
	print('F1 score: ', np.mean(results['test_f1_score']), np.std(results['test_f1_score']))
	print('AP: ', np.mean(results['test_average_precision']), np.std(results['test_average_precision']))
	

print("Random forest")
from sklearn.ensemble import RandomForestClassifier
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
cvClassification(allFeatures, labels, rfClassifier)

print("Shuffled:")
shuffle(labels)
cvClassification(allFeatures, labels, rfClassifier)
