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

#The truth data, DEG pairs
degData = np.loadtxt(sys.argv[1], dtype='object')
pathwayAnnotation = np.loadtxt(sys.argv[1] + '_pathwayAnnotation.txt', dtype='object')

#Other comparisons: non-DEG pairs, coding SV-gene pairs, germline pairs, shuffled pairs.
nonDegData = np.loadtxt(sys.argv[2], dtype='object')
codingData = np.loadtxt(sys.argv[3], dtype='object')
germlineData = np.loadtxt(sys.argv[4], dtype='object')
shuffledData = np.loadtxt(sys.argv[5], dtype='object')

#Filter by SV type if required.  
svType = ''
typeLabel = 'all SV types'
nonDEGs = []

for degPair in nonDegData:
	
	sv = degPair[0].split("_")
	if svType != '':
		typeMatch = re.search(svType, sv[8], re.IGNORECASE)
		if typeMatch is None:
			continue
	nonDEGs.append(degPair)

codingPairs = []
for pair in codingData:
	
	sv = pair[0].split("_")
	if svType != '':
		typeMatch = re.match(svType, sv[8], re.IGNORECASE)
		if typeMatch is None:
			continue
	codingPairs.append(pair)

germlinePairs = []
for pair in germlineData:
	
	sv = pair[0].split("_")

	if svType != '':
		typeMatch = re.match(svType, sv[8], re.IGNORECASE)
		if typeMatch is None:
			continue
	germlinePairs.append(pair)

shuffledPairs = []
for pair in shuffledData:
	
	sv = pair[0].split("_")
	if svType != '':
		typeMatch = re.match(svType, sv[8], re.IGNORECASE)
		if typeMatch is None:
			continue
	shuffledPairs.append(pair)

degs = []
for degPair in degData:
	
	sv = degPair[0].split("_")

	if svType != '':
		typeMatch = re.match(svType, sv[8], re.IGNORECASE)
		if typeMatch is None:
			continue

	degs.append(degPair)

degs = np.array(degs)
nonDEGs = np.array(nonDEGs)
codingPairs = np.array(codingPairs)
germlinePairs = np.array(germlinePairs)
shuffledPairs = np.array(shuffledPairs)

print(degs)
print(nonDEGs)
print(codingPairs)
print(germlinePairs)
print(shuffledPairs)

print(degs.shape)
print(nonDEGs.shape)
print(codingPairs.shape)
print(germlinePairs.shape)
print(shuffledPairs.shape)

pathwayDEGs = []
unAnnotatedDEGs = []
for degPair in degs:

	features = degPair[1:]
	
	if degPair[0] in pathwayAnnotation[:,0]:
		pathwayDEGs.append(features)
	else:
		unAnnotatedDEGs.append(features)
	
#pathwayDEGs = np.array(pathwayDEGs)
unAnnotatedDEGs = np.array(unAnnotatedDEGs)

print(unAnnotatedDEGs.shape[0])
print(nonDEGs.shape[0])


#To make the bar plots, add SNV/pathway annotation as stacked bars.
#First add the SNVs, then the pathwyas, a set with both pathways and SNVs, and the remaining are the pairs that are in neither. 

def getLossData(features, totalFeatures): #provide the features for the set of pairs that we are interested in. E.g. for SNVs only, or for pathways only
	
	if len(features) < 1:
		return np.array([])
		
	leftLosses = []
	for i in range(0, 23): #+1 because this goes until 22
		leftLosses.append(features[:,i].astype(float))
	
	#for i in range(46, 51):
	#	leftLosses.append(features[:,i].astype(float))
	
	leftLosses = np.array(leftLosses)
	
	lossData = []
	for loss in leftLosses:
		lossData.append(np.sum(loss))
	
	lossData = np.array(lossData)
	lossData = (lossData / float(totalFeatures.shape[0]))
	
	return np.array(lossData)

#pathwayLosses = getLossData(pathwayDEGs, degData)
unAnnotatedLosses = getLossData(unAnnotatedDEGs, unAnnotatedDEGs)

#non-DEGs
nonDegLosses = getLossData(nonDEGs[:,1:], nonDEGs) ####mind here when looking at SV-types individually, the normalization is then off
#normalize
unAnnotatedLossesNormND = unAnnotatedLosses / nonDegLosses
unAnnotatedLossesNormGL = []
unAnnotatedLossesNormC = []
unAnnotatedLossesNormS = []

if len(germlinePairs) > 0:
	germlineLosses = getLossData(germlinePairs[:,1:], germlinePairs)
	unAnnotatedLossesNormGL = unAnnotatedLosses / germlineLosses

if len(codingPairs) > 0:	
	codingLosses = getLossData(codingPairs[:,1:], codingPairs)
	unAnnotatedLossesNormC = unAnnotatedLosses / codingLosses
	
if len(shuffledPairs) > 0:
	shuffledLosses = getLossData(shuffledPairs[:,1:], shuffledPairs)
	unAnnotatedLossesNormS = unAnnotatedLosses / shuffledLosses





#Stack the bars for the DEG data

#Everything in the same plot

# width = 0.25
# 
# if len(unAnnotatedLossesNormND) > 0:
# 	plt.barh(np.arange(len(unAnnotatedLossesNormND)), unAnnotatedLossesNormND, width, label='DEG pairs / non-DEG pairs', color='red')
# 
# if len(unAnnotatedLossesNormGL) > 0:
# 	plt.barh(np.arange(len(unAnnotatedLossesNormGL))+width, unAnnotatedLossesNormGL, width, label='DEG pairs / germline pairs', color='blue')
# 	
# if len(unAnnotatedLossesNormS) > 0:
# 	plt.barh(np.arange(len(unAnnotatedLossesNormS))+width*2, unAnnotatedLossesNormS, width, label='DEG pairs / shuffled pairs', color='green')
# 	
# if len(unAnnotatedLossesNormC) > 0:
# 	plt.barh(np.arange(len(unAnnotatedLossesNormC))+width*3, unAnnotatedLossesNormC, width, label='DEG pairs / coding pairs', color='purple')
# 	
# plt.yticks(([p + 2.5 * width for p in np.arange(len(unAnnotatedLossesNormND))]), ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
# 												  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'CTCF', 'CTCF+Enhancer',
# 												  'CTCF+Promoter', 'chromHMM Enhancer', 'heterochromatin', 'poised promoter',
# 												  'chromHMM Promoter', 'repeat', 'repressed', 'transcribed'])
# #plt.xlim([-2,2])
# plt.legend(loc='best')
# plt.tight_layout()
# plt.show()
# #plt.savefig('Output/degNonDeg_losses.svg')
# plt.clf()


#Gains

def getGainData(features, totalFeatures): #provide the features for the set of pairs that we are interested in. E.g. for SNVs only, or for pathways only
	
	if len(features) < 1:
		return np.array([])
			
		
	leftGains = []
	for i in range(23, 46): #+1 because this goes until 45
		leftGains.append(features[:,i].astype(float))
	
	#for i in range(46, 51):
	#	leftGains.append(features[:,i].astype(float))
	
	leftGains = np.array(leftGains)
	
	gainData = []
	for loss in leftGains:
		gainData.append(np.sum(loss))
	
	gainData = np.array(gainData)
	gainData = (gainData / float(totalFeatures.shape[0]))

	return np.array(gainData)

#pathwayGains = getGainData(pathwayDEGs, degData)
unAnnotatedGains = getGainData(unAnnotatedDEGs, unAnnotatedDEGs)

nonDegGains = getGainData(nonDEGs[:,1:], nonDEGs)
#normalize
unAnnotatedGainsNormND = unAnnotatedGains / nonDegGains #np.log
unAnnotatedGainsNormC = []
unAnnotatedGainsNormS = []
unAnnotatedGainsNormGL = []

if len(codingPairs) > 0:	
	codingGains = getGainData(codingPairs[:,1:], codingPairs)
	unAnnotatedGainsNormC = unAnnotatedGains / codingGains #np.log
if len(shuffledPairs) > 0:	
	shuffledGains = getGainData(shuffledPairs[:,1:], shuffledPairs)
	unAnnotatedGainsNormS = unAnnotatedGains / shuffledGains #np.log
if len(germlinePairs) > 0:	
	germlineGains = getGainData(germlinePairs[:,1:], germlinePairs)
	unAnnotatedGainsNormGL = unAnnotatedGains / germlineGains #np.log

#Stack the bars for the DEG data

# width = 0.25
# 
# if len(unAnnotatedGainsNormND) > 0:
# 	plt.barh(np.arange(len(unAnnotatedGainsNormND)), unAnnotatedGainsNormND, width, label='DEG pairs / non-DEG pairs', color='red')
# 
# if len(unAnnotatedGainsNormGL) > 0:
# 	plt.barh(np.arange(len(unAnnotatedGainsNormGL))+width, unAnnotatedGainsNormGL, width, label='DEG pairs / germline pairs', color='blue')
# 	
# if len(unAnnotatedGainsNormS) > 0:
# 	plt.barh(np.arange(len(unAnnotatedGainsNormS))+width*2, unAnnotatedGainsNormS, width, label='DEG pairs / shuffled pairs', color='green')
# 	
# if len(unAnnotatedGainsNormC) > 0:
# 	plt.barh(np.arange(len(unAnnotatedGainsNormC))+width*3, unAnnotatedGainsNormC, width, label='DEG pairs / coding pairs', color='purple')
# 
# plt.yticks(([p + 2.5 * width for p in np.arange(len(unAnnotatedGainsNormND))]), ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
# 												  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'CTCF', 'CTCF+Enhancer',
# 												  'CTCF+Promoter', 'chromHMM Enhancer', 'heterochromatin', 'poised promoter',
# 												  'chromHMM Promoter', 'repeat', 'repressed', 'transcribed'])
# #plt.xlim([-2,2])
# plt.legend(loc='best')
# plt.tight_layout()
# plt.show()
# #plt.savefig('Output/degNonDeg_losses.svg')
# plt.clf()

### 4 separate plots, but gains and losses in one

#DEG vs non-DEG
fig = plt.figure()

def plotGainsLossesSamePlot(losses,gains, label, typeLabel, xlabel,  figInd):
	
	#plt.subplot(2, 2, figInd)
	
	width = 0.25
	
	losses = np.log(losses)
	gains = np.log(gains)
	
	print(losses)
	
	if len(losses) > 0:
		plt.barh(np.arange(len(losses)), losses, width, label='Losses', color='blue')
	
	if len(gains) > 0:
		plt.barh(np.arange(len(gains))+width, gains, width, label='Gains', color='red')
	
	plt.yticks(([p + 1.0 * width for p in np.arange(len(losses))]), ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
													  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'CTCF', 'CTCF+enhancer',
													  'CTCF+promoter', 'chromHMM enhancer', 'heterochromatin', 'poised promoter',
													  'chromHMM promoter', 'repeat', 'repressed', 'transcribed'])
	
	plt.xlim([-2,3.5])
	plt.xlabel(xlabel)
	plt.title(label + ': ' + typeLabel)
	plt.legend(loc='best')
	plt.tight_layout()
	plt.savefig('gains_losses_' + label + '_' + svType + '.svg')
	plt.clf()
	
	

plotGainsLossesSamePlot(unAnnotatedLossesNormND, unAnnotatedGainsNormND, 'DEG pairs vs. non-DEG pairs', typeLabel, 'log(% of DEG pairs / % of non-DEG pairs)', 1)
plotGainsLossesSamePlot(unAnnotatedLossesNormGL, unAnnotatedGainsNormGL, 'DEG pairs vs. germline pairs', typeLabel, 'log(% of DEG pairs / % of germline pairs)', 2)
plotGainsLossesSamePlot(unAnnotatedLossesNormC, unAnnotatedGainsNormC, 'DEG pairs vs. coding pairs', typeLabel, 'log(% of DEG pairs / % of coding pairs)', 3)
plotGainsLossesSamePlot(unAnnotatedLossesNormS, unAnnotatedGainsNormS, 'DEG pairs vs. shuffled pairs', typeLabel, 'log(% of DEG pairs / % of shuffled pairs)', 4)
#plt.tight_layout()
#plt.savefig('gains_losses_' + svType + '.svg')
#plt.show()

exit()
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
