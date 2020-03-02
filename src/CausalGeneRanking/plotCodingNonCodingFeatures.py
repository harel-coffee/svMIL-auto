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
from scipy.stats import chi2_contingency
import random

#The truth data, DEG pairs
degData = np.loadtxt(sys.argv[1], dtype='object')

#pathwayAnnotation = np.loadtxt(sys.argv[1] + '_pathwayAnnotation.txt', dtype='object')

#Other comparisons: non-DEG pairs, coding SV-gene pairs, germline pairs, shuffled pairs.
nonDegData = np.loadtxt(sys.argv[2], dtype='object')
#codingData = np.loadtxt(sys.argv[3], dtype='object')
germlineData = np.loadtxt(sys.argv[3], dtype='object')
shuffledData = np.loadtxt(sys.argv[4], dtype='object')

codingData = nonDegData

#Filter by SV type if required.  
svType = 'INV'
typeLabel = 'inversions'
nonDEGs = []

for degPair in nonDegData:
	
	sv = degPair[0].split("_")
	if svType != '':
		typeMatch = re.search(svType, sv[12], re.IGNORECASE)
		if typeMatch is None:
			continue
	nonDEGs.append(degPair)

codingPairs = []
for pair in codingData:
	
	sv = pair[0].split("_")
	if svType != '':
		typeMatch = re.match(svType, sv[12], re.IGNORECASE)
		if typeMatch is None:
			continue
	codingPairs.append(pair)

germlinePairs = []
for pair in germlineData:
	
	sv = pair[0].split("_")

	if svType != '':
		typeMatch = re.match(svType, sv[12], re.IGNORECASE)
		if typeMatch is None:
			continue
	germlinePairs.append(pair)

shuffledPairs = []
for pair in shuffledData:
	
	sv = pair[0].split("_")
	if svType != '':
		typeMatch = re.match(svType, sv[12], re.IGNORECASE)
		if typeMatch is None:
			continue
	shuffledPairs.append(pair)

degs = []
for degPair in degData:
	
	sv = degPair[0].split("_")

	if svType != '':
		typeMatch = re.match(svType, sv[12], re.IGNORECASE)
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
	
	#if degPair[0] in pathwayAnnotation[:,0]:
	#	pathwayDEGs.append(features)
	#else:
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
	for i in range(0, 26): #+1 because this goes until 22
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

#test chi2

def getLossSignificance(features, shuffledFeatures):
	
	#for each loss feature, check if more/less than by random chance.
	lossSignificances = []
	for i in range(0, 26):

		degLosses = len(np.where(features[:,i] == '1.0')[0])
		degNoLosses = len(np.where(features[:,i] == '0.0')[0])
		
		shuffledLosses = len(np.where(shuffledFeatures[:,i] == '1.0')[0])
		shuffledNoLosses = len(np.where(shuffledFeatures[:,i] == '0.0')[0])
		
		
		if degLosses == 0 or degNoLosses == 0 or shuffledLosses == 0 or shuffledNoLosses == 0:
			lossSignificances.append(1)
			continue
		
		obs = np.array([[degLosses, degNoLosses], [shuffledLosses, shuffledNoLosses]])
		
		result = chi2_contingency(obs)
		p = result[1]

		lossSignificances.append(p)

	return lossSignificances

lossSignificancesTmp = [0]*26

if svType == 'INV':
	lossSignificances = getLossSignificance(unAnnotatedDEGs, nonDEGs[:,1:])
	lossSignificancesGL = getLossSignificance(unAnnotatedDEGs, germlinePairs[:,1:])
	lossSignificancesS = getLossSignificance(unAnnotatedDEGs, shuffledPairs[:,1:])
elif svType == 'ITX':
	lossSignificances = getLossSignificance(unAnnotatedDEGs, nonDEGs[:,1:])
	lossSignificancesGL = lossSignificancesTmp
	lossSignificancesS = getLossSignificance(unAnnotatedDEGs, shuffledPairs[:,1:])
else:
	lossSignificances = lossSignificancesTmp
	lossSignificancesGL = lossSignificancesTmp
	lossSignificancesS = lossSignificancesTmp


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
	for i in range(26, 51): #+1 because this goes until 45
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

def getGainSignificance(features, shuffledFeatures):
	
	#for each gain feature, check if more/less than by random chance.
	gainSignificances = []
	for i in range(26, 52):

		degLosses = len(np.where(features[:,i] == '1.0')[0])
		degNoLosses = len(np.where(features[:,i] == '0.0')[0])
		
		shuffledLosses = len(np.where(shuffledFeatures[:,i] == '1.0')[0])
		shuffledNoLosses = len(np.where(shuffledFeatures[:,i] == '0.0')[0])
		
		if degLosses == 0 or degNoLosses == 0 or shuffledLosses == 0 or shuffledNoLosses == 0:
			gainSignificances.append(1)
			continue
			
		obs = np.array([[degLosses, degNoLosses], [shuffledLosses, shuffledNoLosses]])
		
		result = chi2_contingency(obs)
		p = result[1]

		gainSignificances.append(p)

	return gainSignificances

gainSignificancesTmp = [0]*26

if svType == 'ITX':
	gainSignificances = getGainSignificance(unAnnotatedDEGs, nonDEGs[:,1:])
	gainSignificancesGL = gainSignificancesTmp
	gainSignificancesS = getGainSignificance(unAnnotatedDEGs, shuffledPairs[:,1:])
elif svType == 'DUP':
	gainSignificances = getGainSignificance(unAnnotatedDEGs, nonDEGs[:,1:])
	gainSignificancesGL = gainSignificancesTmp
	gainSignificancesS = getGainSignificance(unAnnotatedDEGs, shuffledPairs[:,1:])
else:
	
	gainSignificances = getGainSignificance(unAnnotatedDEGs, nonDEGs[:,1:])
	#gainSignificancesGL = getGainSignificance(unAnnotatedDEGs, germlinePairs[:,1:])
	#gainSignificancesS = getGainSignificance(unAnnotatedDEGs, shuffledPairs[:,1:])
exit()


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
def plotGainsLossesSamePlot(losses,gains, lossSignificances, gainSignificances, label, typeLabel, xlabel,  figInd):
	
	significances = []
	for sign in range(0, len(lossSignificances)):
	
		if lossSignificances[sign] < 0.05:
			significances.append('*')
		else:
			significances.append('')
			
		if gainSignificances[sign] < 0.05:
			significances.append('*')
		else:
			significances.append('')
	
	#plt.subplot(2, 2, figInd)
	fig, ax = plt.subplots()
	
	width = 0.25
	
	losses = np.log(losses)
	gains = np.log(gains)
	
	print(losses)
	
	if len(losses) > 0:
		plt.barh(np.arange(len(losses)), losses, width, label='Losses', color='blue')
	
	if len(gains) > 0:
		plt.barh(np.arange(len(gains))+width, gains, width, label='Gains', color='red')
	
	if len(significances) > 0:
		rects = ax.patches
		
		# For each bar: Place a label
		ind = 0
		for rect in rects:
			# Get X and Y placement of label from rect.
			x_value = rect.get_width()
			y_value = rect.get_y() + rect.get_height() / 2
		
			# Number of points between bar and label. Change to your liking.
			space = 5
			# Vertical alignment for positive values
			ha = 'left'
		
			# If value of bar is negative: Place label left of bar
			if x_value < 0:
				# Invert space to place label to the left
				space *= -1
				# Horizontally align label at right
				ha = 'right'
		
			# Use X value as label and format number with one decimal place
			
			# Create annotation
			plt.annotate(
				significances[ind],                      # Use `label` as label
				(x_value, y_value),         # Place label at end of the bar
				xytext=(space, 0),          # Horizontally shift label by `space`
				textcoords="offset points", # Interpret `xytext` as offset in points
				va='center',                # Vertically center label
				ha=ha,
				fontsize=5)                      # Horizontally align label differently for
											# positive and negative values.
			ind += 1
		
	
	plt.yticks(([p + 1.0 * width for p in np.arange(len(losses))]), ['eQTLs', 'enhancers', 'promoters', 'CpG', 'TF', 'HiC', 'h3k9me3', 'h3k4me3',
													  'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3','DNAseI', 'RNA pol II', 'CTCF', 'CTCF+enhancer',
													  'CTCF+promoter', 'chromHMM enhancer', 'heterochromatin', 'poised promoter',
													  'chromHMM promoter', 'repeat', 'repressed', 'transcribed', 'super enhancer', 'ctcf'])
	
	plt.xlim([-2,3.5])
	plt.xlabel(xlabel)
	plt.title(label + ': ' + typeLabel)
	plt.legend(loc='best')
	plt.tight_layout()
	#plt.savefig('gains_losses_' + label + '_' + svType + '.svg')
	plt.show()
	plt.clf()

	
	

plotGainsLossesSamePlot(unAnnotatedLossesNormND, unAnnotatedGainsNormND, lossSignificances, gainSignificances, 'Positive pairs vs. negative pairs', typeLabel, 'log(% of positive pairs / % of negative pairs)', 1)
plotGainsLossesSamePlot(unAnnotatedLossesNormGL, unAnnotatedGainsNormGL, lossSignificancesGL, gainSignificancesGL, 'Positive pairs vs. germline pairs', typeLabel, 'log(% of positive pairs / % of germline pairs)', 2)
#plotGainsLossesSamePlot(unAnnotatedLossesNormC, unAnnotatedGainsNormC, 'DEG pairs vs. coding pairs', typeLabel, 'log(% of DEG pairs / % of coding pairs)', 3)
plotGainsLossesSamePlot(unAnnotatedLossesNormS, unAnnotatedGainsNormS, lossSignificancesS, gainSignificancesS, 'Positive pairs vs. shuffled pairs', typeLabel, 'log(% of positive pairs / % of shuffled pairs)', 3)
#plt.tight_layout()
#plt.savefig('gains_losses_' + svType + '.svg')
#plt.show()
exit()
np.random.seed(0)
positive = degs[:,1:]

indices = np.arange(nonDEGs.shape[0])
rnd_indices = np.random.choice(indices, size=positive.shape[0])

negative = nonDEGs[rnd_indices][:,1:]

allFeatures = np.concatenate((positive, negative), axis=0).astype(float)

#normalize
normalized = np.zeros(allFeatures.shape)
for col in range(0, allFeatures.shape[1]):
	
	if np.min(allFeatures[:,col]) == np.max(allFeatures[:,col]):
		continue
	
	print(col)
	print(np.min(allFeatures[:,col]))
	print(np.max(allFeatures[:,col]))

	normalized[:,col] = (allFeatures[:,col] - np.min(allFeatures[:,col])) / (np.max(allFeatures[:,col] - np.min(allFeatures[:,col])))
	
#normalized = (allFeatures-np.min(allFeatures))/(np.max(allFeatures)-np.min(allFeatures))
allFeatures = normalized
print(allFeatures)



#remove some features
# allFeaturesFiltered = []
# for row in range(0, allFeatures.shape[0]):
# 	rowFeatures = []	
# 	for featureInd in range(0, allFeatures.shape[1]):
# 	
# 		if featureInd > 70 and featureInd < 75:
# 			continue
# 		else:
# 			rowFeatures.append(allFeatures[row][featureInd])
# 
# 	allFeaturesFiltered.append(rowFeatures)
# 	
# allFeatures = np.array(allFeaturesFiltered)	
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

colorLabels = ['r']*positive.shape[0] + ['b']*negative.shape[0]
labels = [1]*positive.shape[0] + [0]*negative.shape[0]

fig,ax=plt.subplots(figsize=(7,5))
plt.scatter(projected[:, 0], projected[:, 1], edgecolors=colorLabels, facecolors='none')
plt.show()

#rasterize the PCA plot and make a density heatmap
# import math
# 
# #
# colorLabels = np.array(colorLabels)
# 
# #Get the minimum and maximum to determine the bounds of the plot.
# xmin = np.min(projected[:,0])
# xmax = np.max(projected[:,0])
# ymin = np.min(projected[:,1])
# ymax = np.max(projected[:,1])
# 
# #Define the box size and how many boxes we should make
# print(xmin, xmax, ymin, ymax)
# 
# #round the values to get covering boxes
# xmin = round(xmin)
# xmax = round(xmax)
# ymin = round(ymin)
# ymax = round(ymax)
# 
# boxWidth = 0.2
# #Take the ceil to get the maximum possible without leaving out points
# xBoxNum = int(math.ceil((xmax - xmin) / boxWidth))
# yBoxNum = int(math.ceil((ymax - ymin) / boxWidth))
# 
# #Placeholder for smoothed data
# plotGrid = np.zeros([xBoxNum, yBoxNum])
# 
# #Loop through the data and show the data in the boxes
# yBoxStart = ymin
# yBoxEnd = ymin + boxWidth
# xBoxStart = xmin
# xBoxEnd = xmin + boxWidth
# for yInd in range(0, yBoxNum):
# 	for xInd in range(0, xBoxNum):
# 		
# 		#Find all data points that are within the current box
# 		xStartMatches = projected[:,0] >= xBoxStart
# 		xEndMatches = projected[:,0] <= xBoxEnd
# 		
# 		xMatches = xStartMatches * xEndMatches
# 		
# 		yStartMatches = projected[:,1] >= yBoxStart
# 		yEndMatches = projected[:,1] <= yBoxEnd
# 		
# 		yMatches = yStartMatches * yEndMatches
# 		
# 		dataInBox = projected[xMatches * yMatches]
# 		boxLabels = colorLabels[xMatches * yMatches]
# 		
# 		if len(dataInBox) > 0:
# 			#print dataInBox
# 			
# 			posCount = len(np.where(boxLabels == 'r')[0]) + 0.01
# 			negCount = len(np.where(boxLabels == 'b')[0]) + 0.01
# 			
# 			#Normalize for the total count of that label
# 			posCount = posCount / len(np.where(colorLabels == 'r')[0])
# 			negCount = negCount / len(np.where(colorLabels == 'b')[0])
# 			
# 			if negCount > 0:
# 				plotGrid[xInd,yInd] = np.log(posCount / float(negCount))
# 			
# 
# 		#Move the box along x
# 		xBoxStart += boxWidth
# 		xBoxEnd += boxWidth
# 	
# 	yBoxStart += boxWidth
# 	yBoxEnd += boxWidth
# 	#Reset the box on x
# 	xBoxStart = xmin
# 	xBoxEnd = xmin + boxWidth
# 
# plotGrid = np.ma.masked_where(plotGrid == 0, plotGrid)
# cmap = plt.cm.seismic
# cmap.set_bad(color='white')
# print(plotGrid)
# plt.imshow(plotGrid, cmap=cmap, interpolation='nearest')		
# plt.show()

#very simple ml test
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import auc, precision_recall_curve
from sklearn import model_selection
from sklearn.metrics import average_precision_score, plot_roc_curve
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

from random import shuffle

def cvClassification(similarityMatrix, bagLabels, clf):
	
	
	scoring = {'accuracy' : make_scorer(accuracy_score), 
			   'precision' : make_scorer(precision_score),
			   'recall' : make_scorer(recall_score), 
			   'f1_score' : make_scorer(f1_score),
			   'average_precision' : make_scorer(average_precision_score),
			   'auc' : make_scorer(roc_auc_score)}
	
	kfold = model_selection.StratifiedKFold(n_splits=10, random_state=42)
	
	results = model_selection.cross_validate(estimator=clf,
											  X=similarityMatrix,
											  y=bagLabels,
											  cv=kfold,
											  scoring=scoring)
	
	#print('accuracy: ', np.mean(results['test_accuracy']), np.std(results['test_accuracy']))
	#print('precision: ', np.mean(results['test_precision']), np.std(results['test_precision']))
	#print('recall: ', np.mean(results['test_recall']), np.std(results['test_recall']))
	print('F1 score: ', np.mean(results['test_f1_score']), np.std(results['test_f1_score']))
	print('AP: ', np.mean(results['test_average_precision']), np.std(results['test_average_precision']))
	print('AUC: ', np.mean(results['test_auc']), np.std(results['test_auc']))
	
	bagLabels = np.array(bagLabels)
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)
	
	fig, ax = plt.subplots()
	for i, (train, test) in enumerate(kfold.split(similarityMatrix, bagLabels)):
		clf.fit(similarityMatrix[train], bagLabels[train])
		viz = plot_roc_curve(clf, similarityMatrix[test], bagLabels[test],
							 name='ROC fold {}'.format(i),
							 alpha=0.3, lw=1, ax=ax)
		interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(viz.roc_auc)
	print(aucs)
	print(np.mean(aucs))
	print(np.std(aucs))

	ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
			label='Chance', alpha=.8)
	
	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	# ax.plot(mean_fpr, mean_tpr, color='b',
	# 		label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
	# 		lw=2, alpha=.8)
	ax.plot(mean_fpr, mean_tpr, color='b',
			label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (np.mean(aucs), np.std(aucs)),
			lw=2, alpha=.8)
	
	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
					label=r'$\pm$ 1 std. dev.')
	
	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
		   title="Receiver operating characteristic: deletions")
	ax.legend(loc="lower right")
	plt.show()
	exit()
	
	return np.mean(results['test_f1_score']), np.mean(results['test_average_precision'])
	
from sklearn.ensemble import RandomForestClassifier

rfClassifier = RandomForestClassifier(max_depth=100, n_estimators=2)

from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

names = ["Random Forest",]

classifiers = [
    #KNeighborsClassifier(3),
    #SVC(kernel="linear", C=0.025),
    #SVC(gamma=2, C=1),
    #GaussianProcessClassifier(1.0 * RBF(1.0)),
    #DecisionTreeClassifier(max_depth=100),
    RandomForestClassifier(max_depth=100, n_estimators=2)]
    #MLPClassifier(alpha=1, max_iter=1000),
    #AdaBoostClassifier(),
    #GaussianNB(),
    #QuadraticDiscriminantAnalysis()]

for classifierInd in range(0, len(classifiers)):
	print(names[classifierInd])
	rfClassifier = classifiers[classifierInd]
	
	featureCount = allFeatures.shape[1]-1
	#featureCount = 1
	f1s = []
	aps = []
	for featureInd in range(featureCount, allFeatures.shape[1]):

		f1, ap = cvClassification(allFeatures[:,0:featureInd], labels, rfClassifier)
		f1s.append(f1)
		aps.append(ap)

print('shuffled: ')

for classifierInd in range(0, len(classifiers)):
	print(names[classifierInd])
	rfClassifier = classifiers[classifierInd]
	f1Shuffled = []
	apShuffled = []
	for featureInd in range(featureCount, allFeatures.shape[1]):

			
		shuffle(labels)
		f1, ap = cvClassification(allFeatures[:,0:featureInd], labels, rfClassifier)
		f1Shuffled.append(f1)
		apShuffled.append(ap)
	
	
	