"""
	Use a MIL-based representation to determine if an SV overlaps a TAD boundary or not. 

"""

import random
import numpy as np
import sys
from inputParser import InputParser
import settings

## Describe each window with the elements that are in that window.
# Get the SV-gene pairs that are DEG and non-DEG
#For the 4 MB windows, collect each element that is inside those windows.
#List: TADs, enhancers, SV, genes

pairs, pairLabels = InputParser().processSVGenePairs(sys.argv[1], sys.argv[2])
enhancerData = InputParser().getEnhancersFromFile(settings.files['enhancerFile'])
tadData = InputParser().getTADsFromFile(settings.files['tadFile'])
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
svGenePairsRules = np.loadtxt(sys.argv[3], dtype='object')

#Combine the genes into one set. 
geneData = np.concatenate((causalGenes, nonCausalGenes), axis=0)

bags = []
labels = []
currentChr = None
window = 1000000 #look in a 2 MB window
pairInd = -1
posBags = []
negBags = []
for pair in pairs:
	pairInd += 1
	pairStr = "_".join([str(i) for i in pair])
	if pairStr not in svGenePairsRules[:,0]:
		continue
	#look in a 1 MB window around the SV of the pair.
	#find enhancers, genes and TADs that are inside this window

	
	sv = pair[1:len(pair)]
	
	
	instances = []
	
	#only get a new chr subset if the chr changes
	if sv[0] != currentChr:
		#Get the subset of eQTls on this chromosome
		chr1SubsetEnh = enhancerData[np.where(enhancerData[:,0] == sv[0])]
		print(chr1SubsetEnh.shape)
		currentChr = sv[0]
		
		#Get the subset of eQTls on this chromosome
		chr1SubsetTAD = tadData[np.where(tadData[:,0] == sv[0])]
		print(chr1SubsetTAD.shape)
		
		chr1SubsetGene = geneData[np.where(geneData[:,0] == sv[0])]
		print(chr1SubsetGene.shape)
	
	startMatches = sv[1] <= chr1SubsetTAD[:,2]
	endMatches = sv[5] >= chr1SubsetTAD[:,1]
			
	tadMatches = chr1SubsetTAD[startMatches * endMatches]
	
	leftTAD = tadMatches[0]
	rightTAD = tadMatches[len(tadMatches)-1]
		
	
	#The only eQTLs we consider are the ones that are starting or ending within the window, and are not within the SV.
	#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
	startMatches = (chr1SubsetEnh[:,1] <= sv[1]) * (chr1SubsetEnh[:,1] >= leftTAD[1]) * (chr1SubsetEnh[:,2] <= sv[1])
			
	#The reverse for the end of the SV.
	endMatches = (chr1SubsetEnh[:,2] >= sv[5]) * (chr1SubsetEnh[:,2] <= rightTAD[2]) * (chr1SubsetEnh[:,1] >= sv[5])
	
	matchingStart = chr1SubsetEnh[startMatches] # must be either matching on the left or right.
	matchingEnd = chr1SubsetEnh[endMatches] # must be either matching on the left or right.
	
	#For every position, determine the position relative to the SV and add it at the right place in the feature vector.
	for match in matchingStart:
		instances.append(np.array([match[1], match[2]]))
		
	
	#TADs
	
	instances.append(np.array([rightTAD[2], leftTAD[1]]))

	
	#Genes
	
	#Match for TADs that start in the start window
	startStartMatches = (chr1SubsetGene[:,1] <= sv[1]) * (chr1SubsetGene[:,1] >= leftTAD[1])
	#Match for TADs that end in the start window		
	startEndMatches = (chr1SubsetGene[:,2] <= sv[1]) * (chr1SubsetGene[:,2] >= leftTAD[1])		
			
	#The reverse for the end ofchr1SubsetGenethe SV.
	endStartMatches = (chr1SubsetGene[:,1] >= sv[5]) * (chr1SubsetGene[:,1] <= rightTAD[2])
	endEndMatches = (chr1SubsetGene[:,2] >= sv[5]) * (chr1SubsetGene[:,2] <= rightTAD[2])
	
	matchingStart = chr1SubsetGene[startStartMatches + startEndMatches] # must be either matching on the left or right.
	matchingEnd = chr1SubsetGene[endStartMatches + endEndMatches] # must be either matching on the left or right.
	
	matches = np.concatenate((matchingStart, matchingEnd))
	
	for match in matches:
		instances.append(np.array([match[1], match[2]]))


	#SV
	instances.append(np.array([sv[1], sv[5]]))
	if pairLabels[pairInd] == 0:
		negBags.append(np.array(instances))
	else:
		posBags.append(np.array(instances))
	
posBags = np.array(posBags)
negBags = np.array(negBags)
negBagsSubsampled = np.random.choice(negBags, posBags.shape[0])	

bags = np.concatenate((posBags, negBagsSubsampled))
instances = np.vstack(bags)
labels = [1]*posBags.shape[0] + [0]*negBagsSubsampled.shape[0]

from random import shuffle
shuffle(labels)

labels = np.array(labels)

print(bags)
print(instances.shape)

#Run MILES

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

#Train random forest

print(similarityMatrix)

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits=10)
np.random.seed(500)

accs = []
aucs = []
coeffs = []
predDiffs = []
for train, test in cv.split(similarityMatrix, labels):
	
	rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
	rfClassifier.fit(similarityMatrix[train], labels[train]) #Use the bag labels, not the instance labels

	predictions = rfClassifier.predict(similarityMatrix[test])
	precision, recall, thresholds = precision_recall_curve(labels[test], predictions)
	aucScore = auc(recall, precision)
	predsDiff = np.average(labels[test] == np.sign(predictions))
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


def plotData(featureMatrix, windowLabels):
		
	from sklearn.decomposition import PCA
	import matplotlib.pyplot as plt
	import math
	from tsne import bh_sne
	
	pca = PCA(n_components=2)
	projected = pca.fit_transform(featureMatrix)
	
	projectedNoise = np.zeros([projected.shape[0], projected.shape[1]])
	for col in range(0, projected.shape[1]):
		for row in range(0, projected.shape[0]):
	
			projectedNoise[row][col] = projected[row][col] + np.random.normal(0, 0.05)
	
	projected = projectedNoise
	
	colorLabels = []
	
	for label in windowLabels:
		
		if label == 1:
			colorLabels.append('r')
		else:
			colorLabels.append('b')
	
	fig,ax=plt.subplots(figsize=(7,5))
	plt.scatter(projected[:, 0], projected[:, 1], c=colorLabels)
	plt.show()
	exit()
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
	
	boxWidth = 100000
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
	plt.imshow(plotGrid, cmap=cmap, interpolation='nearest')		
	plt.show()
	
	
	vis_data = bh_sne(featureMatrix)
	
	# plot the result
	vis_x = vis_data[:, 0]
	vis_y = vis_data[:, 1]
	
	plt.scatter(vis_x, vis_y, c=colorLabels)
	plt.show()

plotData(similarityMatrix, labels)
