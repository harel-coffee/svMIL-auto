"""
	Use a MIL-based representation to determine if an SV overlaps a TAD boundary or not. 

"""

import random
import numpy as np

#let's start with a very simple simulated case to see if our idea works.
#We make bags representing windows, and instances representing the elements inside this window (for now, SVs and TAD boundaries)

#For a positive window, we model a SV overlapping a TAD boundary
#For a negative window, we model an SV never overlapping any TAD boundary

#The instances will for now have simple features: their positions. A start position, and an end postion.
#First, we can use 1 TAD boundary, so there will be 2 instances for each window: a TAD boundary, and an SV.

### For boundary overlap prediction

# 
# window = 2000000 #test with our 2 MB window
# svSize = 500000 #use a 500 kb size for the SV for now. 
# 
# examples = 1000
# 
# bags = []
# labels = []
# 
# for exampleInd in range(0, examples):
# 	#Make 1000 examples where we have a random SV start and SV end, where the TAD boundary position is in the middle (any position).
# 	
# 	#The SV can start max window - size
# 	svStart = random.randint(0, (window-svSize))
# 	svEnd = svStart + svSize
# 	
# 	svInstance = [svStart, svEnd]
# 	
# 	#randomly choose a TAD boundary within the SV.
# 	tadPosPositive = random.randint(svStart, svEnd)
# 	
# 	tadInstance = [tadPosPositive, tadPosPositive] #TAD boundaries start and end at the same bp.
# 	
# 	
# 	#keep these windows,store the instances in the bag.
# 	instances = [svInstance, tadInstance]
# 	bags.append(instances)
# 	labels.append(1)
# 	
# 	#Make 1000 examples where we have a random SV start and SV end, where the TAD boundary position is not inside this SV. 
# 	
# 	svStart = random.randint(0, (window-svSize))
# 	svEnd = svStart + svSize
# 	
# 	svInstance = [svStart, svEnd]
# 	
# 	#Randomly choose a TAD boundary outside of the SV.
# 	#coin flip for which side of the SV
# 	ind = random.randint(0,1)
# 	
# 	if ind == 0: #left side
# 		tadPosNegative = random.randint(0,svStart)
# 	else:
# 		tadPosNegative = random.randint(svEnd, window)
# 
# 	tadInstance = [tadPosNegative, tadPosNegative] #TAD boundaries start and end at the same bp.
# 
# 	instances = [svInstance, tadInstance]
# 	bags.append(instances)
# 	labels.append(0)

###For the full case with enhancer disruptions

windowCount = 2500
windowSize = 2000000
binSize = 1000
tadSize = 10000 #10kb
geneSize = 1000
svSize = 1000
enhancerCount = 5

bags = []
labels = []

for window in range(0, windowCount):

	#2. Randomly place 2 TADs within this window
	#The minimum place to start is at the beginning, and we need to fit 2*tad size (back to back), so the maximum start is window size - 2*tad size
	minimumStart = 0
	maximumStart = windowSize - (2*tadSize)
	
	tad1Start = random.randint(minimumStart, maximumStart)
	tad1End = tad1Start + tadSize
	tad2Start = tad1End
	tad2End = tad2Start + tadSize

	#3. Set a negative SV in 1 window, and a positive in the other
	#the SV must be within either TAD, do a coin flip for that.
	tadNumber = random.randint(0,1)
	
	#Then select a random position within the TAD, but not disrupting the boundary.
	if tadNumber == 0:
		negSvStart = random.randint(tad1Start + 1, tad1End - 1 - svSize)
		negSvEnd = negSvStart + svSize
	else:
		negSvStart = random.randint(tad2Start + 1, tad2End - 1 - svSize)
		negSvEnd = negSvStart + svSize
	
	#The positive SV overlaps the boundary. 
	posSvStart = tad1End - (svSize/2)
	posSvEnd = posSvStart + svSize
	
	#4. Decide based on SV placement where the gene and enhancers will be placed
	#It may be easy to fix the gene placement based on where the SV is not. That may make it quite easy for the classifier, but we can fix that later.
	if tadNumber == 0:
		#place the gene in the second TAD.
		geneStart = random.randint(tad2Start, tad2End)
		geneEnd = geneStart + geneSize
		#place the enhancers randomly around the negative SV in the first TAD, make sure that there is no overlap for now.
		enhancerWindow = list(np.arange(tad1Start, negSvStart))
		enhancerWindow += list(np.arange(negSvEnd, posSvStart))
		enhancersPos = np.random.choice(enhancerWindow, enhancerCount, replace=False)
	else:
		#place the gene in the first TAD.
		geneStart = random.randint(tad1Start, tad1End)
		geneEnd = geneStart + geneSize
		#place the enhancers randomly around the negative SV in the second TAD, make sure that there is no overlap for now.
		enhancerWindow = list(np.arange(posSvEnd, negSvStart)) #exclude the positive SV such that the enhancers never overlap with the SV. 
		enhancerWindow += list(np.arange(negSvEnd, tad2End))
		enhancersPos = np.random.choice(enhancerWindow, enhancerCount, replace=False)
	
	#5. Make feature descriptions like we do for the other windows (same bin size)
	
	posInstances = []
	negInstances = []
	
	geneInstance = [geneStart, geneEnd]
	tad1Instance = [tad1Start, tad1End]
	tad2Instance = [tad2Start, tad2End]
	posSvInstance = [posSvStart, posSvEnd]
	negSvInstance = [negSvStart, negSvEnd]
	
	posInstances = [geneInstance, tad1Instance, tad2Instance, posSvInstance]
	negInstances = [geneInstance, tad1Instance, tad2Instance, negSvInstance]
	
	for enhancerPos in enhancersPos:
		posInstances.append([enhancerPos, enhancerPos])
		negInstances.append([enhancerPos, enhancerPos])
	
	bags.append(posInstances)
	labels.append(1)
	bags.append(negInstances)
	labels.append(0)


bags = np.array(bags)
instances = np.vstack(bags)


from random import shuffle
#shuffle(labels)

labels = np.array(labels)

print(bags)
print(instances)

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
