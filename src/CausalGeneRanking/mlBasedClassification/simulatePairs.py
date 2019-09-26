"""

	Make a fabriated dataset which we expect would be easy for our classifier. No noise, everything easy and perfect. 

"""

import random
import numpy as np

#1. Make windows, one will be negative, the other will be positive

windowCount = 2500
windowSize = 2000000
binSize = 1000
tadSize = 10000 #10kb
geneSize = 1000
svSize = 1000
enhancerCount = 5

#Make a bin map, where each possible position within the window is assigned a bin.
binMap = dict()
currentBin = 0

for pos in range(0, (windowSize*2)+1):
	
	binMap[pos] = currentBin
	if pos % binSize == 0: #new bin
		currentBin += 1

windowFeaturesSVs = []
windowFeaturesGenes = []
windowFeaturesTads = []
windowFeaturesEnhancers = []
windowLabels = []
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
	
	pairFeaturesGenes = np.zeros(int(windowSize*2/binSize+1))
	pairFeaturesTads = np.zeros(int(windowSize*2/binSize+1))
	pairFeaturesEnhancers = np.zeros(int(windowSize*2/binSize+1))
	
	#Gene	
	binRange = range(binMap[geneStart], binMap[geneEnd])
	for newBin in binRange:	
		pairFeaturesGenes[newBin] = 1
	
	#TADs
	tad1StartBin = binMap[tad1Start]
	tad1EndBin = binMap[tad1End]
	tad2EndBin = binMap[tad2End]
	pairFeaturesTads[tad1StartBin] = 1
	pairFeaturesTads[tad1EndBin] = 1
	pairFeaturesTads[tad2EndBin] = 1
	
	#Enhancers
	for pos in enhancersPos:
		enhBin = binMap[pos]
		pairFeaturesEnhancers[enhBin] = 1

	#positive SV
	#First for the positive pair
	pairFeaturesSVs = np.zeros(int(windowSize*2/binSize+1))
	
	#small issue: in our other approach, we are looking at 2 mb around the SV, here the SV is within that window. So how to deal? Go to bp approach immediately? 
	
	binRange = range(binMap[posSvStart], binMap[posSvEnd])
	for newBin in binRange:	
		pairFeaturesSVs[newBin] = 1
		
	windowFeaturesSVs.append(pairFeaturesSVs)
	windowFeaturesGenes.append(pairFeaturesGenes)
	windowFeaturesTads.append(pairFeaturesTads)
	windowFeaturesEnhancers.append(pairFeaturesEnhancers)
	windowLabels.append(1)
	
	#Repeat for the negative pair
	pairFeaturesSVs = np.zeros(int(windowSize*2/binSize+1))
	
	binRange = range(binMap[negSvStart], binMap[negSvEnd])
	for newBin in binRange:	
		pairFeaturesSVs[newBin] = 1
		
	windowFeaturesSVs.append(pairFeaturesSVs)
	windowFeaturesGenes.append(pairFeaturesGenes)
	windowFeaturesTads.append(pairFeaturesTads)
	windowFeaturesEnhancers.append(pairFeaturesEnhancers)
	windowLabels.append(0)

#6. Test some classifiers, with standard looCV. 
#Aggregate the features of the pairs depending on which classifier to use
windowFeaturesSVs = np.array(windowFeaturesSVs)
windowFeaturesGenes = np.array(windowFeaturesGenes)
windowFeaturesTads = np.array(windowFeaturesTads)
windowFeaturesEnhancers = np.array(windowFeaturesEnhancers)

#8. CNN
featureMatrix = np.dstack((np.array(windowFeaturesSVs), np.array(windowFeaturesGenes), np.array(windowFeaturesTads), np.array(windowFeaturesEnhancers)))


from sklearn.model_selection import train_test_split
from mcfly import modelgen, find_architecture, storage
from keras.models import load_model
import os

#use just 1 channel for now, later split into 2 and see if it makes a difference
X_train, X_test, y_train_list, y_test_list = train_test_split(featureMatrix, windowLabels, test_size=.4, random_state=42)

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

exit()



#7. Simple classifiers
featureMatrix = np.concatenate((np.array(windowFeaturesSVs), np.array(windowFeaturesGenes), np.array(windowFeaturesTads), np.array(windowFeaturesEnhancers)), axis=1) #stitch columns together

print(len(windowLabels))
print(featureMatrix.shape)
print(windowFeaturesSVs.shape)
print(windowFeaturesTads.shape)

def plotData():
		
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
	plt.imshow(plotGrid, cmap=cmap, interpolation='nearest')		
	plt.show()
	
	
	vis_data = bh_sne(featureMatrix)
	
	# plot the result
	vis_x = vis_data[:, 0]
	vis_y = vis_data[:, 1]
	
	plt.scatter(vis_x, vis_y, c=colorLabels)
	plt.show()

plotData()

#Random forest
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from scipy import interp

#Randomize labels
from random import shuffle

###very slow function
def loocv(featureMatrix, labels, clf):
	
	featureIndices = np.arange(featureMatrix.shape[0])
	labels = np.array(labels)
	
	loo = LeaveOneOut()
	
	scores = []
	aucs = []
	auprcs = []
	for train_index, test_index in loo.split(featureMatrix):
		
		X_train, X_test = featureMatrix[train_index], featureMatrix[test_index]
		y_train, y_test = labels[train_index], labels[test_index]
		
		
		# clf.fit(X_train, y_train)
		# 	
		# predictions = clf.predict(X_test)
		# predsDiff = np.average(y_test == np.sign(predictions))
		# scores.append(predsDiff)
		
		# preds = clf.predict_proba(X_test)[:,1]
		# fpr, tpr, thresholds = metrics.roc_curve(y_test, preds, pos_label=1)
		# aucS = metrics.auc(fpr, tpr)
		# aucs.append(aucS)
		# 
		# precision, recall, thresholds = precision_recall_curve(y_test, predictions)
		# aucScore = auc(recall, precision)
		# auprcs.append(aucScore)
		
	print("score: ", np.mean(scores))
	print("AUC: ", np.mean(aucs))
	print("AUPRC: ", np.mean(auprcs))	

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
	
	
	
	
				
print("Training classifiers")
from sklearn.ensemble import RandomForestClassifier
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
performCV(featureMatrix, windowLabels, rfClassifier)
#score = cross_val_score(rfClassifier, featureMatrix, windowLabels, cv=10)
#print("RF score: ", np.mean(score))
#loocv(featureMatrix, windowLabels, rfClassifier)

exit()

