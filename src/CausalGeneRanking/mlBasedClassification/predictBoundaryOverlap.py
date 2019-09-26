"""
	Test how we need to define features to be able to predict if SVs disrupt TAD boundaries. 

"""

from random import randrange
import random
import numpy as np


window = 2000000 #test with our 2 MB window
svSize = 500000 #use a 500 kb size for the SV for now. 

examples = 1000

features = []
labels = []

for exampleInd in range(0, examples):
	#Make 1000 examples where we have a random SV start and SV end, where the TAD boundary position is in the middle (any position).
	
	#The SV can start max window - size
	svStart = random.randint(0, (window-svSize))
	svEnd = svStart + svSize
	
	#randomly choose a TAD boundary within the SV.
	tadPosPositive = random.randint(svStart, svEnd)
	
	#keep these windows
	features.append([svStart, svEnd, tadPosPositive])
	labels.append(1)
	
	#Make 1000 examples where we have a random SV start and SV end, where the TAD boundary position is not inside this SV. 
	
	svStart = random.randint(0, (window-svSize))
	svEnd = svStart + svSize
	
	#Randomly choose a TAD boundary outside of the SV.
	#coin flip for which side of the SV
	ind = random.randint(0,1)
	
	if ind == 0: #left side
		tadPosNegative = random.randint(0,svStart)
	else:
		tadPosNegative = random.randint(svEnd, window)

	features.append([svStart, svEnd, tadPosNegative])
	labels.append(0)
	

features = np.array(features)

#Train a classifier

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
performCV(features, labels, rfClassifier)


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

plotData(features, labels)