"""
	Use machine learning approaches to find the most informative features for detecting causal non-coding SVs. 

"""

import sys
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import math
from tsne import bh_sne

#Take windowed SV-pairs and DEG pairs as input
svGenePairs = np.loadtxt(sys.argv[1], dtype='object')
svGeneDegPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')


#Split the data into positive (all DEG pairs) and negative (pairs that are not in the DEG set)
positivePairs = svGeneDegPairs[:,0]
negativePairs = np.setdiff1d(svGenePairs, svGeneDegPairs[:,0])

#Randomly subsample the SV-pairs to be balanced with the DEG pairs
np.random.seed(0)
negativePairsSubsampled = np.random.choice(negativePairs, positivePairs.shape[0])

print positivePairs.shape
print negativePairsSubsampled.shape

#Assign features to each pair (gains/losses of elements yes/no)
#Load the rule-based SV-gene pairs, and get the gains & losses from there.
svGenePairsRules = np.loadtxt(sys.argv[3], dtype='object')
# 
# missing = 0
# present = 0
# trans = 0
# for pair in svGenePairsRules[:,0]:
# 	if pair not in svGenePairs:
# 		
# 		missing += 1
# 	else:
# 		present += 1
# 
# print missing
# print trans
# print present
# exit()

positivePairsFeatures = []
positiveWithFeatures = 0
positiveWithoutFeatures = 0
for pair in positivePairs:
	if pair in svGenePairsRules[:,0]:
		positiveWithFeatures += 1
		#get these features
		features = svGenePairsRules[svGenePairsRules[:,0] == pair,:][0]
		#skip the pair name and also the total score
		features = features[1:len(features)-1]
		features = [float(feature) for feature in features]
		positivePairsFeatures.append(features)
	else: #if not, assign all features as 0
		positiveWithoutFeatures += 1
		positivePairsFeatures.append([0]*26)

#repeat for negative pairs
negativePairsFeatures = []
negativeWithFeatures = 0
negativeWithoutFeatures = 0
for pair in negativePairsSubsampled:
	if pair in svGenePairsRules[:,0]:
		negativeWithFeatures += 1
		#get these features
		features = svGenePairsRules[svGenePairsRules[:,0] == pair,:][0]
		#skip the pair name and also the total score
		features = features[1:len(features)-1]
		features = [float(feature) for feature in features]
		negativePairsFeatures.append(features)
	else: #if not, assign all features as 0
		negativeWithoutFeatures += 1
		negativePairsFeatures.append([0]*26)

positivePairsFeatures = np.array(positivePairsFeatures)
negativePairsFeatures = np.array(negativePairsFeatures)

print positivePairsFeatures
print negativePairsFeatures

print positiveWithFeatures
print positiveWithoutFeatures
print negativeWithFeatures
print negativeWithoutFeatures

allFeatures = np.concatenate((positivePairsFeatures, negativePairsFeatures))
labels = [1]*positivePairsFeatures.shape[0] + [0]*negativePairsFeatures.shape[0]

#First make a PCA and tSNE
# pca = PCA(n_components=2)
# projected = pca.fit_transform(allFeatures)
# 
# colorLabels = []
# 
# for label in labels:
# 	
# 	if label == 1:
# 		colorLabels.append('r')
# 	else:
# 		colorLabels.append('b')
# 
# fig,ax=plt.subplots(figsize=(7,5))
# plt.scatter(projected[:, 0], projected[:, 1], c=colorLabels)
# plt.show()
# 
# 
# 
# vis_data = bh_sne(allFeatures)
# 
# # plot the result
# vis_x = vis_data[:, 0]
# vis_y = vis_data[:, 1]
# 
# plt.scatter(vis_x, vis_y, c=colorLabels)
# plt.show()
# 
# exit()

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
# boxWidth = 1
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
# plt.imshow(plotGrid, cmap=cmap, interpolation='nearest')		
# plt.show()


#### Try classification

#Make random training/test subset, later use cross-validation
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import Lasso
from sklearn.metrics import auc, precision_recall_curve
from sklearn.svm import LinearSVC

#Randomize labels
from random import shuffle
shuffle(labels)

print labels

X_train, X_test, y_train, y_test = train_test_split(allFeatures, labels, test_size=0.4, random_state=42)

#1. Random forest
from sklearn.ensemble import RandomForestClassifier
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
rfClassifier.fit(X_train, y_train) #Use the bag labels, not the instance labels

predictions = rfClassifier.predict(X_test)
predsDiff = np.average(y_test == np.sign(predictions))
print("RF score: ", predsDiff)

precision, recall, thresholds = precision_recall_curve(y_test, predictions)
aucScore = auc(recall, precision)
print("RF AUC: ", aucScore)

#2. Lasso

currentAlpha = 1e-2
lasso = Lasso(alpha=currentAlpha)
lasso.fit(X_train,y_train)

test_score=lasso.score(X_test,y_test)
coeff_used = np.sum(lasso.coef_!=0)
preds = lasso.predict(X_test)
predsDiff = np.average(y_test == np.sign(preds))
print("lasso score: ", predsDiff)

precision, recall, thresholds = precision_recall_curve(y_test, preds)
aucScore = auc(recall, precision)
print("lasso AUC: ", aucScore)

#3. SVM


clf = LinearSVC()
clf.fit(X_train, y_train)
score = clf.score(X_test, y_test)
preds = clf.predict(X_test)
predsDiff = np.average(y_test == np.sign(preds))
print("SVM score: ", score)

precision, recall, thresholds = precision_recall_curve(y_test, preds)
aucScore = auc(recall, precision)

print("SVM AUC: ", aucScore)
