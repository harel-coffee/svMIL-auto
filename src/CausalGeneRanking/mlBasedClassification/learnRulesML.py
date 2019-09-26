"""
	Use machine learning approaches to find the most informative features for detecting causal non-coding SVs. 

"""

import sys
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import math
from tsne import bh_sne

from inputParser import InputParser


from cleanlab.classification import LearningWithNoisyLabels
from cleanlab.noise_generation import generate_noise_matrix_from_trace
from cleanlab.noise_generation import generate_noisy_labels
from cleanlab.util import print_noise_matrix

#Take windowed SV-pairs and DEG pairs as input
print("loading pairs")
svGenePairs = np.loadtxt(sys.argv[1], dtype='object')
print("loading deg pairs")
svGeneDegPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')

print(svGeneDegPairs.shape)
print(svGenePairs.shape)


print("splitting pairs")
#Split the data into positive (all DEG pairs) and negative (pairs that are not in the DEG set)
positivePairs = svGeneDegPairs[:,0]
#negativePairs = np.setdiff1d(svGenePairs, svGeneDegPairs[:,0])

#setdiff is very slow in python 3 somehow, so switch to ultrafast dictionaries
posDict = dict()
for pair in positivePairs:
	posDict[pair] = 1
for pair in svGenePairs:
	if pair not in posDict:
		posDict[pair] = 0

negativePairs = []
for pair in posDict:
	if posDict[pair] == 0:
		negativePairs.append(pair)

negativePairs = np.array(negativePairs, dtype='object')

# print("filtering out translocations pos")
# #First focus only on intrachromosomal SVs.
# intraPairsPositive = []
# for pair in positivePairs:
# 	splitPair = pair.split("_")
# 	
# 	chr1 = splitPair[1]
# 	chr2 = splitPair[4]
# 	
# 	if chr1 == chr2:
# 		intraPairsPositive.append(pair)
# intraPairsNegative = []
# print("filtering out translocations neg")
# for pair in negativePairs:
# 	splitPair = pair.split("_")
# 	
# 	chr1 = splitPair[1]
# 	chr2 = splitPair[4]
# 	
# 	if chr1 == chr2:
# 		intraPairsNegative.append(pair)
# 
# positivePairs = np.array(intraPairsPositive)
# negativePairs = np.array(intraPairsNegative)

#Randomly subsample the SV-pairs to be balanced with the DEG pairs
np.random.seed(0)
negativePairsSubsampled = np.random.choice(negativePairs, positivePairs.shape[0])

print(positivePairs.shape)
print(negativePairsSubsampled.shape)

#allPairs = np.concatenate((positivePairs, negativePairsSubsampled))
#labels = np.array([1]*positivePairs.shape[0] + [0]*negativePairsSubsampled.shape[0])

allPairs = np.concatenate((positivePairs, negativePairs))
labels = np.array([1]*positivePairs.shape[0] + [0]*negativePairs.shape[0])

#split the pairs to be able to sort them
splitPairs = []
for pair in allPairs:
	splitPair = pair.split("_")
	splitPairs.append([splitPair[0], splitPair[1], int(splitPair[2]), int(splitPair[3]), splitPair[4], int(splitPair[5]), int(splitPair[6]), splitPair[7]])

splitPairs = np.array(splitPairs, dtype='object')
sortedPairs = splitPairs[splitPairs[:,1].argsort()]
sortedLabels = labels[splitPairs[:,1].argsort()] #TO DO here make labels as well with the right sorting

# 
# #Check features
# enhancerFeatures = np.load('features.npy')
# 
# print enhancerFeatures.shape
# #how many rows contain at least one?
# print np.unique(np.where(enhancerFeatures == 1)[0]).shape[0]
# 
# 
# vis_data = bh_sne(enhancerFeatures)
# 
# # plot the result
# vis_x = vis_data[:, 0]
# vis_y = vis_data[:, 1]
# 
# plt.scatter(vis_x, vis_y, c=sortedLabels)
# plt.show()
# 
# exit()
# 
# #First make a PCA and tSNE
# pca = PCA(n_components=2)
# projected = pca.fit_transform(enhancerFeatures)
# 
# colorLabels = []
# 
# for label in sortedLabels:
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
# exit()
# ###1. Assign features based on position in the nearby window
# 
# window = 2000000 #look in a 2 MB window
# binSize = 1000
# 
# #Make a bin map, where each possible position within the window is assigned a bin.
# binMap = dict()
# currentBin = 0
# 
# for pos in range(0, (window*2)+1):
# 	
# 	binMap[pos] = currentBin
# 	if pos % binSize == 0: #new bin
# 		currentBin += 1
# 
# enhancerData = InputParser().getEnhancersFromFile('../../../data/enhancers/enhancer_gene_GM12878.txt')
# 
# 
# #For every pair, check which eQTLs are in that window.
# windowFeaturesEnhancers = []
# 
# print("Making features for pairs")
# currentChr = None
# count = 0
# for pair in sortedPairs:
# 	
# 	pairFeatures = np.zeros(window*2/binSize+1)
# 
# 	sv = pair[1:len(pair)]
# 	
# 	#only get a new chr subset if the chr changes
# 	if sv[0] != currentChr:
# 		#Get the subset of eQTls on this chromosome
# 		chr1Subset = enhancerData[np.where(enhancerData[:,0] == sv[0])]
# 		print chr1Subset.shape
# 		currentChr = sv[0]
# 	
# 	#The only eQTLs we consider are the ones that are starting or ending within the window, and are not within the SV.
# 	#So, the start must be befor the SV breakpoint, the end after the SV bp-window, but end before the SV breakpoint.
# 	startMatches = (chr1Subset[:,1] <= sv[1]) * (chr1Subset[:,1] >= (sv[1] - window)) * (chr1Subset[:,2] <= sv[1])
# 			
# 	#The reverse for the end of the SV.
# 	endMatches = (chr1Subset[:,2] >= sv[5]) * (chr1Subset[:,2] <= (sv[5] + window)) * (chr1Subset[:,1] >= sv[5])
# 	
# 	matchingStart = chr1Subset[startMatches] # must be either matching on the left or right.
# 	matchingEnd = chr1Subset[endMatches] # must be either matching on the left or right.
# 	
# 	#For every position, determine the position relative to the SV and add it at the right place in the feature vector.
# 	for match in matchingStart:
# 		#the position in the feature vector is 2 MB - abs(the SV start - element pos)
# 		featurePos = window - abs(match[1] - int(sv[1]))
# 		
# 		correctBin = binMap[featurePos]
# 		pairFeatures[correctBin] = 1
# 	for match in matchingEnd:
# 		#Here, we do + 2MB
# 		featurePos = window + abs(match[1] - int(sv[5]))
# 		correctBin = binMap[featurePos]
# 		pairFeatures[correctBin] = 1
# 	
# 	windowFeaturesEnhancers.append(pairFeatures)
# 	
# 
# print np.array(windowFeaturesEnhancers)
# np.save('features.npy', windowFeaturesEnhancers)
# 
# exit()
# 


###2. Way to assign features using the current rules

#Assign features to each pair (gains/losses of elements yes/no)
#Load the rule-based SV-gene pairs, and get the gains & losses from there.
svGenePairsRules = np.loadtxt(sys.argv[3], dtype='object')

print("no of rules: ", svGenePairsRules.shape)

positivePairsFeatures = []
positiveWithFeatures = 0
positiveWithoutFeatures = 0
svFeaturesPos = dict()
for pair in positivePairs:
	if pair in svGenePairsRules[:,0]:
		positiveWithFeatures += 1
		#get these features
		features = svGenePairsRules[svGenePairsRules[:,0] == pair,:][0]
		#skip the pair name and also the total score
		
		positivePairsFeatures.append(features)
		splitPair = pair.split("_")

		sv = "_".join(splitPair[1:len(splitPair)])
		
		if sv not in svFeaturesPos:
			svFeaturesPos[sv] = []
		svFeaturesPos[sv].append(list(features[1:len(features)-1].astype('float')))
		
	# 	features = features[1:len(features)-1]
	# 	features = [float(feature) for feature in features]
	# 	positivePairsFeatures.append(features)
	# else: #if not, assign all features as 0
	# 	positiveWithoutFeatures += 1
	# 	positivePairsFeatures.append([0]*26)

#repeat for negative pairs
negativePairsFeatures = []
negativeWithFeatures = 0
negativeWithoutFeatures = 0
svFeaturesNeg = dict()
for pair in negativePairsSubsampled:
#for pair in negativePairs:
	if pair in svGenePairsRules[:,0]:
		negativeWithFeatures += 1
		#get these features
		features = svGenePairsRules[svGenePairsRules[:,0] == pair,:][0]
		#skip the pair name and also the total score
		
		negativePairsFeatures.append(features)
		
		splitPair = pair.split("_")
		
		sv = "_".join(splitPair[1:len(splitPair)])
		
		if sv not in svFeaturesNeg:
			svFeaturesNeg[sv] = []

		svFeaturesNeg[sv].append(list(features[1:len(features)-1].astype('float')))
		
	# 	features = features[1:len(features)-1]
	# 	features = [float(feature) for feature in features]
	# 	negativePairsFeatures.append(features)
	# else: #if not, assign all features as 0
	# 	negativeWithoutFeatures += 1
	# 	negativePairsFeatures.append([0]*26)

positivePairsFeatures = np.array(positivePairsFeatures, dtype='object')
negativePairsFeatures = np.array(negativePairsFeatures, dtype='object')

#Try MIL

#Make bags
bags = []
labels = []
for sv in svFeaturesPos:
	bags.append(svFeaturesPos[sv])
	labels.append(1)

for sv in svFeaturesNeg:
	bags.append(svFeaturesNeg[sv])
	labels.append(0)

bags = np.array(bags)
instances = np.vstack(bags)
labels = np.array(labels)


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

		currentMinDistance = np.min(summedDistance)
		if currentMinDistance < np.min(minDistance):
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


exit()




#output to files
np.savetxt('degPairsFeatures.txt', positivePairsFeatures, fmt='%s', delimiter='\t')
np.savetxt('nonDegPairsFeatures.txt', negativePairsFeatures, fmt='%s', delimiter='\t')

print(positivePairsFeatures)
print(negativePairsFeatures)

print(positiveWithFeatures)
print(positiveWithoutFeatures)
print(negativeWithFeatures)
print(negativeWithoutFeatures)

#For each pair that does not match, are the SV and gene > 2 MB away from each other?
#Get the genes and positions
from inputParser import InputParser
import settings
causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set. 
geneData = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#which rule-based pairs are not found by the windowed approach? 
farAway = 0
missed = 0
for pair in svGenePairsRules[:,0]:
	if pair not in svGenePairs:
		missed += 1
		splitPair = pair.split("_")
		gene = geneData[geneData[:,3] == splitPair[0]][0]
		
		startStartDist = abs(int(splitPair[2]) - gene[1])
		startEndDist = abs(int(splitPair[2]) - gene[2])
		endStartDist = abs(int(splitPair[6]) - gene[1])
		endEndDist = abs(int(splitPair[6]) - gene[2])
		
		#get the minimum distance
		minVal = np.min([startStartDist, startEndDist, endStartDist, endEndDist])
		if minVal >= 2000000:
			farAway += 1
		else:
			print("why am I here?")
			print(pair)
			print(gene)

print("faraway: ", farAway)
print('missed', missed)


allFeatures = np.concatenate((positivePairsFeatures, negativePairsFeatures))
labels = [1]*positivePairsFeatures.shape[0] + [0]*negativePairsFeatures.shape[0]
# 
# #First make a PCA and tSNE
# pca = PCA(n_components=2)
# projected = pca.fit_transform(allFeatures)
# 
colorLabels = []

for label in labels:
	
	if label == 1:
		colorLabels.append('r')
	else:
		colorLabels.append('b')

# fig,ax=plt.subplots(figsize=(7,5))
# plt.scatter(projected[:, 0], projected[:, 1], c=colorLabels)
# plt.show()



vis_data = bh_sne(allFeatures)

# plot the result
vis_x = vis_data[:, 0]
vis_y = vis_data[:, 1]

plt.scatter(vis_x, vis_y, c=colorLabels)
plt.show()

# exit()
# 
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
# boxWidth = 0.1
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
# 

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
#shuffle(labels)

#print(labels)

labels = np.array(labels)
X_train, X_test, y_train, y_test = train_test_split(allFeatures, labels, test_size=0.4, random_state=42)

#Implement leave-one-patient out CV for all classifiers, make a function
def loopoCV(sortedPairs, sortedLabels, allFeatures, clf, aucBool= True):
	
	#gather all unique patients
	patients = np.unique(sortedPairs[:,7])
	patientScores = []
	patientAucs = []
	patientAuprcs = []
	posZero = 0
	for patient in patients:
		#gather a training set for all but this patient
		X_train = allFeatures[sortedPairs[:,7] != patient]
		y_train = sortedLabels[sortedPairs[:,7] != patient]
		
		X_test = allFeatures[sortedPairs[:,7] == patient]
		y_test = sortedLabels[sortedPairs[:,7] == patient]
		
		
		trainFeatPos = X_train[y_train == 1]
		testFeatPos = X_test[y_test == 1]
		
		#Train classifier on each fold.
		clf.fit(X_train, y_train)

		predictions = rfClassifier.predict(X_test)
		predsDiff = np.average(y_test == np.sign(predictions))
		patientScores.append(predsDiff)
		
		if aucBool == True: #not available for all classifiers
			preds = rfClassifier.predict_proba(X_test)[:,1]
			fpr, tpr, thresholds = metrics.roc_curve(y_test, preds, pos_label=1)
			aucS = metrics.auc(fpr, tpr)
			patientAucs.append(aucS)
		
		precision, recall, thresholds = precision_recall_curve(y_test, predictions)
		aucScore = auc(recall, precision)
		patientAuprcs.append(aucScore)
	
	#report averages and std
	print("Average CV score: ", np.mean(patientScores))
	print("std CV score: ", np.std(patientScores))
	if aucBool == True:
		print("Average CV AUC: ", np.mean(patientAucs))
		print("std CV AUC: ", np.std(patientAucs))
		
	print("Average CV AUPRC: ", np.mean(patientAuprcs))
	print("std CV AUPRC: ", np.std(patientAuprcs))
	
	
#Try each classifier
print("Random forest")
from sklearn.ensemble import RandomForestClassifier
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
loopoCV(sortedPairs, sortedLabels, allFeatures, rfClassifier, True)

print("lasso")

currentAlpha = 1e-2
lasso = Lasso(alpha=currentAlpha)
loopoCV(sortedPairs, sortedLabels, allFeatures, lasso, False)

print("linear SVC")

clf = LinearSVC()
loopoCV(sortedPairs, sortedLabels, allFeatures, clf, False)

print("cleaned labels RF")
rp = LearningWithNoisyLabels(rfClassifier)
loopoCV(sortedPairs, sortedLabels, allFeatures, rp, False)

#print("cleaned labels lasso")
#rp = LearningWithNoisyLabels(lasso)
#loopoCV(sortedPairs, sortedLabels, allFeatures, rp, False)

print("cleaned labels SVC")
rp = LearningWithNoisyLabels(clf)
loopoCV(sortedPairs, sortedLabels, allFeatures, rp, False)


exit()



#1. Random forest
from sklearn.ensemble import RandomForestClassifier
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
rfClassifier.fit(X_train, y_train) #Use the bag labels, not the instance labels

predictions = rfClassifier.predict(X_test)
predsDiff = np.average(y_test == np.sign(predictions))
print("RF score: ", predsDiff)

preds = rfClassifier.predict_proba(X_test)[:,1]
fpr, tpr, thresholds = metrics.roc_curve(y_test, preds, pos_label=1)
aucS = metrics.auc(fpr, tpr)
print("RF AUC: ", aucS)

precision, recall, thresholds = precision_recall_curve(y_test, predictions)
aucScore = auc(recall, precision)
print("RF AUPRC: ", aucScore)

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
print("lasso AUPRC: ", aucScore)

#3. SVM


clf = LinearSVC()
clf.fit(X_train, y_train)
score = clf.score(X_test, y_test)
preds = clf.predict(X_test)
predsDiff = np.average(y_test == np.sign(preds))
print("SVM score: ", score)

precision, recall, thresholds = precision_recall_curve(y_test, preds)
aucScore = auc(recall, precision)

print("SVM AUPRC: ", aucScore)


#test noisy labels

rp = LearningWithNoisyLabels(rfClassifier)
rp.fit(X_train, y_train)
score = rp.clf.score(X_test, y_test)
print("Nosy labels score: ", score)

preds = rp.predict_proba(X_test)[:,1]
fpr, tpr, thresholds = metrics.roc_curve(y_test, preds, pos_label=1)
aucS = metrics.auc(fpr, tpr)
print("Noisy labels AUC: ", aucS)

predictions = rp.predict(X_test)
precision, recall, thresholds = precision_recall_curve(y_test, predictions)
aucScore = auc(recall, precision)

print("Noisy labels AUPRC: ", aucScore)
