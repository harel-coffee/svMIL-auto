"""

	Make features matrices for all features in a 1 MB window.
	Try training a deep learning model

"""

import numpy as np
import sys
from inputParser import InputParser
from featureMatrixDefiner import FeatureMatrixDefiner
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import math
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.linear_model import Lasso
from sklearn.metrics import auc, precision_recall_curve
from sklearn.svm import LinearSVC
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from scipy import interp
import pickle as pkl

#1. Process the data

#Get the bags
with open(sys.argv[1], 'rb') as handle:
	bagDict = pkl.load(handle)

gl = False
if gl == True:
	with open(sys.argv[4], 'rb') as handle:
		glBagDict = pkl.load(handle)



#determine the bag labels given a file of DEG pairs
degPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')
#pathwayAnnotation = np.loadtxt(sys.argv[3], dtype='object')

print("initial bag num: ", len(bagDict))
print('initial deg pairs: ', degPairs.shape[0])

#feature selection, increase number of features
featureCount = len(bagDict[list(bagDict.keys())[0]][0])
featureStart = featureCount #set this to featureCount to run with all features.
similarityMatrices = dict()
bagLabels = []
for featureInd in range(featureStart, featureCount+1):
	
	bags = []
	bagLabels = []
	posCount = 0
	negCount = 0
	positiveBags = []
	negativeBags = []
	removedPathwayPairs = 0
	svType = ''
	for pair in bagDict:
		
		if svType != '':
			splitPair = pair.split("_")
			if splitPair[8] != svType:
				continue
		
		#get the label of the bag by checking if it exists in degPairs
		if pair in degPairs[:,0]:
			#if pair in pathwayAnnotation[:,0]: #skip the ones that have possible pathway effects
			#	removedPathwayPairs += 1
				#continue
			bagLabels.append(1)
			posCount += 1
			#get the right number of features per instance
			instances = []
			for instance in bagDict[pair]:
				instances.append(instance[0:featureInd+1])
				
			positiveBags.append(instances)
			#positiveBags.append(bagDict[pair])
		else:
			bagLabels.append(0)
			negCount += 1
			
			#get the right number of features per instance
			instances = []
			for instance in bagDict[pair]:

				instances.append(instance[0:featureInd+1])
			#print(instances)	
			negativeBags.append(instances)
			
			#negativeBags.append(bagDict[pair])
	
		#bags.append(bagDict[pair])
	
	if gl == True:
		for bag in glBagDict:
			negativeBags.append(glBagDict[bag])
	
	print('removed pathway pairs: ', removedPathwayPairs)
	positiveBags = np.array(positiveBags)
	negativeBags = np.array(negativeBags)
	
	#take a random subset for speed
	#positiveBags = np.random.choice(positiveBags, 2000)
	#negativeBags = np.random.choice(negativeBags, 2000)
	
	#print(positiveBags.shape)
	#print(negativeBags.shape)
	#exit()
	print("positive bags: ", positiveBags.shape[0])
	print("negative bags: ", negativeBags.shape[0])
	
	#np.random.seed(0)
	#negativeBagsSubsampled = np.random.choice(negativeBags, posCount)
	
	#bags = np.array(bags)
	#bags = np.concatenate((positiveBags, negativeBagsSubsampled))
	#bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBagsSubsampled.shape[0])
	bags = np.concatenate((positiveBags, negativeBags))
	bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBags.shape[0])
	
	instances = np.vstack(bags)
	
	print(instances.shape)
	print(bagLabels.shape)

	#Make similarity matrix
	
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
		
		instanceAvg = np.mean(instanceSubset, axis=0)
		
		#compute distance to all other instances
		distance = np.abs(instanceAvg - instances)
		
		summedDistance = np.sum(distance,axis=1)
		similarityMatrix[bagInd,:] = summedDistance
		continue
		
		
		
		#Compute the pairwise distance matrix here
		# minDistance = float("inf")
		# minDistanceInd = 0
		# for instanceInd in range(0, instanceSubset.shape[0]):
		# 	instance = instanceSubset[instanceInd]
		# 	distance = np.abs(instance - instances) #compute the distances to the train instances, otherwise we are not in the same similarity space. 
		# 
		# 	#distance = np.abs(instance - otherInstances)
		# 
		# 	summedDistance = np.sum(distance,axis=1)
		# 
		# 	currentMinDistance = np.mean(summedDistance)
		# 	if currentMinDistance < np.mean(minDistance):
		# 		minDistance = summedDistance
		# 		minDistanceInd = instanceInd
		# 
		#This instance will be used as representative for this bag. We use this value as the similarity to all other instances.  
		similarityMatrix[bagInd] = minDistance
	
	print(similarityMatrix)
	print(bagLabels)
	similarityMatrices[featureInd] = similarityMatrix
	print(featureInd)

pca = PCA(n_components=2)
projected = pca.fit_transform(similarityMatrices[featureCount])
projectedWithOffset = projected

for row in range(0, projected.shape[0]):
	for col in range(0, projected.shape[1]):
		projectedWithOffset[row][col] += np.random.normal(-1, 1) * 0.5
		
projected = projectedWithOffset

colorLabels = []

for label in bagLabels:
	
	if label == 1:
		colorLabels.append('r')
	else:
		colorLabels.append('b')

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

boxWidth = 2
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

from random import shuffle
#shuffle(bagLabels)

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from inspect import signature
from sklearn import model_selection
from sklearn.metrics import average_precision_score
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score

f, axes = plt.subplots()	
def cvClassification(similarityMatrix, bagLabels, clf, color, plot, aps, f1s):
	
	scoring = {'accuracy' : make_scorer(accuracy_score), 
			   'precision' : make_scorer(precision_score),
			   'recall' : make_scorer(recall_score), 
			   'f1_score' : make_scorer(f1_score),
			   'average_precision' : make_scorer(average_precision_score)}
	
	kfold = model_selection.StratifiedKFold(n_splits=10)
	
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
	
	print('APs: ', results['test_average_precision'])
	print('Precisions: ', results['test_precision'])
	print('Recalls: ', results['test_recall'])
	return 0, 0
	#plot PR curves
	
	#f, axes = plt.subplots()	
	cv = StratifiedKFold(n_splits=10)
	
	accs = []
	aucs = []
	coeffs = []
	predDiffs = []
	thresholds = []
	precisions = []
	recalls = []
	i = -1
	for train, test in cv.split(similarityMatrix, bagLabels):
		#f, axes = plt.subplots()
		i+=1
		clf.fit(similarityMatrix[train], bagLabels[train]) #Use the bag labels, not the instance labels
	
		predictions = clf.predict_proba(similarityMatrix[test])
		precision, recall, thresholds = precision_recall_curve(bagLabels[test], predictions[:,1])
		preds = clf.predict(similarityMatrix[test])

		print(i)
		print(thresholds)
		print(len(precision))
		print(len(recall))
		recalls.append(recall)
		precisions.append(precision)
		
		
		#lab = 'Fold %d AP=%.4f' % (i+1, average_precision_score(bagLabels[test], predictions[:,1]))
		# print('fold: ', i, ':', average_precision_score(bagLabels[test], predictions[:,1]))
		# print(precision)
		# print(recall)
		#plt.plot(recall, precision)
		#axes.step(recall, precision, label=lab)
		#axes.set_xlabel('Recall')
		#axes.set_ylabel('Precision')
		#axes.legend(loc='lower right', fontsize='small')
		#y_real.append(bagLabels[test])
		#y_proba.append(predictions[:,1])
		#y_pred.append(preds)

		# step_kwargs = ({'step': 'post'} if 'step' in signature(plt.fill_between).parameters else {})
		# plt.step(recall, precision, color='b', alpha=0.2, where='post')
		# plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)
		#  
		# plt.xlabel('Recall')
		# plt.ylabel('Precision')
		# plt.ylim([0.0, 1.05])
		# plt.xlim([0.0, 1.0])
		# plt.title('2-class Precision-Recall curve: AP={0:0.2f}'.format(average_precision_score(bagLabels[test], predictions[:,1])))
		# plt.show()
		
	
	
	# y_real = np.concatenate(y_real)
	# y_proba = np.concatenate(y_proba)
	# y_pred = np.concatenate(y_pred)
	# print("overall CV AP: ", average_precision_score(y_real, y_proba))
	# print("overall CV F1: ", f1_score(y_real, y_pred))
	# aps.append(average_precision_score(y_real, y_proba))
	# f1s.append( f1_score(y_real, y_pred))

	# if plot == True:
	# 	precision, recall, thresholds = precision_recall_curve(y_real, y_proba)
	# 	lab = 'Overall AP=%.4f' % average_precision_score(y_real, y_proba)
	# 	axes.step(recall, precision, label=lab, lw=2, color=color)
	# 	axes.set_xlabel('Recall')
	# 	axes.set_ylabel('Precision')
	# 	axes.legend(loc='lower right', fontsize='small')
	
	return aps, f1s


rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)

#For each feature round, save the f1 and AP scores
#Then at the end, show the curves only for the last round (all features)
#show a separate figure for the feature selection scores
aps = []
f1s = []
for featureInd in range(featureStart, featureCount+1):
	
	plot = False
	if featureInd == featureCount:
		plot = True

	aps, f1s = cvClassification(similarityMatrices[featureInd], bagLabels, rfClassifier, 'black', plot, aps, f1s)

#print("All APS: ", aps)
#print("All F1S: ", f1s)

#repeat for shuffled labels, do this separately to make sure that labels are only shuffled from here
shuffle(bagLabels)
shuffledAps = []
shuffledF1s = []
print("shuffled")
for featureInd in range(featureStart, featureCount+1):

	plot = False
	if featureInd == featureCount:
		plot = True
	shuffledAps, shuffledF1s = cvClassification(similarityMatrices[featureInd], bagLabels, rfClassifier, 'red', plot, shuffledAps, shuffledF1s)

#print("Shuffled APs: ", shuffledAps)
#print("Shuffled F1s: ", shuffledF1s)
exit()
plt.show()

plt.clf()
plt.title('Non-shuffled APs (black) and F1s (blue)')
plt.plot(aps, color='black')
plt.plot(f1s, color='blue')
plt.ylim([0.5, 0.9])
plt.show()
plt.clf()

plt.title('Shuffled APs (black) and F1s (blue)')
plt.plot(shuffledAps, color='black')
plt.plot(shuffledF1s, color='blue')
plt.ylim([0.5, 0.9])
plt.show()
plt.clf()

exit()
# 
# 
# 
# vis_data = bh_sne(matrix)
# 
# # plot the result
# vis_x = vis_data[:, 0]
# vis_y = vis_data[:, 1]
# 
# plt.scatter(vis_x, vis_y, c=colorLabels)
# plt.show()
# 
# exit()


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
		
		if len(testFeatPos) == 0 or len(trainFeatPos) == 0:
			continue #not enough positive samples for classification
		
		
		#Train classifier on each fold.
		clf.fit(X_train, y_train)

		predictions = clf.predict(X_test)
		predsDiff = np.average(y_test == np.sign(predictions))
		patientScores.append(predsDiff)
		
		if aucBool == True: #not available for all classifiers
			preds = clf.predict_proba(X_test)[:,1]
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
	
	
	

#Try each classifier
print("Random forest")
from sklearn.ensemble import RandomForestClassifier
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
#loopoCV(pairs, labels, matrix, rfClassifier, True)
performCV(matrix, labels, rfClassifier)

print("lasso")

currentAlpha = 1e-2
lasso = Lasso(alpha=currentAlpha)
#loopoCV(pairs, labels, matrix, lasso, False)
#performCV(matrix, labels, lasso)
print("linear SVC")

#clf = LinearSVC()
#loopoCV(pairs, labels, matrix, clf, False)

print("cleaned labels RF")
rp = LearningWithNoisyLabels(rfClassifier)
#loopoCV(pairs, labels, matrix, rp, False)
performCV(matrix, labels, rp)

#print("cleaned labels lasso")
#rp = LearningWithNoisyLabels(lasso)
#loopoCV(sortedPairs, sortedLabels, allFeatures, rp, False)
# 
# print("cleaned labels SVC")
# rp = LearningWithNoisyLabels(clf)
# loopoCV(pairs, labels, matrix, rp, False)

exit()

# 
# 
# from sklearn.model_selection import train_test_split
# from sklearn import metrics
# from sklearn.metrics import roc_curve, auc
# from sklearn.linear_model import Lasso
# from sklearn.metrics import auc, precision_recall_curve
# from sklearn.svm import LinearSVC
# 
# #Randomize labels
# from random import shuffle
# #shuffle(labels)
# 
# print labels
# 
# X_train, X_test, y_train, y_test = train_test_split(matrix, labels, test_size=0.4, random_state=42)
# 
# #1. Random forest
# from sklearn.ensemble import RandomForestClassifier
# rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
# rfClassifier.fit(X_train, y_train) #Use the bag labels, not the instance labels
# 
# predictions = rfClassifier.predict(X_test)
# predsDiff = np.average(y_test == np.sign(predictions))
# print("RF score: ", predsDiff)
# 
# precision, recall, thresholds = precision_recall_curve(y_test, predictions)
# aucScore = auc(recall, precision)
# print("RF AUC: ", aucScore)
# 
# #2. Lasso
# 
# currentAlpha = 1e-2
# lasso = Lasso(alpha=currentAlpha)
# lasso.fit(X_train,y_train)
# 
# test_score=lasso.score(X_test,y_test)
# coeff_used = np.sum(lasso.coef_!=0)
# preds = lasso.predict(X_test)
# predsDiff = np.average(y_test == np.sign(preds))
# print("lasso score: ", predsDiff)
# 
# precision, recall, thresholds = precision_recall_curve(y_test, preds)
# aucScore = auc(recall, precision)
# print("lasso AUC: ", aucScore)
# 
# #3. SVM
# 
# 
# clf = LinearSVC()
# clf.fit(X_train, y_train)
# score = clf.score(X_test, y_test)
# preds = clf.predict(X_test)
# predsDiff = np.average(y_test == np.sign(preds))
# print("SVM score: ", score)
# 
# precision, recall, thresholds = precision_recall_curve(y_test, preds)
# aucScore = auc(recall, precision)
# 
# print("SVM AUC: ", aucScore)
# 
# 
# 
# 
# 
# exit()
#3. Run deep learning model

from sklearn.model_selection import train_test_split
from mcfly import modelgen, find_architecture, storage
from keras.models import load_model
import os

#use just 1 channel for now, later split into 2 and see if it makes a difference
X_train, X_test, y_train_list, y_test_list = train_test_split(matrix, labels, test_size=.4, random_state=42)

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

