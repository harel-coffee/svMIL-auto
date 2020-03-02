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

svType = 'DEL'

#Get the bags
with open(sys.argv[1], 'rb') as handle:
	bagDict = pkl.load(handle)

degPairs = np.loadtxt(sys.argv[2], dtype='object')


#normalize
currentMax = [0]*30
currentMin = [float('inf')]*30
normalizedBagDict = dict()
unnormalizedInstances = []
for pair in bagDict:
	
	for instance in bagDict[pair]:
		
		if instance[0] == '0' and instance[1] == '0':
			continue
		
		for featureInd in range(0, len(instance)):
			feature = instance[featureInd]
			if feature > currentMax[featureInd]:
				currentMax[featureInd] = feature
			if feature < currentMin[featureInd]:
				currentMin[featureInd] = feature

typeValues = []
for pair in bagDict:
	
	normalizedBagDict[pair] = []
	
	for instance in bagDict[pair]:
		
		normInstance = []
		for featureInd in range(0, len(instance)):
			
			feature = instance[featureInd]
			
			if currentMin[featureInd] == 0 and currentMax[featureInd] == 0:
				normInstance.append(0)
				continue 
			
			normFeature = (feature-currentMin[featureInd])/(currentMax[featureInd]-currentMin[featureInd])
			normInstance.append(normFeature)
			
			if featureInd == 28:
				if feature not in typeValues:
					typeValues.append(feature)
		
			
		normalizedBagDict[pair].append(normInstance)
print(typeValues)
bagDict = normalizedBagDict

#determine the bag labels given a file of DEG pairs
#degPairs = np.load(sys.argv[2], allow_pickle=True, encoding='latin1')

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
	instanceLabels = []
	posCount = 0
	negCount = 0
	positiveBags = []
	negativeBags = []
	removedPathwayPairs = 0
	
	for pair in bagDict:
		
		splitPair = pair.split("_")
		shortPair = splitPair[7] + '_' + splitPair[0]
		
		if svType != '':	
			if splitPair[12] != svType:
				continue
		
		#get the label of the bag by checking if it exists in degPairs
		if shortPair in degPairs[:,0]: #only add this gene if it has a z-score, so is not something with mutations. 
			
			#check if this is true or not.
			degPairInfo = degPairs[degPairs[:,0] == shortPair][0]
			
			#if degPairInfo[3] == 'True':
			if float(degPairInfo[5]) > 2 or float(degPairInfo[5]) < -2:	
				#if pair in pathwayAnnotation[:,0]: #skip the ones that have possible pathway effects
				#	removedPathwayPairs += 1
					#continue
				bagLabels.append(1)
				posCount += 1
				#get the right number of features per instance
				instances = []
				for instance in bagDict[pair]:
					
					#instances.append(instance[0:featureInd+1])
					if instance[0] == 0 and instance[1] == 0:
						continue
					
					instances.append(instance[0:featureInd])
					
					instanceLabels.append(pair)
					
				if len(instances) < 1:
					continue
				
				positiveBags.append(instances)
				
			else:
				
				bagLabels.append(0)
				negCount += 1
				
				#get the right number of features per instance
				instances = []
				for instance in bagDict[pair]:
					if instance[0] == 0 and instance[1] == 0:
						continue
					#instances.append(instance[0:featureInd+1])
					instances.append(instance[0:featureInd])
					instanceLabels.append(pair)
					
				#print(instances)
				if len(instances) < 1:
					continue
			
				negativeBags.append(instances)
				
				#negativeBags.append(bagDict[pair])
	
		#bags.append(bagDict[pair])
	
	positiveBags = np.array(positiveBags)
	negativeBags = np.array(negativeBags)
	
	print(positiveBags.shape)
	print(negativeBags.shape)
	#take a random subset for speed
	#positiveBags = np.random.choice(positiveBags, 2000)
	#negativeBags = np.random.choice(negativeBags, 2000)
	
	np.random.seed(0)
	negativeBagsSubsampled = np.random.choice(negativeBags, posCount)
	
	bags = np.array(bags)
	bags = np.concatenate((positiveBags, negativeBagsSubsampled))
	bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBagsSubsampled.shape[0])
	#bags = np.concatenate((positiveBags, negativeBags))
	#bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBags.shape[0])
	
	instances = np.vstack(bags)
	
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

# 
# pca = PCA(n_components=2)
# projected = pca.fit_transform(similarityMatrices[featureCount])
# projectedWithOffset = projected
# 
# for row in range(0, projected.shape[0]):
# 	for col in range(0, projected.shape[1]):
# 		projectedWithOffset[row][col] += np.random.normal(-1, 1) * 0.5
# 		
# projected = projectedWithOffset
# 
# colorLabels = []
# 
# for label in bagLabels:
# 	
# 	if label == 1:
# 		colorLabels.append('r')
# 	else:
# 		colorLabels.append('b')
# 
# fig,ax=plt.subplots(figsize=(7,5))
# plt.scatter(projected[:, 0], projected[:, 1], edgecolors=colorLabels, facecolors='none')
# #plt.show()


# #rasterize the PCA plot and make a density heatmap
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
# boxWidth = 2
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

from random import shuffle
#shuffle(bagLabels)

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold
from inspect import signature
from sklearn import model_selection
from sklearn.metrics import average_precision_score
from sklearn.metrics import plot_roc_curve
from sklearn.metrics import make_scorer, accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve

#f, axes = plt.subplots()	
def cvClassification(similarityMatrix, bagLabels, clf, color, plot):
	
	shuffle(bagLabels)
	
	scoring = {'accuracy' : make_scorer(accuracy_score), 
			   'precision' : make_scorer(precision_score),
			   'recall' : make_scorer(recall_score), 
			   'f1_score' : make_scorer(f1_score),
			   'average_precision' : make_scorer(average_precision_score),
			   'auc' : make_scorer(roc_auc_score)}
	
	kfold = model_selection.StratifiedKFold(n_splits=10, shuffle=True, random_state=10)
	
	results = model_selection.cross_validate(estimator=clf,
											  X=similarityMatrix,
											  y=bagLabels,
											  cv=kfold,
											  scoring=scoring,
											  return_estimator=True)

	
	print('accuracy: ', np.mean(results['test_accuracy']), np.std(results['test_accuracy']))
	#print('precision: ', np.mean(results['test_precision']), np.std(results['test_precision']))
	#print('recall: ', np.mean(results['test_recall']), np.std(results['test_recall']))
	print('F1 score: ', np.mean(results['test_f1_score']), np.std(results['test_f1_score']))
	#print('AP: ', np.mean(results['test_average_precision']), np.std(results['test_average_precision']))
	print('AUC: ', np.mean(results['test_auc']), np.std(results['test_auc']))
	
	print(results['test_auc'])
		
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

	
	# instanceImportances = dict()
	# for idx,estimator in enumerate(results['estimator']):
	# 	print("Features sorted by their score for estimator {}:".format(idx))
	# 	feature_importances = estimator.feature_importances_
	# 	print(feature_importances)
	# 	
	# 	for i in range(0, len(feature_importances)):
	# 		if i not in instanceImportances:
	# 			instanceImportances[i] = []
	# 		instanceImportances[i].append(feature_importances[i])
	# 
	# scores = []	
	# for instance in instanceImportances:
	# 	scores.append([instance, np.mean(instanceImportances[instance]), np.std(instanceImportances[instance])])
	# 
	# scores = np.array(scores)
	# scores = scores[scores[:,1].argsort()][::-1]
	# 
	# print(scores)
	# 
	# #plt.bar(np.arange(scores.shape[0]), scores[:,1])
	# #plt.show()
	# 
	# for i in range(0,10):
	# 	score = scores[i]
	# 	
	# 	print(i)
	# 	
	# 	print(instanceLabels[int(score[0])])
	# 	print(instances[int(score[0])])
	# 	
	exit()
	
	return np.mean(results['test_average_precision']), np.mean(results['test_f1_score'])
	
rfClassifier = RandomForestClassifier(max_depth=100, n_estimators=2)

from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import train_test_split

#Number of trees in random forest
n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
# Number of features to consider at every split
max_features = ['auto', 'sqrt']
# Maximum number of levels in tree
max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
max_depth.append(None)
# Minimum number of samples required to split a node
min_samples_split = [2, 5, 10]
# Minimum number of samples required at each leaf node
min_samples_leaf = [1, 2, 4]
# Method of selecting samples for training each tree
bootstrap = [True, False]# Create the random grid
random_grid = {'n_estimators': n_estimators,
               'max_features': max_features,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf,
               'bootstrap': bootstrap}
# rf_random = RandomizedSearchCV(estimator = rfClassifier, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)
# 
# X_train, X_test, y_train, y_test = train_test_split(similarityMatrices[featureStart], bagLabels, test_size=0.33, random_state=42)
# rf_random.fit(X_train, y_train)
# print(rf_random.best_params_)
# print(rf_random.score(X_test, y_test))
# 
# print('base:')
# rfClassifier.fit(X_train, y_train)
# print(rfClassifier.score(X_test, y_test))
# exit()
#best DEL: n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True
#best ITX: n_estimators= 1000, min_samples_split=2, min_samples_leaf=1, max_features='sqrt', max_depth=110, bootstrap=True
#best inv: n_estimators= 400, min_samples_split=10, min_samples_leaf=2, max_features='auto', max_depth=40, bootstrap=False
#best DUP: n_estimators= 1200, min_samples_split=2, min_samples_leaf=4, max_features='sqrt', max_depth=70, bootstrap=False

from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

names = [ "Random Forest"]

classifiers = [
	#KNeighborsClassifier(3),
	#SVC(kernel="linear", C=0.025),
	#SVC(gamma=2, C=1),
	#GaussianProcessClassifier(1.0 * RBF(1.0)),
	#DecisionTreeClassifier(max_depth=100),
	#RandomForestClassifier(n_estimators=200),
	RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)]
	#MLPClassifier(alpha=1, max_iter=1000),
	#AdaBoostClassifier(),
	#GaussianNB(),
	#QuadraticDiscriminantAnalysis()]

for classifierInd in range(0, len(classifiers)):
	print(names[classifierInd])
	rfClassifier = classifiers[classifierInd]
	#For each feature round, save the f1 and AP scores
	#Then at the end, show the curves only for the last round (all features)
	#show a separate figure for the feature selection scores
	aps = []
	f1s = []
	for featureInd in range(featureStart, featureCount+1):
		
		plot = False
		if featureInd == featureCount:
			plot = True
	
		ap, f1 = cvClassification(similarityMatrices[featureInd], bagLabels, rfClassifier, 'black', plot)
		aps.append(ap)
		f1s.append(f1)
	


for classifierInd in range(0, len(classifiers)):
	print(names[classifierInd])
	rfClassifier = classifiers[classifierInd]
	#repeat for shuffled labels, do this separately to make sure that labels are only shuffled from here
	shuffle(bagLabels)
	shuffledAps = []
	shuffledF1s = []
	print("shuffled")
	for featureInd in range(featureStart, featureCount+1):
	
		plot = False
		if featureInd == featureCount:
			plot = True
		shuffledAp, shuffledF1 = cvClassification(similarityMatrices[featureInd], bagLabels, rfClassifier, 'red', plot)
		shuffledAps.append(shuffledAp)
		shuffledF1s.append(shuffledF1)
	
exit()