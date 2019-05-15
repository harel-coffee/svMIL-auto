import sys
import numpy as np

from sklearn import metrics
from sklearn.metrics import roc_curve, auc

## Running MILES with somatic only, using the > 3 patient DEG labels

somaticScores = np.loadtxt(sys.argv[1], dtype="object")
degLabels = np.loadtxt(sys.argv[2], dtype="object")

#Make bags

#Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted

#First assign the feature vectors to the right SV in order

svBagContents = dict()
genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
for geneSVPair in somaticScores:
	
	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
	score = False
	for feature in features:
		if feature > 0:
			score = True
			break
	if score == True:	
		splitGeneSVPairName = geneSVPair[0].split("_")
		geneName = splitGeneSVPairName[0]
		
		#The first element will be the gene name, the rest is the SV information
		splitGeneSVPairName.pop(0) #remove first element
		sv = "_".join(splitGeneSVPairName)
		
		if sv not in genesPerBag:
			genesPerBag[sv] = []
		genesPerBag[sv].append(geneName)

		if sv not in svBagContents:
			svBagContents[sv] = []
		
		svBagContents[sv].append(features)
	

#Then construct the bags and add the right labels

bags = []
labels = []
pairNames = []
for sv in svBagContents:
	
	bagLabel = 0
	for gene in genesPerBag[sv]:	
		pairNames.append(sv + "_" + gene)
		
		#The bag can only be negative if all instances are negative. 
		
		#get the patient ID
		splitSVName = sv.split("_")
		patientId = splitSVName[len(splitSVName)-1]
		
		if gene in degLabels[:,0]: #we already make bags for sv-gene pairs, so if the gene is not affected by an SV, the bag will not be positive. 
			bagLabel = 1

		
	bags.append(svBagContents[sv])
	labels.append(bagLabel)
	
# # Running MILES with somatic only, but using per patient-gene pair DEG labels
# 
# somaticScores = np.loadtxt(sys.argv[1], dtype="object")
# degLabels = np.loadtxt(sys.argv[2], dtype="object")
# 
# 
# #Make bags
# 
# #Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted
# 
# #First assign the feature vectors to the right SV in order
# 
# svBagContents = dict()
# genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
# for geneSVPair in somaticScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 	
# 
# #Then construct the bags and add the right labels
# 
# bags = []
# labels = []
# pairNames = []
# for sv in svBagContents:
# 	
# 	bagLabel = 0
# 	for gene in genesPerBag[sv]:	
# 		pairNames.append(sv + "_" + gene)
# 		
# 		#The bag can only be negative if all instances are negative. 
# 		
# 		#get the patient ID
# 		splitSVName = sv.split("_")
# 		patientId = splitSVName[len(splitSVName)-1]
# 		
# 		if gene + "_" + patientId in degLabels[:,0]:
# 			bagLabel = 1
# 
# 		
# 	bags.append(svBagContents[sv])
# 	labels.append(bagLabel)
# 	
# # Constructing MILES for somatic vs germline bags
# 
# somaticScores = np.loadtxt(sys.argv[1], dtype="object")
# germlineScores = np.loadtxt(sys.argv[2], dtype="object")
# 
# #Make bags
# 
# #Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted
# #A bag is positive if it comes from the somatic set, and negative if it comes from the germline set
# 
# #First assign the feature vectors to the right SV in order
# 
# svBagContents = dict()
# genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
# for geneSVPair in somaticScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		# samples = somaticAllSampleScores[somaticAllSampleScores[:,0] == geneName, 31][0]
# 		# splitSamples = samples.split(",")
# 		# if len(splitSamples) < 7:
# 		# 	continue
# 		
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 	
# 
# #Then construct the bags
# 
# bags = []
# labels = []
# pairNames = []
# for sv in svBagContents:
# 	
# 	for gene in genesPerBag[sv]:	
# 		pairNames.append(sv + "_" + gene)
# 	bags.append(svBagContents[sv])
# 	labels.append(1)
# 	
# svBagContents = dict()
# genesPerBag = dict()
# for geneSVPair in germlineScores:
# 	
# 	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
# 	score = False
# 	for feature in features:
# 		if feature > 0:
# 			score = True
# 			break
# 	if score == True:	
# 		splitGeneSVPairName = geneSVPair[0].split("_")
# 		geneName = splitGeneSVPairName[0]
# 		
# 		# samples = germlineAllSampleScores[germlineAllSampleScores[:,0] == geneName, 31][0]
# 		# splitSamples = samples.split(",")
# 		# if len(splitSamples) < 7:
# 		# 	continue
# 		
# 		#The first element will be the gene name, the rest is the SV information
# 		splitGeneSVPairName.pop(0) #remove first element
# 		sv = "_".join(splitGeneSVPairName)
# 		
# 		if sv not in genesPerBag:
# 			genesPerBag[sv] = []
# 		genesPerBag[sv].append(geneName)
# 		
# 		if sv not in svBagContents:
# 			svBagContents[sv] = []
# 		
# 		svBagContents[sv].append(features)
# 
# 
# #Then construct the bags
# 
# #bags = []
# for sv in svBagContents:
# 	for gene in genesPerBag[sv]:	
# 		pairNames.append(sv + "_" + gene)
# 	bags.append(svBagContents[sv])
# 	labels.append(0)

bags = np.array(bags)
bags = bags[1:500]
labels = labels[1:500]


from sklearn.model_selection import train_test_split
trainBags, testBags, trainLabels, testLabels = train_test_split(bags, labels, test_size=.4, random_state=42) #60%/40%


#MILES classifier


#3. Generate a similarity matrix

trainInstances = np.vstack(trainBags)
print "Number of training instances: ", trainInstances.shape[0]

print "generating similarity matrix"
similarityMatrix = np.zeros([len(trainBags), trainInstances.shape[0]])

print "size of the similarity matrix: "
print similarityMatrix.shape

print "size of the training bags: "
print trainBags.shape


#Instead of the distance from the center, compute the distance from every instance to every instance in the other bags, and take the smallest distance as the similarity

#1. First compute the instance x instance distance matrix
instanceDistMatrix = np.zeros([trainInstances.shape[0], trainInstances.shape[0]])

from scipy.spatial.distance import pdist, squareform

#instanceDistMatrix = squareform(pdist(trainInstances))


#2. Make a map where we know which instances belong to which bag
bagMap = dict()
reverseBagMap = dict()
geneInd = 0
for bagInd in range(0, trainBags.shape[0]):
	reverseBagMap[bagInd] = []
	for gene in trainBags[bagInd]:
		bagMap[geneInd] = bagInd
		reverseBagMap[bagInd].append(geneInd)
		
		geneInd += 1

#3. Go through the bags, and find which instance in this bag has the smallest distance to all other instances (make sure to exclude itself)
#4. Make a bag x instance similarity matrix
bagInstanceSimilarityTrain = np.zeros([trainBags.shape[0], trainInstances.shape[0]])
print "Number of training bags: ", trainBags.shape[0]
for bagInd in range(0, trainBags.shape[0]):
	
	print bagInd
	
	#Get the indices of the instances that are in this bag
	instanceIndices = reverseBagMap[bagInd]
	
	trainInstanceSubset = trainInstances[instanceIndices,:]
	
	#Compute the pairwise distance matrix here
	minDistance = float("inf")
	minDistanceInd = 0
	for instanceInd in range(0, trainInstanceSubset.shape[0]):
		instance = trainInstanceSubset[instanceInd]
		distance = np.abs(instance - trainInstances)

		summedDistance = np.sum(distance,axis=1)

		currentMinDistance = np.min(summedDistance)
		if currentMinDistance < np.min(minDistance):
			minDistance = summedDistance
			minDistanceInd = instanceInd
		
		# currentMinDistance = np.median(summedDistance)
		# if currentMinDistance < np.median(minDistance):
		# 	minDistance = summedDistance
		# 	minDistanceInd = instanceInd

	#This instance will be used as representative for this bag. We use this value as the similarity to all other instances.  
	bagInstanceSimilarityTrain[bagInd] = minDistance

	
#4. Train a classifier on the similarity space
from sklearn.ensemble import RandomForestClassifier
print "training the classifier in similarity space"
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
rfClassifier.fit(bagInstanceSimilarityTrain, trainLabels) #Use the bag labels, not the instance labels



#5. Repeat the same steps for the test bags
bagMap = dict()
reverseBagMap = dict()
geneInd = 0
for bagInd in range(0, testBags.shape[0]):
	reverseBagMap[bagInd] = []
	for gene in testBags[bagInd]:
		bagMap[geneInd] = bagInd
		reverseBagMap[bagInd].append(geneInd)
		
		geneInd += 1

testInstances = np.vstack(testBags)

#3. Go through the bags, and find which instance in this bag has the smallest distance to all other instances (make sure to exclude itself)
#4. Make a bag x instance similarity matrix
bagInstanceSimilarityTest = np.zeros([testBags.shape[0], trainInstances.shape[0]])
print "Number of test bags: ", testBags.shape[0]
for bagInd in range(0, testBags.shape[0]):
	
	print bagInd
	
	#Get the indices of the instances that are in this bag
	instanceIndices = reverseBagMap[bagInd]
	
	testInstanceSubset = testInstances[instanceIndices,:]
	
	#Compute the pairwise distance matrix here
	minDistance = float("inf")
	minDistanceInd = 0
	for instanceInd in range(0, testInstanceSubset.shape[0]):
		instance = testInstanceSubset[instanceInd]
		distance = np.abs(instance - trainInstances) #compute the distances to the train instances, otherwise we are not in the same similarity space. 

		summedDistance = np.sum(distance,axis=1)
		
		
		currentMinDistance = np.min(summedDistance)
		if currentMinDistance < np.min(minDistance):
			minDistance = summedDistance
			minDistanceInd = instanceInd
		
		# currentMinDistance = np.median(summedDistance)
		# if currentMinDistance < np.median(minDistance):
		# 	minDistance = summedDistance
		# 	minDistanceInd = instanceInd

	#This instance will be used as representative for this bag. We use this value as the similarity to all other instances.  
	bagInstanceSimilarityTest[bagInd] = minDistance

#6. Predict bag labels and instance labels
# 
# print "scoring classifier"
# score = rfClassifier.score(bagInstanceSimilarityTest, testLabels) #First try with classifying bags, later we can do instance-based classification.
# print score
# 
# predictions = rfClassifier.predict(bagInstanceSimilarityTest)
# print np.average(testLabels == np.sign(predictions))
# # positiveInd = np.where(np.array(predictions) == 1)[0]
# # print positiveInd
# 
# predictionProbs = rfClassifier.predict_proba(bagInstanceSimilarityTest)[:,1]
# 
# fpr, tpr, thresholds = metrics.roc_curve(testLabels, predictionProbs, pos_label=1.)
# roc_auc = auc(fpr, tpr)
# 
# print "AUC: ", roc_auc
# 
# print "SVM performance: "
# 
# from sklearn.svm import SVC
# clf = SVC(kernel='linear', probability=True)
# clf.fit(bagInstanceSimilarityTrain, trainLabels) 
# score = clf.score(bagInstanceSimilarityTest, testLabels)
# print score
# 
# predictionProbs = clf.predict_proba(bagInstanceSimilarityTest)[:,1]
# 
# fpr, tpr, thresholds = metrics.roc_curve(testLabels, predictionProbs, pos_label=1.)
# roc_auc = auc(fpr, tpr)
# 
# print "AUC: ", roc_auc
# print len(clf.coef_[0])
# print "Used features: ", len(clf.coef_[0] != 0)

print "lasso performance and AUC: "

from sklearn.linear_model import Lasso
from sklearn.metrics import auc, precision_recall_curve
# 
# alphas = [1e-3]
# 
# 
# for currentAlpha in alphas:
# 	print "alpha: ", currentAlpha
# 
# 	lasso = Lasso(alpha=currentAlpha)
# 	lasso.fit(bagInstanceSimilarityTrain,trainLabels)
# 	
# 	train_score=lasso.score(bagInstanceSimilarityTrain,trainLabels)
# 	test_score=lasso.score(bagInstanceSimilarityTest,testLabels)
# 	coeff_used = np.sum(lasso.coef_!=0)
# 	preds = lasso.predict(bagInstanceSimilarityTest)
# 	
# 	print "train score: ", train_score
# 	print "test score: ", test_score
# 	print "number of features used: ", coeff_used
# 
# 	precision, recall, thresholds = precision_recall_curve(testLabels, preds)
# 	aucScore = auc(recall, precision)
# 	print "AUC: ", aucScore
		 

	
	
	
	
# 	
# 	accuracy = np.average(labels[test].T == np.sign(predictions))
# 	accuracies.append(accuracy)
# 	
# 	fpr, tpr, thresholds = metrics.roc_curve(labels[test], predictionProbs, pos_label=1.)
# 	
# 	tprs.append(interp(mean_fpr, fpr, tpr))
# 	roc_auc = auc(fpr, tpr)
# 	aucs.append(roc_auc)


# 
# print "lasso performance: "
# #Train lasso instead of RF
# from sklearn.linear_model import Lasso
# alphas = [1e-15, 1e-10, 1e-8, 1e-4, 1e-3, 1e-2, 1, 5, 10, 20]
# 
# for currentAlpha in alphas:
# 	print "alpha: ", currentAlpha
# 
# 	lasso = Lasso(alpha=currentAlpha)
# 	lasso.fit(bagInstanceSimilarityTrain,trainLabels)
# 	
# 	train_score=lasso.score(bagInstanceSimilarityTrain,trainLabels)
# 	test_score=lasso.score(bagInstanceSimilarityTest,testLabels)
# 	coeff_used = np.sum(lasso.coef_!=0)
# 	
# 	print "train score: ", train_score
# 	print "test score: ", test_score
# 	print "number of features used: ", coeff_used


# #Get the probabilities for AUC
# coef = lasso.coef_
# data_loss = 0.5 * ((bagInstanceSimilarityTest.dot(coef) - testLabels) ** 2).sum()
# n_samples, n_features = bagInstanceSimilarityTest.shape
# penalty = n_samples * lasso.alpha * np.abs(coef).sum()
# likelihood = np.exp(-(data_loss + penalty))
# 
# fpr, tpr, thresholds = metrics.roc_curve(testLabels, likelihood, pos_label=1.)
# roc_auc = auc(fpr, tpr)
# 
# print "AUC: ", roc_auc

exit()


positivePairs = []
positiveGenes = dict()
for ind in positiveInd:
	positivePairs.append(pairNames[ind])
	splitPairName = pairNames[ind].split("_")
	geneName = splitPairName[len(splitPairName)-1]
	if geneName not in positiveGenes:
		positiveGenes[geneName] = 0
	positiveGenes[geneName] += 1

print positivePairs
print positiveGenes.keys()
print len(positiveGenes.keys())
exit()
