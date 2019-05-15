import sys
import numpy as np
import matplotlib.pyplot as plt


# Running MILES with somatic only, using the > 3 patient DEG labels
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
# 		if gene in degLabels[:,0]: #we already make bags for sv-gene pairs, so if the gene is not affected by an SV, the bag will not be positive. 
# 			bagLabel = 1
# 
# 		
# 	bags.append(svBagContents[sv])
# 	labels.append(bagLabel)


#Running MILES with somatic only, but using per patient-gene pair DEG labels


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


### Germline vs somatic miles
# 
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
# 

# bags = bags[1:100]
# # labels = labels[1:100]
# bags = np.array(bags)
# instances = np.vstack(bags)
# 
# 
# print "generating similarity matrix"
# 
# #Unfold the training bags so that we can compute the distance matrix at once to all genes
# bagMap = dict()
# reverseBagMap = dict()
# geneInd = 0
# for bagInd in range(0, bags.shape[0]):
# 	reverseBagMap[bagInd] = []
# 	for gene in bags[bagInd]:
# 		bagMap[geneInd] = bagInd
# 		reverseBagMap[bagInd].append(geneInd)
# 		
# 		geneInd += 1
# 
# 
# similarityMatrix = np.zeros([bags.shape[0], instances.shape[0]])
# print "Number of bags: ", bags.shape[0]
# for bagInd in range(0, bags.shape[0]):
# 	
# 	print bagInd
# 	
# 	#Get the indices of the instances that are in this bag
# 	instanceIndices = reverseBagMap[bagInd]
# 	
# 	instanceSubset = instances[instanceIndices,:]
# 	
# 	#Compute the pairwise distance matrix here
# 	minDistance = float("inf")
# 	minDistanceInd = 0
# 	for instanceInd in range(0, instanceSubset.shape[0]):
# 		instance = instanceSubset[instanceInd]
# 		distance = np.abs(instance - instances) #compute the distances to the train instances, otherwise we are not in the same similarity space. 
# 
# 		summedDistance = np.sum(distance,axis=1)
# 
# 		currentMinDistance = np.min(summedDistance)
# 		if currentMinDistance < np.min(minDistance):
# 			minDistance = summedDistance
# 			minDistanceInd = instanceInd
# 
# 	#This instance will be used as representative for this bag. We use this value as the similarity to all other instances.  
# 	similarityMatrix[bagInd] = minDistance

# #4. Train a classifier on the similarity space
# from sklearn.ensemble import RandomForestClassifier
# print "training the classifier in similarity space"
# np.random.seed(500)
# rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
# rfClassifier.fit(similarityMatrix, labels) #Use the bag labels, not the instance labels
# 
# print similarityMatrix
# print labels
# 
# predictions = rfClassifier.predict(similarityMatrix)
# print predictions
# print np.average(labels == np.sign(predictions))
# 
# #Now select the most important features with random forest
# importances = rfClassifier.feature_importances_
# std = np.std([tree.feature_importances_ for tree in rfClassifier.estimators_],
#              axis=0)
# indices = np.argsort(importances)[::-1]
# 
# 
# nonZeroIndices = []
# for index in indices:
# 	if importances[index] > 0:
# 		nonZeroIndices.append(index)
# 
# #Get the genes at these indices
# positiveGenes = dict()
# positivePairs = []
# for index in nonZeroIndices:
# 	pairName = pairNames[index]
# 	positivePairs.append(pairName)
# 	splitPairName = pairName.split("_")
# 	geneName = splitPairName[len(splitPairName)-1]
# 	if geneName not in positiveGenes:
# 		positiveGenes[geneName] = 0
# 	positiveGenes[geneName] += 1

### Loading pre-made data to save time

#To save time, bags and labels have been stored on disk already and can be re-loaded
bags = np.load("3DEGs/bags.txt.npy")
labels = np.load("3DEGs/labels.txt.npy")
pairNames = np.load("3DEGs/pairNames.txt.npy") #the sv-gene pair names of each bag entry
similarityMatrix = np.load("3DEGs/similarityMatrix.txt.npy")

#Shuffle the labels
np.random.shuffle(labels)


# similarityMatrix = similarityMatrix[1:100, 1:100]
# labels = labels[1:100]

### Using SVM
# from sklearn.metrics import auc, precision_recall_curve
# print "SVM performance: "
# 
# from sklearn.svm import LinearSVC
# clf = LinearSVC()
# print "Fitting classifier: "
# clf.fit(similarityMatrix, labels)
# print "Scoring classifier: "
# score = clf.score(similarityMatrix, labels)
# print score
# 
# #predictionProbs = clf.predict_proba(similarityMatrix)[:,1]
# preds = clf.predict(similarityMatrix)
# predsDiff = np.average(labels == np.sign(preds))
# print "mean score: ", predsDiff
# 
# precision, recall, thresholds = precision_recall_curve(labels, preds)
# aucScore = auc(recall, precision)
# 
# print "AUC: ", aucScore
# print "Used features: ", len(clf.coef_[0] != 0)
# 
# 
# #collect the pairs in this order
# pairsRanking = np.empty([len(clf.coef_[0]), 2], dtype="object")
# for coefInd in range(0, len(clf.coef_[0])):
# 	pairName = pairNames[coefInd]
# 	pairsRanking[coefInd, 0] = pairName
# 	pairsRanking[coefInd, 1] = clf.coef_[0][coefInd]
# 
# #Sort the features by importance/weights
# pairsRanking = pairsRanking[np.argsort(np.abs(pairsRanking[:,1]))[::-1]]	
# 	
# #Output the pairs in this order
# pairsRankingOut = "svmSomaticGermline/pairsRanking.txt"
# with open(pairsRankingOut, 'w') as outF:
# 	for pairInd in range(0, pairsRanking.shape[0]):
# 		outF.write(pairsRanking[pairInd, 0] + "\t" + str(pairsRanking[pairInd, 1]) + "\n")

# ###Using SVM with cross validation
# print "SVM with RFECV performance: "
# from sklearn.feature_selection import RFECV
# from sklearn.metrics import auc, precision_recall_curve
# from sklearn.svm import LinearSVC
# clf = LinearSVC()
# selector = RFECV(clf, step=1, cv=5)
# print "Fitting classifier: "
# selector = selector.fit(similarityMatrix, labels)
# print "score: ", selector.score(similarityMatrix, labels)
# preds = selector.predict(similarityMatrix)
# predsDiff = np.average(labels == np.sign(preds))
# 
# print "mean score: ", predsDiff
# 
# precision, recall, thresholds = precision_recall_curve(labels, preds)
# aucScore = auc(recall, precision)
# 
# print "AUC: ", aucScore
# 
# selectedGenesInd = selector.support_
# selectedPairs = pairNames[selectedGenesInd]
# 
# print "Number of selected pairs: ", selectedPairs
# pairsRankingOut = "svmRFECVPerPatient/pairsRanking.txt"
# with open(pairsRankingOut, 'w') as outF:
# 	for pairInd in range(0, selectedPairs.shape[0]):
# 		outF.write(selectedPairs[pairInd] + "\n")

#output the similarity matrix and also the labels and bag indices so that we can do analysis after running once
#np.save("svmRFECVPerPatient/bags.txt", bags)
#np.save("svmRFECVPerPatient/pairNames.txt", pairNames) #the sv-gene pair names of each bag entry
#np.save("svmRFECVPerPatient/labels.txt", labels)
#np.save("svmSomaticGermline/similarityMatrix.txt", similarityMatrix) #no need to store this every time
#np.save("svmRFECVPerPatient/coef.txt", clf.coef_[0])



### Using Lasso
from sklearn.linear_model import Lasso
from sklearn.metrics import auc, precision_recall_curve
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits=10)
labels = np.array(labels)
alphas = [1e-15, 1e-10, 1e-8, 1e-4, 1e-3, 1e-2, 1, 5, 10, 20]
#alphas = [1e-4, 1e-3, 1e-2, 1]
#alphas = [1e-2]

accs = dict()
aucs = dict()
coeffs = dict()
predDiffs = dict()

# # Lasso but then with CV
# for currentAlpha in alphas:
# 	print "alpha: ", currentAlpha
# 	accs[currentAlpha] = []
# 	aucs[currentAlpha] = []
# 	coeffs[currentAlpha] = []
# 	predDiffs[currentAlpha] = []
# 	for train, test in cv.split(similarityMatrix, labels):
# 		
# 		lasso = Lasso(alpha=currentAlpha)
# 		lasso.fit(similarityMatrix[train],labels[train])
# 		
# 		#train_score=lasso.score(bagInstanceSimilarityTrain,trainLabels)
# 		test_score=lasso.score(similarityMatrix[test],labels[test])
# 		coeff_used = np.sum(lasso.coef_!=0)
# 		preds = lasso.predict(similarityMatrix[test])
# 		predsDiff = np.average(labels[test] == np.sign(preds))
# 		
# 		precision, recall, thresholds = precision_recall_curve(labels[test], preds)
# 		aucScore = auc(recall, precision)
# 		
# 		accs[currentAlpha].append(test_score)
# 		aucs[currentAlpha].append(aucScore)
# 		coeffs[currentAlpha].append(coeff_used)
# 		predDiffs[currentAlpha].append(predsDiff)
# 
# 	#Report the averages per alpha
# 	
# 	print "Actual acc: ", np.mean(predDiffs[currentAlpha])
# 	print "Mean acc: ", np.mean(accs[currentAlpha])
# 	print "Mean AUC: ", np.mean(aucs[currentAlpha])
# 	print "Mean coeffs: ", np.mean(coeffs[currentAlpha])
# 	
# 	np.save("lassoPerPatient/acc.txt", np.mean(predDiffs[currentAlpha]))
# 	np.save("lassoPerPatient/preds.txt", np.mean(accs[currentAlpha]))
# 	np.save("lassoPerPatient/auc.txt", np.mean(aucs[currentAlpha]))
# 	np.save("lassoPerPatient/coeffs.txt", np.mean(coeffs[currentAlpha]))
# 	
# 	

# #####Getting the concept genes without cross validation
currentAlpha = 1e-2
lasso = Lasso(alpha=currentAlpha)
lasso.fit(similarityMatrix,labels)

#train_score=lasso.score(bagInstanceSimilarityTrain,trainLabels)
test_score=lasso.score(similarityMatrix,labels)
coeff_used = np.sum(lasso.coef_!=0)
preds = lasso.predict(similarityMatrix)
predsDiff = np.average(labels == np.sign(preds))

precision, recall, thresholds = precision_recall_curve(labels, preds)
aucScore = auc(recall, precision)

print "lasso 2 patient alpha 0.01"
print "acc: ", test_score
print "predsDiff: ", predsDiff
print "auc: ", aucScore
print "coeffs: ", coeff_used

geneIndices = np.where(lasso.coef_ !=0)[0]
positiveGenes = dict()
positivePairs = []
for index in geneIndices:
	pairName = pairNames[index]
	positivePairs.append(pairName)
	splitPairName = pairName.split("_")
	geneName = splitPairName[len(splitPairName)-1]
	if geneName not in positiveGenes:
		positiveGenes[geneName] = 0
	positiveGenes[geneName] += 1

print len(positivePairs)
print len(positiveGenes)

milesConceptGenesOut = "lasso2Patients/milesConceptGenes.txt"
with open(milesConceptGenesOut, 'w') as outF:
	for gene in positiveGenes:
		outF.write(gene + "\t" + str(positiveGenes[gene]) + "\n")


#output the similarity matrix and also the labels and bag indices so that we can do analysis after running once
# np.save("lasso2Patients/alpha.txt", currentAlpha)
# np.save("lasso2Patients/bags.txt", bags)
# np.save("lasso2Patients/pairNames.txt", pairNames) #the sv-gene pair names of each bag entry
# np.save("lasso2Patients/labels.txt", labels)
# np.save("lasso2Patients/similarityMatrix.txt", similarityMatrix)
# np.save("lasso2Patients/conceptIndices.txt", geneIndices)

#Instead save the scores to file
# np.save("lassoPerPatient/acc_random.txt", test_score)
# np.save("lassoPerPatient/preds_random.txt", predsDiff)
# np.save("lassoPerPatient/auc_random.txt", aucScore)
# np.save("lassoPerPatient/coeffs_random.txt", coeff_used)

exit()



