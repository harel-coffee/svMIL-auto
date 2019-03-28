import sys
import numpy as np


somaticScores = np.loadtxt(sys.argv[1], dtype="object")
germlineScores = np.loadtxt(sys.argv[2], dtype="object")

#Use these to filter the gene-SV pairs by the number of smaples affecting them in total
somaticAllSampleScores = np.loadtxt(sys.argv[3], dtype="object")
germlineAllSampleScores = np.loadtxt(sys.argv[4], dtype="object")

#Make bags

#Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted
#A bag is positive if it comes from the somatic set, and negative if it comes from the germline set

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
		
		# samples = somaticAllSampleScores[somaticAllSampleScores[:,0] == geneName, 31][0]
		# splitSamples = samples.split(",")
		# if len(splitSamples) < 7:
		# 	continue
		
		
		#The first element will be the gene name, the rest is the SV information
		splitGeneSVPairName.pop(0) #remove first element
		sv = "_".join(splitGeneSVPairName)
		
		if sv not in genesPerBag:
			genesPerBag[sv] = []
		genesPerBag[sv].append(geneName)

		if sv not in svBagContents:
			svBagContents[sv] = []
		
		svBagContents[sv].append(features)
	

#Then construct the bags

bags = []
labels = []
pairNames = []
for sv in svBagContents:
	
	for gene in genesPerBag[sv]:	
		pairNames.append(sv + "_" + gene)
	bags.append(svBagContents[sv])
	labels.append(1)
	
svBagContents = dict()
genesPerBag = dict()
for geneSVPair in germlineScores:
	
	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
	score = False
	for feature in features:
		if feature > 0:
			score = True
			break
	if score == True:	
		splitGeneSVPairName = geneSVPair[0].split("_")
		geneName = splitGeneSVPairName[0]
		
		# samples = germlineAllSampleScores[germlineAllSampleScores[:,0] == geneName, 31][0]
		# splitSamples = samples.split(",")
		# if len(splitSamples) < 7:
		# 	continue
		
		#The first element will be the gene name, the rest is the SV information
		splitGeneSVPairName.pop(0) #remove first element
		sv = "_".join(splitGeneSVPairName)
		
		if sv not in genesPerBag:
			genesPerBag[sv] = []
		genesPerBag[sv].append(geneName)
		
		if sv not in svBagContents:
			svBagContents[sv] = []
		
		svBagContents[sv].append(features)


#Then construct the bags

#bags = []
for sv in svBagContents:
	for gene in genesPerBag[sv]:	
		pairNames.append(sv + "_" + gene)
	bags.append(svBagContents[sv])
	labels.append(-1)

bags = np.array(bags)
from sklearn.model_selection import train_test_split
trainBags, testBags, trainLabels, testLabels = train_test_split(bags, labels, test_size=.4, random_state=42) #60%/40%


#MILES classifier

#1. Obtain centers in each bag. Here I will use random points
print "computing centers"
centers = []
for bag in trainBags:
	
	#Number of instances in the bag
	bag = np.array(bag)
	instanceNum = bag.shape[0]
	randomInstanceInd = np.random.choice(range(0, instanceNum), 1)[0]
	
	randomInstance = bag[randomInstanceInd,:]
	centers.append(randomInstance)

print len(centers)


#2. Compute the distance from each cluster center to other points and take the minimum
def absoluteDistance(center, instances): #compute distance simply based on integer values
	#distances = np.abs(center - instances)
	#Then take the sum across all differences, and report one column
	#return np.sum(distances, axis=1)
	return 1

#3. Generate a similarity matrix

print "generating similarity matrix"
similarityMatrix = np.zeros([len(trainBags), len(centers)])

print "size of the similarity matrix: "
print similarityMatrix.shape

print "size of the training bags: "
print trainBags.shape

#The goal here is:
#- we want to have a matrix where the rows are the instances, and the columns are the bags.
#- To compute a distance from the bag to the instances, we want a bag center to all other instances, and take the smallest distance. 
#- Can we compute the distances quicker than now? At the moment, we first take a center for each bag, and then compute the distance to all other elements.
#- Currently, we have 350 bags. Then there will be 350x350 distances (bag to smallest instance in each bag).
#- Can we compute the distance to all instances at once and then take the minimum?

#Unfold the training bags so that we can compute the distance matrix at once to all genes
genes = []
bagMap = dict()
geneInd = 0
for bagInd in range(0, len(trainBags)):
	for gene in trainBags[bagInd]:
		genes.append(gene)
		bagMap[geneInd] = bagInd
		geneInd += 1

genes = np.array(genes)


for centerInd in range(0, len(centers)):
	print "center: ", centerInd

	distances = np.abs(centers[centerInd] - genes)
	summedDistances = np.sum(distances, axis=1)
	
	smallestDistanceInd = np.argmin(summedDistances)
	bagInd = bagMap[smallestDistanceInd]
	if bagInd == centerInd:
		continue
	
	#get the bag ind
	
	similarityMatrix[bagInd,centerInd] = summedDistances[smallestDistanceInd]
	
	
	# for bagInd in range(0, len(trainBags)):
	# 	
	# 	#Skip this bag if the center is in the current bag
	# 	if centerInd == bagInd:
	# 		continue
	# 	
	# 	#Compute the distance from center to all instances for each set of features with the same distance function at once
	# 	#allDistances = absoluteDistance(centers[centerInd], np.array(trainBags[bagInd]))
	# 	#distances = np.abs(centers[centerInd] - trainBags[bagInd])
	# 	distances = np.abs(centers[centerInd] - np.array([1]))
	# 	continue
	# 	#Here we should find the row (or the instance) that is closest to this center
	# 	
	# 	#Sum across the rows
	# 	#Identify the row (instance) that has the smallest total sum
	# 	summedDistances = np.sum(allDistances)
	# 
	# 	smallestDistance = np.min(summedDistances)
	# 	
	# 	similarityMatrix[bagInd, centerInd] = smallestDistance


#4. Train a classifier on the similarity space
from sklearn.ensemble import RandomForestClassifier
print "training the classifier in similarity space"
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
rfClassifier.fit(similarityMatrix, trainLabels) #Use the bag labels, not the instance labels

#5. Test the classifier on a new point (should this also be in the similarity space? Yes, right?)

#Convert test bags to a format that we can use
testInstances = np.vstack(testBags)

#Do we now compute the distance to the training data to get the same number of features?

testSimilarityMatrix = np.zeros([len(testBags), len(centers)])

#Computing the distance from every instance to all other instances is computationally expensive
#we should compute the distance from a random center
#But should the distance be computed to the center points of the training data? or to centers of the test data? training, right? otherwise the space is not the same.

#Here we should not use the distance from center to bags, but from center to each instance. Otherwise we cannot do instance-based classification!

print "computing test data similarity matrix"
genes = []
bagMap = dict()
geneInd = 0
for bagInd in range(0, len(testBags)):
	for gene in testBags[bagInd]:
		genes.append(gene)
		bagMap[geneInd] = bagInd
		geneInd += 1

genes = np.array(genes)

for centerInd in range(0, len(centers)):
	print "center: ", centerInd
	
	distances = np.abs(centers[centerInd] - genes)
	summedDistances = np.sum(distances, axis=1)
	smallestDistanceInd = np.argmin(summedDistances)
	#get the bag ind
	bagInd = bagMap[smallestDistanceInd]
	smallestDistance = summedDistances[smallestDistanceInd]
	
	if bagInd == centerInd:
		continue
	
	testSimilarityMatrix[bagInd, centerInd] = smallestDistance



print "similarity matrix shape: "
print testSimilarityMatrix.shape
print "test labels shape: "
print len(testLabels)
print np.where(np.array(testLabels) == 1)[0].shape
print np.where(np.array(testLabels) == -1)[0].shape
		
print "scoring classifier"
score = rfClassifier.score(testSimilarityMatrix, testLabels) #First try with classifying bags, later we can do instance-based classification.
print score

predictions = rfClassifier.predict(testSimilarityMatrix)
print predictions
print np.average(testLabels == np.sign(predictions))
positiveInd = np.where(np.array(predictions) == 1)[0]
print positiveInd

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
# 
# 
# #open the scores across samples, and write the scores to a new outfile subset
# totalScores = np.loadtxt(sys.argv[3], dtype='object')
# 
# subsetScores = np.empty([len(positiveGenes.keys()), totalScores.shape[1]], dtype="object")
# for geneInd in range(0, len(positiveGenes.keys())):
# 	
# 	geneName = positiveGenes.keys()[geneInd]
# 	
# 	totalScoresGene = totalScores[np.where(totalScores[:,0] == geneName)[0],:]
# 	subsetScores[geneInd,:] = totalScoresGene
# 	subsetScores[geneInd,2:31] = subsetScores[geneInd,2:31].astype(float)
# 	
# print subsetScores
# subsetScores = subsetScores[subsetScores[:,30].argsort()[::-1]] #Select the column  to rank by
# np.savetxt(sys.argv[4], subsetScores, delimiter='\t', fmt='%s')	
# 
# #Read the gene-patient pair DEG p-values
# 
# pairPValues = np.loadtxt(sys.argv[5], dtype="object")
# 
# #How many genes in this pvalue set are DEG?
# allSignGenes = []
# for pair in pairPValues:
# 	splitPairName = pair[0].split("_")
# 	geneName = splitPairName[0]
# 	if float(pair[1]) == 0 or float(pair[1]) <= 0.05:
# 		allSignGenes.append(geneName)
# 
# print np.unique(allSignGenes).shape
# #exit()
# 
# 
# print pairPValues
# signGenes = []
# for pair in positivePairs:
# 	splitPairName = pair.split("_")
# 	genePatientPairName = splitPairName[len(splitPairName)-1] + "_" + splitPairName[len(splitPairName)-2]
# 	
# 	pValue = pairPValues[np.where(pairPValues[:,0] == genePatientPairName)[0],1]
# 	if pValue.shape[0] > 0:
# 		if float(pValue[0]) == 0 or float(pValue[0]) <= 0.05:
# 			signGenes.append(splitPairName[len(splitPairName)-1])
# 			
# signGenes = np.unique(signGenes)	
# print len(signGenes)
# print len(signGenes) / float(len(positiveGenes.keys()))


#Simple classifier

print "Training the bag classifier: "

bags = np.array(bags)

#Shuffle the labels
import random
import os
from sklearn import metrics
sys.path.append(os.path.realpath('../../../../Methods/'))
from MILpy.Algorithms.simpleMIL import simpleMIL
from MILpy.functions.mil_cross_val import mil_cross_val
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from scipy import interp
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.utils.fixes import signature

#Shuffle 100 random times

tprs = []
aucs = []
totalAccuracy = 0
mean_fpr = np.linspace(0, 1, 100)

testReal = []
testProbs = []

for i in range(0,100):
	
	random.shuffle(labels)
	labels = np.array(labels)
	
	cv = StratifiedKFold(n_splits=10)
	SMILa = simpleMIL()
	
	iterationAucs = []
	iterationTprs = []
	accuracies = []
	testReal = []
	testProbs = []
	for train, test in cv.split(bags, labels):
		print len(train)
		print len(test)
		SMILa = simpleMIL() 
		SMILa.fit(bags[train], labels[train], type='average')
		predictions = SMILa.predict(bags[test])
		
		accuracy = np.average(labels[test].T == np.sign(predictions))
		accuracies.append(accuracy)
		
		fpr, tpr, thresholds = metrics.roc_curve(labels[test], predictions, pos_label=1.)
		
		iterationTprs.append(interp(mean_fpr, fpr, tpr))
		iterationTprs[-1][0] = 0.0
		roc_auc = auc(fpr, tpr)
		iterationAucs.append(roc_auc)
		
		# #precision-recall
		# precision, recall, thresholds = precision_recall_curve(labels[test], predictions)
		# print precision
		# print recall
		testReal.append(labels[test])
		testProbs.append(predictions)
	
		
	# precisions.append(precision)
	# recalls.append(recall)
	# 
	totalAccuracy += np.mean(accuracies)	
	tprs.append(np.mean(iterationTprs, axis=0))
	aucs.append(np.mean(iterationAucs))

#Plot precision-recall
testReal = np.concatenate(testReal)
testProbs = np.concatenate(testProbs)

precision, recall, _ = precision_recall_curve(testReal, testProbs)


# In matplotlib < 1.5, plt.fill_between does not have a 'step' argument
step_kwargs = ({'step': 'post'}
               if 'step' in signature(plt.fill_between).parameters
               else {})

plt.step(recall, precision, color='b', alpha=0.2,
         where='post')
plt.fill_between(recall, precision, alpha=0.2, color='b', **step_kwargs)

plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
#plt.title('2-class Precision-Recall curve')
plt.show()


plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
plt.plot(mean_fpr, mean_tpr, color='b',
         label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
         lw=2, alpha=.8)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                 label=r'$\pm$ 1 std. dev.')

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
title = 'Receiver operating characteristic, average accuracy: ' + str(totalAccuracy / float(2))
plt.title(title)
plt.legend(loc="lower right")
plt.show()

exit()	
# 	
# print "Average accuracy: ", totalAccuracy / float(100)
# print "Average AUC: ", totalAuc / float(100)

	# 
	# 
	# from sklearn.model_selection import train_test_split
	# trainBags, testBags, trainLabels, testLabels = train_test_split(bags, labels, test_size=.4, random_state=42) #60%/40%
	# 
	# import os
	# from sklearn import metrics
	# sys.path.append(os.path.realpath('../../../../Methods/'))
	# from MILpy.Algorithms.simpleMIL import simpleMIL
	# from MILpy.functions.mil_cross_val import mil_cross_val
	# 
	# testLabels = np.array(testLabels)
	# 
	# SMILa = simpleMIL() 
	# SMILa.fit(trainBags, trainLabels, type='average')
	# predictions = SMILa.predict(testBags)
	# accuracie = np.average(testLabels.T == np.sign(predictions))
	# print '\n Accuracy: %.2f%%' % (100 * accuracie)
	# fpr, tpr, thresholds = metrics.roc_curve(testLabels, predictions, pos_label=1.)
	# print metrics.auc(fpr, tpr)
	# 
	# totalAccuracy += 100*accuracie
	# totalAuc += metrics.auc(fpr,tpr)

# SMILa = simpleMIL()
# parameters_smil = {'type': 'average'}
# #En este me funciono maxDD porque no tiene problem con parametros 
# accuracie, results_accuracie, auc,results_auc,elapsed  = mil_cross_val(bags=bags,labels=np.array(labels), model=SMILa, folds=10, parameters=parameters_smil,timer=True)
# print accuracie
# print results_accuracie
# print auc
# print results_auc
# 
# exit()
# 
# # 
# # print np.array(bags[0])
# # exit()
# 
# #Make labels
# 
# #Run the MILES classifier
# import misvm
# #classifier = misvm.MISVM(kernel='linear', C=1.0, max_iters=50)
# classifier = misvm.SIL(kernel='linear', C=1.0)
# print "training classifier: "
# classifier.fit(trainBags, trainLabels)
# print "testing classifier: "
# predLabels = classifier.predict(testBags)
# 
# print np.average(testLabels == np.sign(predLabels))

#Visualize the results