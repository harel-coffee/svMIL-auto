#Run 1 iteration of the SMIL classifier on randomized labels


#Shuffle the labels
import random
import os
import sys
import numpy as np
from sklearn import metrics
sys.path.append(os.path.realpath('../../../../Methods/'))
from MILpy.Algorithms.simpleMIL import simpleMIL
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from scipy import interp
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.utils.fixes import signature


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
	
	bagLabel = -1
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
	

# # Somatic vs germline
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
# 	labels.append(-1)

bags = np.array(bags)

import misvm
#classifier = misvm.MISVM(kernel='linear', C=1.0, max_iters=50)
classifier = misvm.sbMIL(kernel='linear', eta=0.1, C=1e2)

classifier.fit(bags, labels)
bag_labels, instance_labels = classifier.predict(bags, instancePrediction=True)

print len(labels)
print len(np.where(instance_labels == 1)[0])
exit()



#Plot ROC & precision-recall
# 
# tprs = []
# aucs = []
# totalAccuracy = 0
# mean_fpr = np.linspace(0, 1, 100)
# 
# testReal = []
# testProbs = []
# 
# #random.shuffle(labels)
# labels = np.array(labels)
# 
# cv = StratifiedKFold(n_splits=10)
# SMILa = simpleMIL()
# 
# #Each CV split will output 1 AUC, 1 accuracy for the ROC
# #For precision-recall, each split will output precitions & we have real labels. This is what we store. 
# 
# aucs = []
# accuracies = []
# testReal = []
# testProbs = []
# tprs = []
# for train, test in cv.split(bags, labels):
# 
# 	SMILa = simpleMIL() 
# 	SMILa.fit(bags[train], labels[train], type='average')
# 	predictions = SMILa.predict(bags[test])
# 	predictionProbs = SMILa.predict_probs(bags[test])[:,1]
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
# 	# #precision-recall
# 	# precision, recall, thresholds = precision_recall_curve(labels[test], predictions)
# 	# print precision
# 	# print recall
# 	testReal.append(labels[test])
# 	testProbs.append(predictionProbs)
# # 
# # tprs.append(np.mean(iterationTprs, axis=0))
# # aucs.append(np.mean(iterationAucs))
# 
# #For ROC, we simply average the AUC and accuracies.
# avgAuc = np.mean(aucs)
# avgAccuracy = np.mean(accuracies)
# meanTpr = np.mean(tprs, axis=0)
# 
# #For precision-recall, we can only make the curve at the end, because toherwise there are different thresholds. SO we need to write the labels and predictions to the file. 
# testReal = np.concatenate(testReal)
# testProbs = np.concatenate(testProbs)
# 
# precision, recall, _ = precision_recall_curve(testReal, testProbs)
# 
# with open(sys.argv[3] + '/' + sys.argv[4] + ".txt", 'w') as outF:
# 	outF.write(str(avgAuc) + "\n")
# 	outF.write(str(avgAccuracy) + "\n")
# 	outF.write(",".join([str(i) for i in meanTpr]) + "\n")
# 	outF.write(",".join([str(i) for i in testReal]) + "\n")
# 	outF.write(",".join([str(i) for i in testProbs]) + "\n")










