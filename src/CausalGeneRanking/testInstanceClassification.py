from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
from six.moves import range


somaticScores = np.loadtxt(sys.argv[1], dtype="object")
germlineScores = np.loadtxt(sys.argv[2], dtype="object")

#Use these to filter the gene-SV pairs by the number of smaples affecting them in total
somaticAllSampleScores = np.loadtxt(sys.argv[3], dtype="object")
germlineAllSampleScores = np.loadtxt(sys.argv[4], dtype="object")

#Here already divide the data into train & test. Keep a set separate for both the somatic and germline that where we do not do the selection on, e.g. 40%.
#Then these pairs are removed from the scores, and we use the rest to make training bags, which we filter for samples > 6. 

#Take a random sample with 40% of both the somatic and germline scores, keep this data separate

#Make separate test bags for this data

#Make bags as usual for the rest of the data

#Do not include genes that do not have > 6 samples

#Train the classifier on the training bags

#Test on the left out data



#Make bags

#Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted
#A bag is positive if it comes from the somatic set, and negative if it comes from the germline set

#First assign the feature vectors to the right SV in order

svBagContents = dict()
genesPerBag = dict() #genes per bag, to later link back to instances in the correct order
bagSamples = dict()
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
		
		samples = somaticAllSampleScores[somaticAllSampleScores[:,0] == geneName, 31][0]
		splitSamples = samples.split(",")
		if len(splitSamples) > 6:
			bagSamples[geneName] = True
		else:
			bagSamples[geneName] = False
	

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
		
		samples = somaticAllSampleScores[somaticAllSampleScores[:,0] == geneName, 31][0]
		splitSamples = samples.split(",")
		if len(splitSamples) > 6:
			bagSamples[geneName] = True
		else:
			bagSamples[geneName] = False


#Then construct the bags

#bags = []
for sv in svBagContents:
	for gene in genesPerBag[sv]:	
		pairNames.append(sv + "_" + gene)
	bags.append(svBagContents[sv])
	labels.append(-1)

bags = np.array(bags)
from sklearn.model_selection import train_test_split
indices = list(range(0,len(bags)))



# import misvm
# classifier = misvm.SIL(kernel='linear', C=1.0)
# classifier.fit(trainBags, trainLabels)
# bag_labels, instance_labels = classifier.predict(testBags, instancePrediction=True)
# print instance_labels
import random
import os
from sklearn import metrics
sys.path.append(os.path.realpath('../../../../Methods/'))
from MILpy.Algorithms.simpleMIL import simpleMIL

X_train, X_test, y_train, y_test = train_test_split(bags, labels, test_size=.4, random_state=42) #60%/40%
instances = []
instanceLabels = []
for bagInd in range(0, testBags.shape[0]):
	bag = testBags[bagInd]
	
	for instance in bag:
		instances.append([instance])
		instanceLabels.append(testLabels[bagInd])

instances = np.array(instances)

SMILa = simpleMIL() 
SMILa.fit(trainBags, trainLabels, type='average')
predictions = SMILa.predict(testBags)
predictions = SMILa.predict(instances)


print(np.average(instanceLabels == np.sign(predictions)))
print(np.where(np.array(predictions) == 1)[0].shape)
print(np.where(np.array(predictions) == -1)[0].shape)

positiveInd = np.where(np.array(predictions) == 1)[0]
print(positiveInd)

positivePairs = []
positiveGenes = dict()
for ind in positiveInd:
	positivePairs.append(pairNames[ind])
	splitPairName = pairNames[ind].split("_")
	geneName = splitPairName[len(splitPairName)-1]
	if geneName not in positiveGenes:
		positiveGenes[geneName] = 0
	positiveGenes[geneName] += 1

print(len(positiveGenes))

fpr, tpr, thresholds = metrics.roc_curve(instanceLabels, predictions, pos_label=1.)
print(metrics.auc(fpr, tpr))


