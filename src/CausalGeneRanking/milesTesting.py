import sys
import numpy as np


somaticScores = np.loadtxt(sys.argv[1], dtype="object")
germlineScores = np.loadtxt(sys.argv[2], dtype="object")

#Make bags

#Each SV gets a new bag. In this bag are all the feature vectors of the genes that are disrupted
#A bag is positive if it comes from the somatic set, and negative if it comes from the germline set

#First assign the feature vectors to the right SV in order
svBagContents = dict()
for geneSVPair in somaticScores:
	
	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
	score = False
	for feature in features:
		if feature > 0:
			score = True
			break
	if score == True:	
		splitGeneSVPairName = geneSVPair[0].split("_")
		
		#The first element will be the gene name, the rest is the SV information
		splitGeneSVPairName.pop(0) #remove first element
		
		sv = "_".join(splitGeneSVPairName)
		if sv not in svBagContents:
			svBagContents[sv] = []
		
		svBagContents[sv].append(features)
	

#Then construct the bags

bags = []
labels = []
for sv in svBagContents:
	
	bags.append(svBagContents[sv])
	labels.append(1)
	
svBagContents = dict()
for geneSVPair in germlineScores:
	
	features = [float(i) for i in geneSVPair[1:len(geneSVPair)-1]]
	score = False
	for feature in features:
		if feature > 0:
			score = True
			break
	if score == True:	
		splitGeneSVPairName = geneSVPair[0].split("_")
		
		#The first element will be the gene name, the rest is the SV information
		splitGeneSVPairName.pop(0) #remove first element
		
		sv = "_".join(splitGeneSVPairName)
		if sv not in svBagContents:
			svBagContents[sv] = []
		
		svBagContents[sv].append(features)


#Then construct the bags

#bags = []
for sv in svBagContents:
	
	bags.append(svBagContents[sv])
	labels.append(-1)

bags = np.array(bags)
from sklearn.model_selection import train_test_split
trainBags, testBags, trainLabels, testLabels = train_test_split(bags, labels, test_size=.4, random_state=42) #60%/40%

import os
from sklearn import metrics
sys.path.append(os.path.realpath('../../../../Methods/'))
from MILpy.Algorithms.simpleMIL import simpleMIL
from MILpy.Algorithms.MILES import MILES
from MILpy.functions.mil_cross_val import mil_cross_val

testLabels = np.array(testLabels)

SMILa = simpleMIL() 
SMILa.fit(trainBags, trainLabels, type='average')
predictions = SMILa.predict(testBags)
accuracie = np.average(testLabels.T == np.sign(predictions))
print '\n Accuracy: %.2f%%' % (100 * accuracie)
fpr, tpr, thresholds = metrics.roc_curve(testLabels, predictions, pos_label=1.)
print metrics.auc(fpr, tpr)

SMILa = simpleMIL()
parameters_smil = {'type': 'max'}
#En este me funciono maxDD porque no tiene problem con parametros 
accuracie, results_accuracie, auc,results_auc,elapsed  = mil_cross_val(bags=bags,labels=np.array(labels), model=SMILa, folds=10, parameters=parameters_smil,timer=True)
print accuracie
print results_accuracie
print auc
print results_auc

exit()

# 
# print np.array(bags[0])
# exit()

#Make labels

#Run the MILES classifier
import misvm
#classifier = misvm.MISVM(kernel='linear', C=1.0, max_iters=50)
classifier = misvm.SIL(kernel='linear', C=1.0)
print "training classifier: "
classifier.fit(trainBags, trainLabels)
print "testing classifier: "
predLabels = classifier.predict(testBags)

print np.average(testLabels == np.sign(predLabels))

#Visualize the results