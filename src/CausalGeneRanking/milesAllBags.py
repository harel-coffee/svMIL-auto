import sys
import numpy as np
import matplotlib.pyplot as plt


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


## Running MILES with somatic only, but using per patient-gene pair DEG labels
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


#### Germline vs somatic miles

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
# 	labels.append(-1)


# bags = bags[1:100]
# labels = labels[1:100]
bags = np.array(bags)
instances = np.vstack(bags)


print "generating similarity matrix"

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


similarityMatrix = np.zeros([bags.shape[0], instances.shape[0]])
print "Number of bags: ", bags.shape[0]
for bagInd in range(0, bags.shape[0]):
	
	print bagInd
	
	#Get the indices of the instances that are in this bag
	instanceIndices = reverseBagMap[bagInd]
	
	instanceSubset = instances[instanceIndices,:]
	
	#Compute the pairwise distance matrix here
	minDistance = float("inf")
	minDistanceInd = 0
	for instanceInd in range(0, instanceSubset.shape[0]):
		instance = instanceSubset[instanceInd]
		distance = np.abs(instance - instances) #compute the distances to the train instances, otherwise we are not in the same similarity space. 

		summedDistance = np.sum(distance,axis=1)

		currentMinDistance = np.min(summedDistance)
		if currentMinDistance < np.min(minDistance):
			minDistance = summedDistance
			minDistanceInd = instanceInd

	#This instance will be used as representative for this bag. We use this value as the similarity to all other instances.  
	similarityMatrix[bagInd] = minDistance

#4. Train a classifier on the similarity space
from sklearn.ensemble import RandomForestClassifier
print "training the classifier in similarity space"
np.random.seed(500)
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
rfClassifier.fit(similarityMatrix, labels) #Use the bag labels, not the instance labels

print similarityMatrix
print labels

predictions = rfClassifier.predict(similarityMatrix)
print predictions
print np.average(labels == np.sign(predictions))

#Now select the most important features with random forest
importances = rfClassifier.feature_importances_
std = np.std([tree.feature_importances_ for tree in rfClassifier.estimators_],
             axis=0)
indices = np.argsort(importances)[::-1]


nonZeroIndices = []
for index in indices:
	if importances[index] > 0:
		nonZeroIndices.append(index)

#Get the genes at these indices
positiveGenes = dict()
positivePairs = []
for index in nonZeroIndices:
	pairName = pairNames[index]
	positivePairs.append(pairName)
	splitPairName = pairName.split("_")
	geneName = splitPairName[len(splitPairName)-1]
	if geneName not in positiveGenes:
		positiveGenes[geneName] = 0
	positiveGenes[geneName] += 1

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
# print "Used features: ", len(clf.coef_[0] != 0)


print len(positivePairs)
print len(positiveGenes)

milesConceptGenesOut = "SomaticGermline/milesConceptGenes.txt"
with open(milesConceptGenesOut, 'w') as outF:
	for gene in positiveGenes:
		outF.write(gene + "\t" + str(positiveGenes[gene]) + "\n")

#output the similarity matrix and also the labels and bag indices so that we can do analysis after running once

np.save("SomaticGermline/bags.txt", bags)
np.save("SomaticGermline/pairNames.txt", pairNames) #the sv-gene pair names of each bag entry
np.save("SomaticGermline/labels.txt", labels)
np.save("SomaticGermline/similarityMatrix.txt", similarityMatrix)
np.save("SomaticGermline/conceptIndices.txt", nonZeroIndices)


exit()


# Print the feature ranking
print("Feature ranking:")

for f in range(500):
    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))
exit()

# Plot the feature importances of the forest
plt.figure()
plt.title("Feature importances")
plt.bar(range(instances.shape[0]), importances[indices],
       color="r", yerr=std[indices], align="center")
plt.xticks(range(instances.shape[0]), indices)
plt.xlim([-1, instances.shape[0]])
plt.show()



