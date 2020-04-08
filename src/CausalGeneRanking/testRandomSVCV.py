"""

	Compare the positive bags to bags of random SVs. Any better?

"""


import sys
import os
import numpy as np
import pickle as pkl
import random


import matplotlib
#matplotlib.use('Agg')

svTypes = ['DEL']
outDir = sys.argv[3]

#input the normalized bags
with open(sys.argv[1], 'rb') as handle:
	bagDict = pkl.load(handle)
with open(sys.argv[2], 'rb') as handle:
	bagDictRand = pkl.load(handle)

#load the bags of the randmo SVs

#get the information for the bag labels
degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object') #labels

print("initial number of bags: ", len(bagDict))
print('initial deg pairs: ', degPairs.shape[0])

#function to get the similarity matrix
def getSimilarityMatrix(bags, instances, reverseBagMap):

	bagIndices = np.arange(bags.shape[0])
	similarityMatrix = np.zeros([bags.shape[0], instances.shape[0]])
	print("Number of bags: ", bags.shape[0])
	for bagInd in range(0, bags.shape[0]):

		#Get the indices of the instances that are in this bag
		instanceIndices = reverseBagMap[bagInd]
		instanceSubset = instances[instanceIndices,:]

		#get the average of all instances in this bag
		instanceAvg = np.mean(instanceSubset, axis=0)

		#compute distance to all other instances from this bag average
		distance = np.abs(instanceAvg - instances)

		#sum the distances to get 1 similarity score
		summedDistance = np.sum(distance,axis=1)
		similarityMatrix[bagInd,:] = summedDistance

	return similarityMatrix

def getSimilarityMatrixTest(testBags, trainInstances, labels):

	similarityMatrix = np.zeros([testBags.shape[0], trainInstances.shape[0]])

	#print(similarityMatrix.shape)

	for bagInd in range(0, testBags.shape[0]):
		#print(labels[bagInd])
		#get the average of all instances in this test patient bag
		testInstances = testBags[bagInd]

		instanceAvg = np.mean(testInstances, axis=0)

		#compute distance to all other instances from this bag average
		distance = np.abs(instanceAvg - trainInstances)

		#sum the distances to get 1 similarity score
		summedDistance = np.sum(distance,axis=1)
		#print(summedDistance)
		similarityMatrix[bagInd,:] = summedDistance

	return similarityMatrix



#then, generate the similarity matrices for the SV types
for svType in svTypes:

	#allow for running with feature selection
	bagLabels = []
	positiveBagPairNames = []
	negativeBagPairNames = []
	positiveInstanceLabels = []
	positiveBags = []
	negativeBags = []
	patientNegativeBags = []
	patientNegativeBagPairNames = []

	#for each SV-gene pair, get the instances
	for pair in bagDict:

		#check if the SV type matches our selection
		splitPair = pair.split("_")
		shortPair = splitPair[7] + '_' + splitPair[0]

		if svType != '':
			if splitPair[12] != svType:
				continue

		#get the label of the bag by checking if it exists in degPairs, some pairs do not have a z-score because the gene is excluded due to mutations.
		if shortPair in degPairs[:,0]:

			#get the z-score of the pair.
			degPairInfo = degPairs[degPairs[:,0] == shortPair][0]

			#if the z-score matches this criterion, the SV-gene pair is positive
			if float(degPairInfo[5]) > 1.5 or float(degPairInfo[5]) < -1.5:
			#if float(degPairInfo[5]) > 1 or float(degPairInfo[5]) < -1:
				#go through the instances of this SV-gene pair, and include only those that have gains and losses, and more than 1 instance. This should in principle not happen, but good to keep a check.
				instances = []
				for instance in bagDict[pair]:

					if instance[0] == 0 and instance[1] == 0:
						continue


					if instance[34] == 1:
						continue

					#feature selection by hand
					featureInd = -1
					newInstance = []
					for feature in instance:
						featureInd += 1
						#if featureInd > 36:

						#	continue
						newInstance.append(feature)


					#instances.append(instance)
					instances.append(newInstance)

				if len(instances) < 1:
					continue


				positiveBagPairNames.append(pair)
				positiveBags.append(instances)

			else:
				instances = []
				for instance in bagDict[pair]:

					if instance[0] == 0 and instance[1] == 0:
						continue

					if instance[34] == 1:
						continue

					#feature selection by hand
					newInstance = []
					featureInd = -1
					for feature in instance:
						featureInd += 1
						#if featureInd > 36:
						#	continue
						newInstance.append(feature)


					#instances.append(instance)
					instances.append(newInstance)



				if len(instances) < 1:
					continue


				patientNegativeBagPairNames.append(pair)
				patientNegativeBags.append(instances)

	for pair in bagDictRand:

		splitPair = pair.split("_")
		shortPair = splitPair[7] + '_' + splitPair[0]

		if shortPair in degPairs[:,0]:

			degPairInfo = degPairs[degPairs[:,0] == shortPair][0]

			#if the z-score matches this criterion, the SV-gene pair is positive
			if float(degPairInfo[5]) > 1.5 or float(degPairInfo[5]) < -1.5:
				continue

		if svType != '':
			if splitPair[12] != svType:
				continue

		#get the right number of features per instance
		instances = []
		for instance in bagDictRand[pair]:

			if instance[0] == 0 and instance[1] == 0:
				continue
			
			

			if instance[34] == 1:
				continue

			newInstance = []
			featureInd = -1
			for feature in instance:
				
				featureInd += 1
				#if featureInd > 36:
				#	continue
				newInstance.append(feature)
				
				

			#instances.append(instance)
			instances.append(newInstance)

		if len(instances) < 1:
			continue

		negativeBags.append(instances)
		negativeBagPairNames.append(pair)


	positiveBags = np.array(positiveBags)
	negativeBags = np.array(negativeBags)
	positiveBagPairNames = np.array(positiveBagPairNames)
	negativeBagPairNames = np.array(negativeBagPairNames)
	patientNegativeBags = np.array(patientNegativeBags)
	patientNegativeBagPairNames = np.array(patientNegativeBagPairNames)

	print(positiveBags.shape)
	print(negativeBags.shape)
	
	instPerBagPos = dict()
	instPerBagNeg = dict()
	instPerBagRand = dict()
	for bag in positiveBags:
		instCount = len(bag)
		if instCount not in instPerBagPos:
			instPerBagPos[instCount] = 0
		instPerBagPos[instCount] += 1

		enhCount = 0
		promCount = 0
		eQTLCount = 0
		seCount = 0
		for instance in bag:
			
			if instance[36] == 1:
				seCount += 1
			if instance[35] == 1:
				eQTLCount += 1
			if instance[34] == 1:
				promCount += 1
			if instance[33] == 1:
				enhCount += 1

		for instance in bag:
			instance.append(instCount / positiveBags.shape[0])
			#instance.append(enhCount / positiveBags.shape[0])
			#instance.append(promCount / positiveBags.shape[0])
			#instance.append(eQTLCount / positiveBags.shape[0])
			#instance.append(seCount / positiveBags.shape[0])

	for bag in negativeBags:
		instCount = len(bag)
		if instCount not in instPerBagNeg:
			instPerBagNeg[instCount] = 0
		instPerBagNeg[instCount] += 1

		for instance in bag:
			instance.append(instCount / negativeBags.shape[0])
		# enhCount = 0
		# promCount = 0
		# eQTLCount = 0
		# seCount = 0
		# for instance in bag:
		# 	if instance[36] == 1:
		# 		seCount += 1
		# 	if instance[35] == 1:
		# 		eQTLCount += 1
		# 	if instance[34] == 1:
		# 		promCount += 1
		# 	if instance[33] == 1:
		# 		enhCount += 1
		#
		# for instance in bag:

			#instance.append(enhCount / negativeBags.shape[0])
			#instance.append(promCount / negativeBags.shape[0])
			#instance.append(eQTLCount / negativeBags.shape[0])
			#instance.append(seCount / negativeBags.shape[0])

	for bag in patientNegativeBags:
		instCount = len(bag)
		if instCount not in instPerBagRand:
			instPerBagRand[instCount] = 0
		instPerBagRand[instCount] += 1

		for instance in bag:
			instance.append(instCount / patientNegativeBags.shape[0])
		# enhCount = 0
		# promCount = 0
		# eQTLCount = 0
		# seCount = 0
		# for instance in bag:
		# 	if instance[36] == 1:
		# 		seCount += 1
		# 	if instance[35] == 1:
		# 		eQTLCount += 1
		# 	if instance[34] == 1:
		# 		promCount += 1
		# 	if instance[33] == 1:
		# 		enhCount += 1
		#
		# for instance in bag:
		# 	#instance.append(instCount / positiveBags.shape[0])
		# 	instance.append(enhCount / patientNegativeBags.shape[0])
		# 	instance.append(promCount / patientNegativeBags.shape[0])
		# 	instance.append(eQTLCount / patientNegativeBags.shape[0])
		# 	instance.append(seCount / patientNegativeBags.shape[0])

	print(instPerBagPos)
	print(instPerBagNeg)
	print(instPerBagRand)
	
	negativeBags = patientNegativeBags
	negativeBagPairNames = patientNegativeBagPairNames
	

	posInstances = np.vstack(positiveBags)
	negInstances = np.vstack(negativeBags)

	allInstances = np.concatenate((posInstances, negInstances))
	#check variances within instances
	from sklearn.feature_selection import VarianceThreshold
	t = (.8 * (1 - .8))
	t = 0
	vt = VarianceThreshold(threshold=t)
	vt.fit(allInstances)
	idx = np.where(vt.variances_ > t)[0]

	newPositiveBags = []
	newNegativeBags = []
	for bag in positiveBags:
		instances = []
		for instance in bag:
			filteredInstance = []
			featureInd = 0
			for feature in instance:
				if featureInd in idx:
					filteredInstance.append(feature)
				featureInd += 1
			instances.append(filteredInstance)

		newPositiveBags.append(instances)

	for bag in negativeBags:
		instances = []
		for instance in bag:
			filteredInstance = []
			featureInd = 0
			for feature in instance:
				if featureInd in idx:
					filteredInstance.append(feature)
				featureInd += 1
			instances.append(filteredInstance)

		newNegativeBags.append(instances)

	positiveBags = np.array(newPositiveBags)
	negativeBags = np.array(newNegativeBags)
	print(positiveBags.shape)
	print(negativeBags.shape)

	posInstances = np.vstack(positiveBags)
	negInstances = np.vstack(negativeBags)

	#what are the feature distributions of these sets?
	goodInstancesSum = np.sum(posInstances, axis=0)
	goodInstancesAvg = goodInstancesSum / posInstances.shape[0]
	badInstancesSum = np.sum(negInstances, axis=0)
	badInstancesAvg = badInstancesSum / negInstances.shape[0]

	print(goodInstancesAvg)
	print(badInstancesAvg)

	#xlabels = np.array(xlabels)
	#xlabels = xlabels[allowedList]

	import matplotlib.pyplot as plt

	barWidth = 0.35
	plt.bar(range(0, len(goodInstancesAvg)), goodInstancesAvg, width=barWidth, color='blue')
	#plt.xticks(range(0, len(xlabels)), xlabels, rotation=90)
	#plt.show()
	pos = np.arange(len(goodInstancesAvg))
	r2 = [i + barWidth for i in pos]
	plt.bar(r2, badInstancesAvg, color='orange', width=barWidth)
	#plt.xticks(r2, xlabels, rotation=90)
	plt.show()

	

#
# ####lopocv tesing
# perPatientPositiveBags = dict()
# for bagInd in range(0, positiveBags.shape[0]):
#
# 	#get the label of this bag
# 	bagPairLabel = positiveBagPairNames[bagInd]
# 	splitLabel = bagPairLabel.split('_')
#
# 	patientId = splitLabel[7]
# 	if patientId not in perPatientPositiveBags:
# 		perPatientPositiveBags[patientId] = []
#
# 	perPatientPositiveBags[patientId].append(positiveBags[bagInd])
#
# perPatientNegativeBags = dict()
# for bagInd in range(0, negativeBags.shape[0]):
#
# 	#get the label of this bag
# 	bagPairLabel = negativeBagPairNames[bagInd]
# 	splitLabel = bagPairLabel.split('_')
#
#
# 	patientId = splitLabel[7]
# 	if patientId not in perPatientNegativeBags:
# 		perPatientNegativeBags[patientId] = []
#
# 	perPatientNegativeBags[patientId].append(negativeBags[bagInd])
#
# #then, for each patient, use it as the test set.
# aucs = []
# performances = []
# posPerformances = []
# for patient in perPatientPositiveBags:
# 	print(patient)
#
# 	if patient not in perPatientNegativeBags:
# 		continue
#
# 	#Get as many negative bags from this patient as there are positive bags.
# 	if len(perPatientNegativeBags[patient]) > len(perPatientPositiveBags[patient]):
# 		randInd = random.sample(range(0, len(perPatientNegativeBags[patient])), len(perPatientPositiveBags[patient]))
#
# 		patientPositiveBags = []
# 		for bag in perPatientPositiveBags[patient]:
#
# 			patientPositiveBags.append(bag)
#
#
# 		#patientPositiveBags = perPatientPositiveBags[patient]
# 		#patientPositiveBags = np.array(patientPositiveBags)
# 		patientNegativeBags = perPatientNegativeBags[patient]
# 		patientNegativeBags = np.array(patientNegativeBags)[randInd]
# 		patientNegativeBags = list(patientNegativeBags)
# 	else:
# 		#get as any positives as there are negatives
# 		randInd = random.sample(range(0, len(perPatientPositiveBags[patient])), len(perPatientNegativeBags[patient]))
#
# 		patientNegativeBags = []
# 		for bag in perPatientNegativeBags[patient]:
#
# 			patientNegativeBags.append(bag)
#
#
# 		#patientPositiveBags = perPatientPositiveBags[patient]
# 		#patientPositiveBags = np.array(patientPositiveBags)
# 		patientPositiveBags = perPatientPositiveBags[patient]
# 		patientPositiveBags = np.array(patientPositiveBags)[randInd]
# 		patientPositiveBags = list(patientPositiveBags)
#
#
# 	#print(patientPositiveBags.shape)
# 	#print(patientNegativeBags.shape)
#
# 	patientBags = patientPositiveBags + patientNegativeBags
# 	print(len(patientBags))
# 	patientBags = np.array(patientBags)
# 	print(patientBags.shape)
#
# 	patientPositiveBags = np.array(patientPositiveBags)
# 	patientNegativeBags = np.array(patientNegativeBags)
# 	#then, get the bags of all other patients.
# 	allPositiveBags = []
# 	allPossibleNegativeBags = []
# 	patientCount = 0
# 	for patient2 in perPatientPositiveBags:
#
# 		if patient2 == patient:
# 			continue
# 		for bag in perPatientPositiveBags[patient2]:
# 			allPositiveBags.append(bag)
#
# 		if patient2 not in perPatientNegativeBags:
# 			continue
# 		for bag in perPatientNegativeBags[patient2]:
# 			allPossibleNegativeBags.append(bag)
#
#
# 		#print(len(perPatientPositiveBags[patient2]))
# 		#print(patientNegativeBags.shape)
#
# 		patientCount += 1
#
# 	allPositiveBags = np.array(allPositiveBags)
# 	print(allPositiveBags.shape)
#
# 	allPossibleNegativeBags = np.array(allPossibleNegativeBags)
# 	randInd = random.sample(range(0, allPossibleNegativeBags.shape[0]), allPositiveBags.shape[0])
# 	allNegativeBags = allPossibleNegativeBags[randInd]
# 	#get the same number of negative bags, sample randomly.
# 	print(allNegativeBags.shape)
#
#
#
# 	trainBags = np.concatenate((allPositiveBags, allNegativeBags))
# 	trainInstances = np.vstack(trainBags)
# 	#testBags = np.concatenate((patientPositiveBags, patientNegativeBags))
# 	testBags = patientBags
#
# 	trainLabels = [1]*allPositiveBags.shape[0] + [0]*allNegativeBags.shape[0]
# 	testLabels = [1]*patientPositiveBags.shape[0] + [0]*patientNegativeBags.shape[0]
#
# 	reverseBagMapOtherPatients = dict() #lookup instance by bag index
# 	instanceInd = 0
# 	for bagInd in range(0, trainBags.shape[0]):
# 		reverseBagMapOtherPatients[bagInd] = []
# 		for instance in trainBags[bagInd]:
# 			reverseBagMapOtherPatients[bagInd].append(instanceInd)
# 			instanceInd += 1
#
# 	#make similarity matrices
# 	similarityMatrixTrain = getSimilarityMatrix(trainBags, trainInstances, reverseBagMapOtherPatients)
# 	print(similarityMatrixTrain.shape)
# 	#now the curent patient bags need to be to the instances of the training set
# 	similarityMatrixTest = getSimilarityMatrixTest(testBags, trainInstances, testLabels)
# 	print(similarityMatrixTest.shape)
#
# 	from sklearn.ensemble import RandomForestClassifier
# 	from sklearn.model_selection import StratifiedKFold
# 	from sklearn import model_selection
# 	from sklearn.metrics import plot_roc_curve, auc
# 	import matplotlib.pyplot as plt
# 	from scipy import interp
#
# 	classifier = RandomForestClassifier(n_estimators= 100)
# 	#then train the classifier
# 	classifier.fit(similarityMatrixTrain, trainLabels)
# 	print(testLabels)
# 	print(classifier.predict(similarityMatrixTest))
#
# 	preds = classifier.predict(similarityMatrixTrain)
# 	diff = np.sum(np.abs(trainLabels - preds)) / len(trainLabels)
# 	print('train diff: ', diff)
#
# 	preds = classifier.predict(similarityMatrixTest)
# 	diff = np.sum(np.abs(testLabels - preds)) / len(testLabels)
# 	print('test diff: ', diff)
#
# 	#how many positives are correct?
# 	correctPos = 0
# 	allPos = 0
# 	for label in range(0, len(testLabels)):
# 		if testLabels[label] == 1 and preds[label] == 1:
# 			correctPos += 1
#
# 		if testLabels[label] == 1:
# 			allPos += 1
#
# 	print('correct positives: ', correctPos / allPos)
# 	posPerformances.append(correctPos/allPos)
#
# 	print('train: ', classifier.score(similarityMatrixTrain, trainLabels))
# 	print('test: ', classifier.score(similarityMatrixTest, testLabels))
# 	performances.append(classifier.score(similarityMatrixTest, testLabels))
#
# 	fig, ax = plt.subplots()
# 	viz = plot_roc_curve(classifier, similarityMatrixTest, testLabels,
# 						 name='roc',
# 						 alpha=0.3, lw=1, ax=ax)
# 	aucs.append(np.mean(viz.roc_auc))
# 	print('auc: ', np.mean(viz.roc_auc))
#
#
#
# print(np.mean(aucs))
# print(np.mean(performances))
# print(np.mean(posPerformances))
#
#
#
# exit()




#clf testing
chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
						   'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
						   'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
						   'chr20', 'chr21', 'chr22']


positiveBagsPerChromosome = dict()
negativeBagsPerChromosome = dict()
for labelInd in range(0, len(positiveBagPairNames)):
	label = positiveBagPairNames[labelInd]
	splitLabel = label.split('_')

	chromosome = splitLabel[1]
	if chromosome not in positiveBagsPerChromosome:
		positiveBagsPerChromosome[chromosome] = []
	positiveBagsPerChromosome[chromosome].append(positiveBags[labelInd])

for labelInd in range(0, len(negativeBagPairNames)):
	label = negativeBagPairNames[labelInd]
	splitLabel = label.split('_')

	chromosome = splitLabel[1]
	if chromosome not in negativeBagsPerChromosome:
		negativeBagsPerChromosome[chromosome] = []
	negativeBagsPerChromosome[chromosome].append(negativeBags[labelInd])


#performance test
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn import model_selection
from sklearn.metrics import plot_roc_curve, auc, average_precision_score
import matplotlib.pyplot as plt
from scipy import interp
#train simple rf, check performance
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
fig, ax = plt.subplots()
performances = []
posPerformances = []
allPreds = []
allTrueLabels = []

classifier = RandomForestClassifier(n_estimators= 200, random_state=42)

trainBags = dict()
testBags = dict()
trainLabels = dict()
testLabels = dict()
for chromosome in chromosomes:
	print(chromosome)
	
	if chromosome not in positiveBagsPerChromosome:
		continue
	if chromosome not in negativeBagsPerChromosome:
		continue

	#if chromosome != 'chr16':
	#	continue


	#make stratified
	testPositiveBags = positiveBagsPerChromosome[chromosome]
	testNegativeBags = negativeBagsPerChromosome[chromosome]

	testPositiveBags = np.array(testPositiveBags)
	testNegativeBags = np.array(testNegativeBags)

	random.seed(785)
	randInd = random.sample(range(0, testNegativeBags.shape[0]), testPositiveBags.shape[0])
	testSubsetNegativeBags = testNegativeBags[randInd]

	allTestBags = []
	for bag in testPositiveBags:
		allTestBags.append(bag)
	for bag in testSubsetNegativeBags:
	#for bag in testNegativeBags:
		allTestBags.append(bag)

	allTestBags = np.array(allTestBags)


	print(allTestBags.shape)


	testBags[chromosome] = allTestBags
	testLabels[chromosome] = [1]*testPositiveBags.shape[0] + [0]*testSubsetNegativeBags.shape[0]
	#testLabels[chromosome] = [1]*testPositiveBags.shape[0] + [0]*testNegativeBags.shape[0]


	labels = [1]*testPositiveBags.shape[0] + [0]*testSubsetNegativeBags.shape[0]
	allTrueLabels += labels

	testPositiveInstances = np.vstack(testPositiveBags)
	testNegativeInstances = np.vstack(testSubsetNegativeBags)
	testPositiveLabels = [1]*testPositiveInstances.shape[0]
	testNegativeLabels = [0]*testNegativeInstances.shape[0]

	allTestInstances = np.concatenate((testPositiveInstances, testNegativeInstances))
	allTestLabels = testPositiveLabels + testNegativeLabels


	#make training set from the rest
	trainingSet = []
	trainingLabels = []
	allTrainInstances = []
	allTrainLabels = []
	for chromosome2 in chromosomes:

		if chromosome == chromosome2:
			continue

		if chromosome2 not in positiveBagsPerChromosome:
			continue
		if chromosome2 not in negativeBagsPerChromosome:
			continue

		#make stratified
		chrPositiveBags = positiveBagsPerChromosome[chromosome2]
		chrNegativeBags = negativeBagsPerChromosome[chromosome2]

		chrPositiveBags = np.array(chrPositiveBags)
		chrNegativeBags = np.array(chrNegativeBags)

		random.seed(785)
		randInd = random.sample(range(0, chrNegativeBags.shape[0]), chrPositiveBags.shape[0])
		subsetNegativeBags = chrNegativeBags[randInd]

		for bag in chrPositiveBags:
			trainingSet.append(bag)

		for bag in subsetNegativeBags:
		#for bag in chrNegativeBags:
			trainingSet.append(bag)

		trainingLabels += [1]*chrPositiveBags.shape[0]
		trainingLabels += [0]*subsetNegativeBags.shape[0]
		#trainingLabels += [0]*chrNegativeBags.shape[0]

		trainPositiveInstances = np.vstack(chrPositiveBags)
		trainNegativeInstances = np.vstack(subsetNegativeBags)
		trainPositiveLabels = [1]*trainPositiveInstances.shape[0]
		trainNegativeLabels = [0]*trainNegativeInstances.shape[0]

		for instance in trainPositiveInstances:
			allTrainInstances.append(instance)
		for instance in trainNegativeInstances:
			allTrainInstances.append(instance)
		allTrainLabels += trainPositiveLabels
		allTrainLabels += trainNegativeLabels

	allTrainInstances = np.array(allTrainInstances)

	allTestInstances = np.concatenate((testPositiveInstances, testNegativeInstances))
	allTestLabels = testPositiveLabels + testNegativeLabels

	trainBags[chromosome] = np.array(trainingSet)
	trainLabels[chromosome] = trainingLabels


	trainInstances = np.vstack(trainBags[chromosome])

	reverseBagMapOtherPatients = dict() #lookup instance by bag index
	instanceInd = 0
	for bagInd in range(0, trainBags[chromosome].shape[0]):
		reverseBagMapOtherPatients[bagInd] = []
		for instance in trainBags[chromosome][bagInd]:
			reverseBagMapOtherPatients[bagInd].append(instanceInd)
			instanceInd += 1

	similarityMatrixTrain = getSimilarityMatrix(trainBags[chromosome], trainInstances, reverseBagMapOtherPatients)
	print(similarityMatrixTrain.shape)
	#now the curent patient bags need to be to the instances of the training set
	similarityMatrixTest = getSimilarityMatrixTest(testBags[chromosome], trainInstances, testLabels)
	print(similarityMatrixTest.shape)


	#then train the classifier
	classifier.fit(similarityMatrixTrain, trainLabels[chromosome])

	trainPreds = classifier.predict(similarityMatrixTrain)
	diff = np.sum(np.abs(trainLabels[chromosome] - trainPreds)) / len(trainLabels[chromosome])
	print('train diff: ', diff)

	preds = classifier.predict(similarityMatrixTest)
	print(testLabels[chromosome])
	print(preds)
	allPreds += list(preds)
	diff = np.sum(np.abs(testLabels[chromosome] - preds)) / len(testLabels[chromosome])
	print('test diff: ', diff)


	predProb = classifier.predict_proba(similarityMatrixTest)
	#print('real auc: ', auc(testLabels[chromosome], list(predProb[:,1])))

	print('train: ', classifier.score(similarityMatrixTrain, trainLabels[chromosome]))
	print('test: ', classifier.score(similarityMatrixTest, testLabels[chromosome]))
	performances.append(classifier.score(similarityMatrixTest, testLabels[chromosome]))

	fig, ax = plt.subplots()
	viz = plot_roc_curve(classifier, similarityMatrixTest, testLabels[chromosome],
						 name='roc',
						 alpha=0.3, lw=1, ax=ax)
	interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
	interp_tpr[0] = 0.0
	tprs.append(interp_tpr)
	aucs.append(np.mean(viz.roc_auc))
	print('auc: ', np.mean(viz.roc_auc))

	print('AP: ', average_precision_score(testLabels[chromosome], predProb[:,1]))

	#exit()
print(np.mean(performances))
print(np.mean(aucs))

