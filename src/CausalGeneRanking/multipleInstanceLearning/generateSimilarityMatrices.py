"""
	The goal of this script is to generate the similarity matrices for MIL

	If we do feature selection, we have similarity matrices for each feature
	output these as well.

"""

import sys
import os
import numpy as np
import pickle as pkl
import random


import matplotlib
#matplotlib.use('Agg')

featureElimination = sys.argv[2]
leaveOnePatientOut = sys.argv[3] #make the similarity matrices for each left out patient
svTypes = ['DEL', 'DUP', 'INV', 'ITX']
svTypes = ['DEL']

outDir = sys.argv[1]
finalOutDir = outDir + '/multipleInstanceLearning/similarityMatrices/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

if featureElimination == "True":
	featureEliminationOutDir = finalOutDir + '/featureSelection'
	if not os.path.exists(featureEliminationOutDir):
		os.makedirs(featureEliminationOutDir)

if leaveOnePatientOut == 'True':
	leaveOnePatientOutDir = finalOutDir + '/leaveOnePatientOut'
	if not os.path.exists(leaveOnePatientOutDir):
		os.makedirs(leaveOnePatientOutDir)

#input the normalized bags
#with open(outDir + '/multipleInstanceLearning/normalizedBags.pkl', 'rb') as handle:
#	bagDict = pkl.load(handle)

with open('./output/RankedGenes/22022020/BRCA/normalizedBags.pkl', 'rb') as handle:
	bagDict = pkl.load(handle)

#get the information for the bag labels
degPairs = np.loadtxt(outDir + '/tadDisruptionsZScores/zScores.txt', dtype='object') #labels

print(degPairs)

print("initial number of bags: ", len(bagDict))
print('initial deg pairs: ', degPairs.shape[0])

#function to get the similarity matrix
def getSimilarityMatrix(bags, instances, reverseBagMap):

	bagIndices = np.arange(bags.shape[0])
	similarityMatrix = np.zeros([bags.shape[0], instances.shape[0]])

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

	#for each SV-gene pair, get the instances
	for pair in bagDict:

		#check if the SV type matches our selection
		splitPair = pair.split("_")
		shortPair = splitPair[7] + '_' + splitPair[0]

		if svType != '':
			if splitPair[12] != svType:
				continue

		xlabels = ['loss', 'gain', 'cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
				   'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
				   'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'enhancerType', 'promoterType', 'eQTLType', 'superEnhancerType']
		#allowedList = [7,8, 9, 10, 11, 12, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
		allowedList = np.arange(len(xlabels))

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


					instances.append(instance)

				if len(instances) < 1:
					continue


				positiveBagPairNames.append(pair)
				positiveBags.append(instances)

			else: #if the z-score is anything else, this bag will be labeled negative.

				#get the right number of features per instance
				instances = []
				for instance in bagDict[pair]:

					if instance[0] == 0 and instance[1] == 0:
						continue


					instances.append(instance)
					
				if len(instances) < 1:
					continue

				negativeBags.append(instances)
				negativeBagPairNames.append(pair)


	positiveBags = np.array(positiveBags)
	negativeBags = np.array(negativeBags)
	positiveBagPairNames = np.array(positiveBagPairNames)
	negativeBagPairNames = np.array(negativeBagPairNames)
	
	#add the number of instances per bag as feature to the instances
	for bag in positiveBags:
		instCount = len(bag)

		for instance in bag:
			instance.append(instCount / positiveBags.shape[0])


	for bag in negativeBags:
		instCount = len(bag)

		for instance in bag:
			instance.append(instCount / negativeBags.shape[0])
	
	#remove instances with no variance
	posInstances = np.vstack(positiveBags)
	negInstances = np.vstack(negativeBags)

	allInstances = np.concatenate((posInstances, negInstances))

	from sklearn.feature_selection import VarianceThreshold
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



	print('Number of positive bags: ', positiveBags.shape)
	print('Number of negative bags: ', negativeBags.shape)

	print('Number of positive instances: ', len(positiveInstanceLabels))

	if positiveBags.shape[0] == 0 or negativeBags.shape[0] == 0:
		continue

	#set a random seed to always subsample the same set
	if leaveOnePatientOut == 'False':
		np.random.seed(0)
		#subsample the negative set to the same number of positives.
		negativeBagsSubsampled = np.random.choice(negativeBags, positiveBags.shape[0])

		negativeBagsSubsampleInd = np.random.choice(np.arange(negativeBags.shape[0]), positiveBags.shape[0])
		negativeBagsSubsampled = negativeBags[negativeBagsSubsampleInd]

		negativeBagPairNamesSubsampled = negativeBagPairNames[negativeBagsSubsampleInd]

		posInstances = np.vstack(positiveBags)
		negInstances = np.vstack(negativeBagsSubsampled)

		#make plot of positive instances vs negative instances
		posInstances = np.array(posInstances)
		negInstances = np.array(negInstances)
		
		print(posInstances.shape)
		print(negInstances.shape)

		#what are the differences?
		goodInstancesSum = np.sum(posInstances, axis=0)
		goodInstancesAvg = goodInstancesSum / posInstances.shape[0]
		badInstancesSum = np.sum(negInstances, axis=0)
		badInstancesAvg = badInstancesSum / negInstances.shape[0]

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
		#exit()


		bagPairLabels = np.concatenate((positiveBagPairNames, negativeBagPairNamesSubsampled))

		#save the bag pair labels for later
		np.save(finalOutDir + '/bagPairLabels_' + svType + '.npy', bagPairLabels)

		#merge the bags so that we can easily get to 1 similarity matrix and do all-to-all computations
		bags = np.concatenate((positiveBags, negativeBagsSubsampled))
		#assign bag labels
		bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBagsSubsampled.shape[0])
	else:
		#in case of leave-one-patient out, we subsample later on
		bagPairLabels = np.concatenate((positiveBagPairNames, negativeBagPairNames))

		#save the bag pair labels for later
		np.save(finalOutDir + '/bagPairLabelsNotSubsampled_' + svType + '.npy', bagPairLabels)

		#merge the bags so that we can easily get to 1 similarity matrix and do all-to-all computations
		bags = np.concatenate((positiveBags, negativeBags))
		#assign bag labels
		bagLabels = np.array([1]*positiveBags.shape[0] + [0]*negativeBags.shape[0])


	#output the bag labels which we can later read with the matrices
	np.save(finalOutDir + '/bagLabels_' + svType + '.npy', bagLabels)

	#stack the instances in the bags so that we can easily compute bag-instance distances
	instances = np.vstack(bags)

	#also output the instances for later
	np.save(finalOutDir + '/instances_' + svType + '.npy', instances)

	print(instances[:,2])

	#Make an index where we can lookup at which position the instances are in the concatenated bag array.
	reverseBagMap = dict() #lookup instance by bag index
	bagMap = dict() #lookup bag by instance index
	instanceInd = 0
	for bagInd in range(0, bags.shape[0]):
		reverseBagMap[bagInd] = []
		for instance in bags[bagInd]:
			reverseBagMap[bagInd].append(instanceInd)
			bagMap[instanceInd] = bagInd

			instanceInd += 1

	#save bagmap for later
	np.save(finalOutDir + '/bagMap_' + svType + '.npy', bagMap)

	#if we do feature selection, randomize the features here
	featureCount = instances.shape[1]

	featureStart = featureCount-1
	if featureElimination == "True":
		featureStart = 0 #set this to featureCount to run with all features. (make setting later)

	from sklearn.ensemble import RandomForestClassifier
	from sklearn.model_selection import StratifiedKFold
	from sklearn import model_selection
	from sklearn.metrics import plot_roc_curve, auc, average_precision_score
	import matplotlib.pyplot as plt
	from scipy import interp
	classifier = RandomForestClassifier(n_estimators= 200, random_state=42)


	#if featureStart is not updated, this will run once
	#otherwise it will randomize a new feature each time
	#for featureInd in range(featureStart, featureCount):
	for featureInd in range(0,featureCount):
		print('current feature: ', featureInd+1)

		if featureElimination == "True":

			#per chromosome, shuffle the features in the training set.
			#then output the original test set
			#so we have per SV type, per chromosome CV, X files for the number of features shuffled
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

			trainBags = dict()
			testBags = dict()
			trainLabels = dict()
			testLabels = dict()
			aucs = []
			performances = []
			tprs = []
			mean_fpr = np.linspace(0, 1, 100)
			fig, ax = plt.subplots()
			for chromosome in chromosomes:

				if chromosome not in positiveBagsPerChromosome:
					continue
				if chromosome not in negativeBagsPerChromosome:
					continue
				

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

				testBags[chromosome] = allTestBags
				testLabels[chromosome] = [1]*testPositiveBags.shape[0] + [0]*testSubsetNegativeBags.shape[0]
				#testLabels[chromosome] = [1]*testPositiveBags.shape[0] + [0]*testNegativeBags.shape[0]

				testPositiveInstances = np.vstack(testPositiveBags)
				testNegativeInstances = np.vstack(testSubsetNegativeBags)
				testPositiveLabels = [1]*testPositiveInstances.shape[0]
				testNegativeLabels = [0]*testNegativeInstances.shape[0]

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


				trainBags[chromosome] = np.array(trainingSet)
				trainLabels[chromosome] = trainingLabels

				trainInstances = np.vstack(trainBags[chromosome])
				
				#shuffle the training instances
				shuffledInstanceValues = trainInstances[:,featureInd]
				randomInd = np.arange(0, shuffledInstanceValues.shape[0])
				np.random.shuffle(randomInd)
				
				#we compute the similarity matrix based on the instances
				#but the instance values need to be reset every iteration
				shuffledInstances = np.zeros(trainInstances.shape)
				for col in range(0, trainInstances.shape[1]):
					if col != featureInd:
						shuffledInstances[:,col] = trainInstances[:,col]
					else:
						shuffledInstances[:,col] = trainInstances[randomInd,col]
				

				reverseBagMapOtherPatients = dict() #lookup instance by bag index
				instanceInd = 0
				for bagInd in range(0, trainBags[chromosome].shape[0]):
					reverseBagMapOtherPatients[bagInd] = []
					for instance in trainBags[chromosome][bagInd]:
						reverseBagMapOtherPatients[bagInd].append(instanceInd)
						instanceInd += 1

				similarityMatrixTrain = getSimilarityMatrix(trainBags[chromosome], shuffledInstances, reverseBagMapOtherPatients)

				#now the curent patient bags need to be to the instances of the training set
				similarityMatrixTest = getSimilarityMatrixTest(testBags[chromosome], shuffledInstances, testLabels)

				classifier.fit(similarityMatrixTrain, trainLabels[chromosome])

				trainPreds = classifier.predict(similarityMatrixTrain)
				diff = np.sum(np.abs(trainLabels[chromosome] - trainPreds)) / len(trainLabels[chromosome])
	

				preds = classifier.predict(similarityMatrixTest)

				diff = np.sum(np.abs(testLabels[chromosome] - preds)) / len(testLabels[chromosome])


				predProb = classifier.predict_proba(similarityMatrixTest)

				performances.append(classifier.score(similarityMatrixTest, testLabels[chromosome]))

				fig, ax = plt.subplots()
				viz = plot_roc_curve(classifier, similarityMatrixTest, testLabels[chromosome],
									 name='roc',
									 alpha=0.3, lw=1, ax=ax)
				interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
				interp_tpr[0] = 0.0
				tprs.append(interp_tpr)
				aucs.append(np.mean(viz.roc_auc))
				plt.close()
				
				#output these to a file
				#write these data to disk so that we can access it later on
				np.save(featureEliminationOutDir + '/' + 'similarityMatrixTrain_' + svType + '_' + chromosome + '_' + str(featureInd) + '.npy', similarityMatrixTrain)
				np.save(featureEliminationOutDir + '/' + 'similarityMatrixTest_' + svType + '_' + chromosome + '_' + str(featureInd) + '.npy', similarityMatrixTest)

				#also save the labels
				np.save(featureEliminationOutDir + '/' + 'bagLabelsTrain_' + svType + '_' + chromosome + '_' + str(featureInd) + '.npy', trainLabels)
				np.save(featureEliminationOutDir + '/' + 'bagLabelsTest_' + svType + '_' + chromosome + '_' + str(featureInd) + '.npy', testLabels)

			print(np.mean(aucs))
			plt.close()
			continue


			#randomize one feature across the bags
			#get all values of this instance
			shuffledInstanceValues = instances[:,featureInd]
			randomInd = np.arange(0, shuffledInstanceValues.shape[0])
			np.random.shuffle(randomInd)
			print(randomInd)

			#we compute the similarity matrix based on the instances
			#but the instance values need to be reset every iteration
			shuffledInstances = np.zeros(instances.shape)
			for col in range(0, instances.shape[1]):
				if col != featureInd:
					shuffledInstances[:,col] = instances[:,col]
				else:
					shuffledInstances[:,col] = instances[randomInd,col]

			print(shuffledInstances[:,featureInd])
			print(instances[:,featureInd])

			#bags are just used for shape, so no need to shuffle those
			similarityMatrix = getSimilarityMatrix(bags, shuffledInstances, reverseBagMap)
			#output this similarity matrix to a file.
			#output to a folder specific for the feature selection data
			np.save(featureEliminationOutDir + '/similarityMatrix_' + svType + '_' + str(featureInd) + '.npy', similarityMatrix)
			
		elif featureElimination == 'False' and leaveOnePatientOut == 'False':
			#Make similarity matrix
			print("generating similarity matrix")
			similarityMatrix = getSimilarityMatrix(bags, instances, reverseBagMap)
			np.save(finalOutDir + '/similarityMatrix_' + svType + '.npy', similarityMatrix)
			print(finalOutDir + '/similarityMatrix_' + svType + '.npy')

			###### again different option, per chromosome

			chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
						   'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
						   'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
						   'chr20', 'chr21'] # 'chr22'


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
			#classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
			from cleanlab.classification import LearningWithNoisyLabels

			classifier = RandomForestClassifier(n_estimators= 200, random_state=42)
			lnl = LearningWithNoisyLabels(clf=classifier)
			trainBags = dict()
			testBags = dict()
			trainLabels = dict()
			testLabels = dict()
			for chromosome in chromosomes:
				print(chromosome)

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
				
				#simple classifier
				# classifier.fit(allTrainInstances, allTrainLabels)
				# print('train: ', classifier.score(allTrainInstances, allTrainLabels))
				# print('test: ', classifier.score(allTestInstances, allTestLabels))
				# 
				# fig, ax = plt.subplots()
				# viz = plot_roc_curve(classifier, allTestInstances, allTestLabels,
				# 					 name='roc',
				# 					 alpha=0.3, lw=1, ax=ax)
				# # interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
				# # interp_tpr[0] = 0.0
				# # tprs.append(interp_tpr)
				# # aucs.append(np.mean(viz.roc_auc))
				# print('auc: ', np.mean(viz.roc_auc))
				# 
				# #print feature importances to see what the classifier is using
				# importances = classifier.feature_importances_
				# 
				# #rank these importances by score
				# indices = np.argsort(importances)[::-1]
				# 
				# topInstances = allTrainInstances[indices[1:100]]
				# 
				# sumInstances = np.sum(topInstances, axis=0)
				# avgInstances = sumInstances / topInstances.shape[1]
				# 
				# plt.bar(range(0, avgInstances.shape[0]), avgInstances)
				# plt.show()
				# exit()
				# 
				# 
				# continue
				
				

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

				#output these to a file
				#write these data to disk so that we can access it later on
				np.save(finalOutDir + '/' + 'similarityMatrixTrain_' + chromosome + '_' + svType + '.npy', similarityMatrixTrain)
				np.save(finalOutDir + '/' + 'similarityMatrixTest_' + chromosome + '_' + svType + '.npy', similarityMatrixTest)

				#also save the labels
				np.save(finalOutDir + '/' + 'bagLabelsTrain_' + chromosome + '_' + svType + '.npy', trainLabels)
				np.save(finalOutDir + '/' + 'bagLabelsTest_' + chromosome + '_' + svType + '.npy', testLabels)

				continue

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
				

				xlabels = ['loss', 'gain', 'cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
						   'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
						   'enhancer_s', 'ctcf_s', 'rnaPol_s', 'h3k9me3_s', 'h3k4me3_s', 'h3k27ac_s', 'h3k27me3_s', 'h3k4me1_s', 'h3k36me3_s', 'enhancerType', 'promoterType', 'eQTLType', 'superEnhancerType']


				#get the instances of the first 3 bags
				goodInstances = []
				goodNegInstances = []
				badInstances = []
				ind = 0
				for bag in trainBags[chromosome]:
					#if trainLabels[chromosome][ind] == 1 and trainPreds[ind] == 1:
					if trainLabels[chromosome][ind] == 1:
						#print('pred diff: ', preds[ind], testLabels[chromosome][ind])
						instances = np.vstack(bag)
						for instance in instances:
							goodInstances.append(instance)

					#elif trainLabels[chromosome][ind] == 0 and trainPreds[ind] == 0:
					elif trainLabels[chromosome][ind] == 0:
						instances = np.vstack(bag)
						for instance in instances:
							goodNegInstances.append(instance)
					#
					# elif ind == 6:
					# 	print('pred diff bad: ', preds[ind], testLabels[chromosome][ind])
					# 	instances = np.vstack(bag)
					# 	print('instances of false but pred true: ')
					# 	for instance in instances:
					# 		badInstances.append(instance)
					# 		featureInd = 0
					# 		for feature in instance:
					# 			print('feature: ', xlabels[featureInd])
					# 			print(feature)
					# 			featureInd += 1
								


							
							
					#elif ind == 0:
					#	print()
						
					ind += 1

				goodInstances = np.array(goodInstances)
				goodNegInstances = np.array(goodNegInstances)
				
				print(goodInstances.shape)
				print(goodNegInstances.shape)

				#what are the differences?
				goodInstancesSum = np.sum(goodInstances, axis=0)
				goodInstancesAvg = goodInstancesSum / goodInstances.shape[1]
				badInstancesSum = np.sum(goodNegInstances, axis=0)
				badInstancesAvg = badInstancesSum / goodNegInstances.shape[1]


				# plt.bar(range(0, len(goodInstancesAvg)), goodInstancesAvg)
				# plt.xticks(range(0, len(xlabels)), xlabels, rotation=90)
				# plt.show()
				# plt.bar(range(0, len(badInstancesAvg)), badInstancesAvg)
				# plt.xticks(range(0, len(xlabels)), xlabels, rotation=90)
				# plt.show()
				#exit()
				
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
			
			# ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
			# 	label='Chance', alpha=.8)
			#
			# mean_tpr = np.mean(tprs, axis=0)
			# mean_tpr[-1] = 1.0
			# mean_auc = auc(mean_fpr, mean_tpr)
			# std_auc = np.std(aucs)
			#
			# ax.plot(mean_fpr, mean_tpr, color='b',
			# 		label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (np.mean(aucs), np.std(aucs)),
			# 		lw=2, alpha=.8)
			#
			# std_tpr = np.std(tprs, axis=0)
			# tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
			# tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
			# ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
			# 				label=r'$\pm$ 1 std. dev.')
			#
			# ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
			# 	   title="Receiver operating characteristic: ")
			# ax.legend(loc="lower right")
			# plt.tight_layout()
			# plt.show()
			exit()
			
			print(allTrueLabels)
			print(allPreds)
			
			allTrueLabels = np.array(allTrueLabels)
			allPreds = np.array(allPreds)

			np.savetxt('true_del.txt', allTrueLabels, fmt='%s', delimiter='\t')
			np.savetxt('pred_del.txt', allPreds, fmt='%s', delimiter='\t')

			#overall scoring
			diff = np.sum(np.abs(allTrueLabels - allPreds)) / len(allTrueLabels)
			print('overall diff: ', diff)
			#auc
			print(auc(allTrueLabels, allPreds))






			exit()

			#################### testing

			#print(positiveBagPairNames)
			#print(negativeBagPairNames)
			#exit()

			#test with folds

			#divide into X bags, regardless of patients
			foldSize = 10
			import math
			bagsPerFold = math.ceil((bags.shape[0] / foldSize) / 2)

			#in each fold, randomly sample positive bags and negative bags of same size
			trainBags = dict()
			testBags = dict()

			trainLabels = dict()
			testLabels = dict()

			testPairLabels = dict()
			#set random bags to use for each fold
			random.seed(785)
			randInd = random.sample(range(0, positiveBags.shape[0]), positiveBags.shape[0])
			randIndNegative = random.sample(range(0, negativeBags.shape[0]), positiveBags.shape[0])
			currentInd = 0
			currentUntil = currentInd + bagsPerFold

			print(positiveBags.shape[0])
			print(negativeBags.shape[0])
			for foldInd in range(0, foldSize):

				#randomly sample x positive and negative bags
				randomPositive = positiveBags[randInd[currentInd:currentUntil]]
				randomNegative = negativeBags[randInd[currentInd:currentUntil]]

				#and labels
				positiveLabels = [1]*randomPositive.shape[0]
				negativeLabels = [0]*randomNegative.shape[0]

				testBags[foldInd] = np.concatenate((randomPositive, randomNegative))
				testLabels[foldInd] = positiveLabels + negativeLabels

				#also get the pair labels
				testPairLabels[foldInd] = positiveBagPairNames[randInd[currentInd:currentUntil]]

				#then the training set will be all other bags

				otherPosInd = []
				for ind in randInd:
					if ind not in randInd[currentInd:currentUntil]:
						otherPosInd.append(ind)
				otherNegInd = []
				for ind in randInd:
					if ind not in randInd[currentInd:currentUntil]:
						otherNegInd.append(ind)

				positiveTrain = positiveBags[otherPosInd]
				negativeTrain = negativeBags[otherPosInd]

				trainPairLabels = positiveBagPairNames[otherPosInd]

				print(foldInd)
				splitTrainPairs = dict()
				for pair in trainPairLabels:
					splitLabel = pair.split('_')
					splitTrainPairs[splitLabel[7] + '_' + splitLabel[0]] = pair
					
				for pair in testPairLabels[foldInd]:
					
					splitPair = pair.split('_')
					if splitPair[7] + '_' + splitPair[0] in splitTrainPairs:
						print('found pair: ')
						print('test: ', pair)
						print('train: ', splitTrainPairs[splitPair[7] + '_' + splitPair[0]])
					
				
				trainBags[foldInd] = np.concatenate((positiveTrain, negativeTrain))
				trainLabels[foldInd] = [1]*len(otherPosInd) + [0]*len(otherNegInd)

				currentInd += bagsPerFold
				if currentUntil + bagsPerFold > positiveBags.shape[0]:
					currentUntil = positiveBags.shape[0]
				else:
					currentUntil += bagsPerFold
			exit()
			#print(trainBags[0])
			#print(testBags[0])
			#print(trainLabels[0])
			#print(testLabels[0])

			#do a quick classifier test here to see if this works or not.

			#make sim matrix per fold
			from sklearn.ensemble import RandomForestClassifier
			from sklearn.model_selection import StratifiedKFold
			from sklearn import model_selection
			from sklearn.metrics import plot_roc_curve, auc
			import matplotlib.pyplot as plt
			from scipy import interp
			#train simple rf, check performance
			aucs = []
			performances = []
			posPerformances = []
			for foldInd in range(0, foldSize):
				print(foldInd)

				#get instances
				trainInstances = np.vstack(trainBags[foldInd])

				#this needs a bag map, which is changed each time we make subsets.
				reverseBagMapOtherPatients = dict() #lookup instance by bag index
				instanceInd = 0
				for bagInd in range(0, trainBags[foldInd].shape[0]):
					reverseBagMapOtherPatients[bagInd] = []
					for instance in trainBags[foldInd][bagInd]:
						reverseBagMapOtherPatients[bagInd].append(instanceInd)
						instanceInd += 1

				#collect all this information as total bags/labels

				similarityMatrixTrain = getSimilarityMatrix(trainBags[foldInd], trainInstances, reverseBagMapOtherPatients)
				print(similarityMatrixTrain.shape)
				#now the curent patient bags need to be to the instances of the training set
				similarityMatrixTest = getSimilarityMatrixTest(testBags[foldInd], trainInstances, testLabels)
				print(similarityMatrixTest.shape)
				
				np.save(finalOutDir + '/similarityMatrixTrain_' + svType + '_' + str(foldInd) + '.npy', similarityMatrixTrain)
				np.save(finalOutDir + '/similarityMatrixTest_' + svType + '_' + str(foldInd) + '.npy', similarityMatrixTest)
				np.save(finalOutDir + '/bagLabelsTrain_' + svType + '_' + str(foldInd) + '.npy', trainLabels[foldInd])
				np.save(finalOutDir + '/bagLabelsTest_' + svType + '_' + str(foldInd) + '.npy', testLabels[foldInd])
				print(finalOutDir + '/similarityMatrix_' + svType + '.npy')
				
				#also about the pair labels per fold, we only need this for the test set
				np.save(finalOutDir + '/bagPairLabelsTest_' + svType + '_' + str(foldInd) + '.npy', testPairLabels[foldInd])


			# 	classifier = RandomForestClassifier(n_estimators= 100)
			# 	#then train the classifier
			# 	classifier.fit(similarityMatrixTrain, trainLabels[foldInd])
			# 	print(classifier.predict(similarityMatrixTest))
			#
			# 	preds = classifier.predict(similarityMatrixTrain)
			# 	diff = np.sum(np.abs(trainLabels[foldInd] - preds)) / len(trainLabels[foldInd])
			# 	print('train diff: ', diff)
			#
			# 	preds = classifier.predict(similarityMatrixTest)
			# 	diff = np.sum(np.abs(testLabels[foldInd] - preds)) / len(testLabels[foldInd])
			# 	print('test diff: ', diff)
			#
			#
			#
			# 	print('train: ', classifier.score(similarityMatrixTrain, trainLabels[foldInd]))
			# 	print('test: ', classifier.score(similarityMatrixTest, testLabels[foldInd]))
			# 	performances.append(classifier.score(similarityMatrixTest, testLabels[foldInd]))
			#
			# 	fig, ax = plt.subplots()
			# 	viz = plot_roc_curve(classifier, similarityMatrixTest, testLabels[foldInd],
			# 						 name='roc',
			# 						 alpha=0.3, lw=1, ax=ax)
			# 	aucs.append(np.mean(viz.roc_auc))
			# 	print('auc: ', np.mean(viz.roc_auc))
			# print(np.mean(performances))
			# print(np.mean(aucs))



		elif featureElimination == 'False' and leaveOnePatientOut == 'True':


			#first, get the bags and labels per patient
			perPatientPositiveBags = dict()
			for bagInd in range(0, positiveBags.shape[0]):

				#get the label of this bag
				bagPairLabel = positiveBagPairNames[bagInd]
				splitLabel = bagPairLabel.split('_')

				shortPair = splitLabel[7] + '_' + splitLabel[0]

				if shortPair in degPairs[:,0]:
					pairInfo = degPairs[degPairs[:,0] == shortPair][0]
					if float(pairInfo[5]) < 1.5 and float(pairInfo[5]) > -1.5:
						print(pairInfo)

				patientId = splitLabel[7]
				if patientId not in perPatientPositiveBags:
					perPatientPositiveBags[patientId] = dict()
					perPatientPositiveBags[patientId]['bags'] = []

				perPatientPositiveBags[patientId]['bags'].append(positiveBags[bagInd])

			perPatientNegativeBags = dict()
			for bagInd in range(0, negativeBags.shape[0]):

				#get the label of this bag
				bagPairLabel = negativeBagPairNames[bagInd]
				splitLabel = bagPairLabel.split('_')


				patientId = splitLabel[7]
				if patientId not in perPatientNegativeBags:
					perPatientNegativeBags[patientId] = dict()
					perPatientNegativeBags[patientId]['bags'] = []

				perPatientNegativeBags[patientId]['bags'].append(negativeBags[bagInd])
	
			#for each patient, randomly subsample as many negative bags as there are positives
			perPatientBags = dict()
			skippedPatients = 0
			for patient in perPatientPositiveBags:

				if patient not in perPatientNegativeBags:
					skippedPatients += 1
					continue

				if patient not in perPatientBags:
					perPatientBags[patient] = dict()
					perPatientBags[patient]['bags'] = []
					perPatientBags[patient]['labels'] = []

				patientNegativeBags = perPatientNegativeBags[patient]['bags']
				patientNegativeBags = np.array(patientNegativeBags)
				#add the same number of positives/negatives

				if len(perPatientPositiveBags[patient]['bags']) > patientNegativeBags.shape[0]:
					sampleCount = patientNegativeBags.shape[0]

					randomInd = random.sample(range(0, patientNegativeBags.shape[0]), sampleCount)

					randomNegativeBags = patientNegativeBags[randomInd]

					for bag in randomNegativeBags:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(0)

					for bag in perPatientPositiveBags[patient]['bags']:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(1)

				else:
					sampleCount = len(perPatientPositiveBags[patient]['bags'])

					randomInd = random.sample(range(0, patientNegativeBags.shape[0]), sampleCount)

					randomNegativeBags = patientNegativeBags[randomInd]

					for bag in randomNegativeBags:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(0)

					for bag in perPatientPositiveBags[patient]['bags']:
						perPatientBags[patient]['bags'].append(bag)
						perPatientBags[patient]['labels'].append(1)
				
			print(skippedPatients)

			#go through each patient, and divide into train/test
			#the training set will be a merge of the bags of all other patients.
			#then for each patient, get the train/test combination
			foldSize = 1

			import math
			folds = math.ceil(len(perPatientBags) / foldSize)
			print(folds)

			testPatients = dict()
			trainPatients = dict()
			ind = 1
			foldInd = 0
			for patient in perPatientBags:
				if foldInd not in testPatients:
					testPatients[foldInd] = []
					trainPatients[foldInd] = []

				testPatients[foldInd].append(patient)
				if ind % foldSize == 0:

					for patient2 in perPatientBags:

						if patient2 not in testPatients[foldInd]:
							trainPatients[foldInd].append(patient2)
					foldInd += 1

				ind += 1

			#add remaining patients
			if foldInd in testPatients and len(trainPatients[foldInd]) < 1:
				for patient2 in perPatientBags:
					if patient2 not in testPatients[foldInd]:
						trainPatients[foldInd].append(patient2)


			for fold in testPatients:

				testBags = []
				trainBags = []
				testLabels = []
				trainLabels = []
				for patient in perPatientBags:

					patientBags = perPatientBags[patient]['bags']
					patientLabels = perPatientBags[patient]['labels']

					if patient in testPatients[fold]:
						testBags += patientBags
						testLabels += patientLabels
					else:
						trainBags += patientBags
						trainLabels += patientLabels

				testBags = np.array(testBags)
				trainBags = np.array(trainBags)

				#get instances
				trainInstances = np.vstack(trainBags)

				#this needs a bag map, which is changed each time we make subsets.
				reverseBagMapOtherPatients = dict() #lookup instance by bag index
				instanceInd = 0
				for bagInd in range(0, trainBags.shape[0]):
					reverseBagMapOtherPatients[bagInd] = []
					for instance in trainBags[bagInd]:
						reverseBagMapOtherPatients[bagInd].append(instanceInd)
						instanceInd += 1

				#collect all this information as total bags/labels

				similarityMatrixTrain = getSimilarityMatrix(trainBags, trainInstances, reverseBagMapOtherPatients)
				#print(similarityMatrixTrain.shape)
				#now the curent patient bags need to be to the instances of the training set
				similarityMatrixTest = getSimilarityMatrixTest(testBags, trainInstances, testLabels)
				print(similarityMatrixTest.shape)

				#write these data to disk so that we can access it later on
				np.save(leaveOnePatientOutDir + '/' + 'similarityMatrixTrain_' + str(fold) + '_' + svType + '.npy', similarityMatrixTrain)
				np.save(leaveOnePatientOutDir + '/' + 'similarityMatrixTest_' + str(fold) + '_' + svType + '.npy', similarityMatrixTest)

				#also save the labels
				np.save(leaveOnePatientOutDir + '/' + 'bagLabelsTrain_' + str(fold) + '_' + svType + '.npy', trainLabels)
				np.save(leaveOnePatientOutDir + '/' + 'bagLabelsTest_' + str(fold) + '_' + svType + '.npy', testLabels)


		else:
			print('Combination of options not implemented')
			exit(1)
			