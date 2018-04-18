#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"

""" Basic test for a MILES-based classifier


#Rather than following all steps in the paper, we can try a slightly different implementation

#Steps:

#1. Obtain cluster centers in each bag

#2. Compute distance from cluster centers to the instances of other bags and take the minimum

#3. Generate a similarity space

#4. Train a classifier on the similarity space (using instance labels)

#5. Classify each instance in a test bag as positive or negative

Current issues:
	- The output file TN_annotations does not look good (nothing will work right now because there are no degree annotations for the TN data)
	- Have a better way of handling features that do not exist for some variants 


This script is currently just a dump of some code to test the MILES classification, it is in no way in it's final form O-(=_= O-)

"""

### Imports ###

import sys
sys.path.insert(0, '../')

import numpy as np
from sklearn.ensemble import RandomForestClassifier
from inputDataParser import InputDataParser

### Code ###

def readAnnotationData(annotationFile):
	"""
		Just some ugly function that will read an annotation file into a numpy matrix.
		Will eventually go to its own class if MILES works out for the project
	
		I will first limit the method to values that can be used in absolute distances.
		Later this can be extended to potentially having a different distance per feature type. 
	
	"""
	
	annotations = []
	
	with open(annotationFile, "r") as f:
		lineCount = 0
		header = []
		for line in f:
			#line = line.strip() #do not do this as it removes entries that are empty if these are at the end of the line
			#print line
			
			if lineCount < 1:
				line = line.strip() #we do not want the newline in the header to be able to match on columns
				splitLine = line.split('\t')
				header = splitLine
				
				lineCount += 1
				continue
			splitLine = line.split('\t')
			#Obtain required data by header information
			chr1Index = header.index("chr1")
			s1Index = header.index("s1")
			e1Index = header.index("e1")
			chr2Index = header.index("chr2")
			s2Index = header.index("s2")
			e2Index = header.index("e2")
			identifierIndex = header.index("identifier")
			noOfGenesInWindowIndex = header.index("noOfGenesInWindow")
			overlappingTadBoundariesIndex = header.index("overlappingTadBoundaries")
			
			#For these features we can maybe take the maximum value.
			hiCDegreeIndex = header.index("hiCDegree")
			hiCBetweennessIndex = header.index("hiCBetweenness")
			
			identifier = splitLine[identifierIndex]
			
			noOfGenesInWindow = splitLine[noOfGenesInWindowIndex]
			overlappingTadBoundaries = splitLine[overlappingTadBoundariesIndex]
			
			hiCDegree = []
			if len(splitLine[hiCDegreeIndex]) > 0:
				#hiCDegree = max(splitLine[hiCDegreeIndex]) #take max for now
				for degree in splitLine[hiCDegreeIndex].split(","):
					if degree != 'NA':
						hiCDegree.append(int(degree))
					else:
						hiCDegree.append(0) #set to 0 if the annotation is NA. In this case, there is no degree so 0 is appropriate. 
				#hiCDegree = [int(degree) for degree in splitLine[hiCDegreeIndex].split(",")] #Split on comma, entry is read as a string rather than a list
				
			else:
				hiCDegree = 0
				
			hiCBetweenness = []
			if len(splitLine[hiCBetweennessIndex]) > 0:
				#hiCBetweenness = max(splitLine[hiCBetweennessIndex])
				for betweenness in splitLine[hiCBetweennessIndex].split(","):
					if betweenness != 'NA':
						hiCBetweenness.append(float(betweenness))
					else:
						hiCBetweenness.append(0)
			else:
				hiCBetweenness = 0
			
			#Currently the list-based distance function is very slow for large sets of data, so I already take the sum of the degrees and betweenness here. This also helps to deal with lists of different sizes, for which we need
			#to use something like sum anway to get to one score. 
			#currentAnnotations = [identifier, int(np.sum(np.array(hiCDegree))), int(np.sum(np.array(hiCBetweenness)))]
			currentAnnotations = [identifier, int(noOfGenesInWindow), int(overlappingTadBoundaries), int(np.sum(np.array(hiCDegree))), int(np.sum(np.array(hiCBetweenness)))]
			#currentAnnotations = [identifier, int(noOfGenesInWindow), int(overlappingTadBoundaries)]
			annotations.append(currentAnnotations)	

	return np.array(annotations, dtype="object")

# # Define data
# patient1P = np.array([[50,150],[60, 200]])
# patient2P = np.array([[60, 200],[100, 50]])
# patient3P = np.array([[100, 75],[80, 80]])
# patient4P = np.array([[120, 50],[90, 100]])
# 
# #Negative bags. 
# patient1N = np.array([[1, 4],[0, 9]])
# patient2N = np.array([[10, 5],[14, 3]])
# #Add more patients later for balance
# 
# bags = np.array([patient1P, patient2P, patient1N, patient3P, patient4P, patient2N])
# labels = np.array([1,1,-1, 1, 1, -1])
# 
# train_bags = bags[3:]
# print "original train_bags: "
# print train_bags
# 
# train_labels = labels[3:]
# test_bags = bags[:3]
# test_labels = labels[:3]

###### Real data test

#Read some test and training data from the real datasets (annotated, the output of main.py)
truePositives = sys.argv[1]
trueNegatives = sys.argv[2]

#Read the annotation files into a list of features

tpAnnotations = readAnnotationData(truePositives)
tnAnnotations = readAnnotationData(trueNegatives)

#Subset the data for training and testing (arbitrary)
#We need to make some bags
#There is patient information in a original file, but this is not parsed right now. In the future, each bag can be a patient, and each instance is a variant of the patient. 
#For now, I will just do that arbitrarily and see what comes out

#Lots of conversions between numpy going on here, it could probably be better

#Make the bags, each sample will go in one bag

#First sort the data by identifier, it is now sorted by chromosome. For the bags, the ordering of the annotations does not matter
idColumn = tpAnnotations[:,0]
sortedTpAnnotationsIndex = np.argsort(idColumn)
sortedTpAnnotations = tpAnnotations[sortedTpAnnotationsIndex,:]

idColumn = tnAnnotations[:,0]
sortedTnAnnotationsIndex = np.argsort(idColumn)
sortedTnAnnotations = tnAnnotations[sortedTnAnnotationsIndex,:]

def makeBags(annotations):
	fullBag = []
		
	previousId = annotations[0,:][0]
	seenIdentifiers = []
	bagStartInd = 0 #keep a start coordinate where the identifier for this bag starts
	currentInd = 0 #increase this each time to keep track of if the bag has ended 
	for annotation in annotations:
		#check if the identifier is still the same as the previous line, otherwise start a new bag
		if annotation[0] != previousId:
			bag = annotations[bagStartInd:currentInd, 1:] #omit the identifier
			fullBag.append(bag)
			bagStartInd = currentInd
			previousId = annotation[0] #update the identifier we are currently looking at
			seenIdentifiers.append(annotation[0])
		currentInd += 1
	
	print len(seenIdentifiers)
	print len(np.unique(np.array(seenIdentifiers)))
	
	return fullBag

allTpBags = makeBags(sortedTpAnnotations)
allTnBags = makeBags(sortedTnAnnotations)
#Convert back to numpy
allTpBags = np.array(allTpBags)
allTnBags = np.array(allTnBags)

#Determine training/test labels

#Is this still the good format?? It does not really look like the other one anymore, there's now types which is a bit weird
#Also I should change the naming conventions in the rest of the code to look the same as in the rest of the program

#Ensure that the labels are added correctly, also the separation of training/test bags needs to be made random

trainRatio = 60 #Use 60% for training, 40% for testing
#There will be a different number of patients in the positive data compared to the negative data.
positiveBagNum = allTpBags.shape[0]
negativeBagNum = allTnBags.shape[0]

#Select training data
trainingTpNum = int(round(positiveBagNum * (trainRatio / float(100.0))))
trainingTnNum = int(round(negativeBagNum * (trainRatio / float(100.0))))

#Randomly select bags that will be in the training data
trainingDataTpInd = np.random.choice(range(0, allTpBags.shape[0]), trainingTpNum, replace=False)
trainingDataTnInd = np.random.choice(range(0, allTnBags.shape[0]), trainingTnNum, replace=False)

trainingBagsTp = allTpBags[trainingDataTpInd]
trainingBagsTn = allTnBags[trainingDataTnInd]

#Use the remaining data as test data
testDataTpInd = np.setdiff1d(range(0,allTpBags.shape[0]), trainingDataTpInd)
testDataTnInd = np.setdiff1d(range(0,allTnBags.shape[0]), trainingDataTnInd)

testBagsTp = allTpBags[testDataTpInd]
testBagsTn = allTnBags[testDataTnInd]

train_bags = np.array([trainingBagsTp, trainingBagsTn])
train_bags = np.concatenate((train_bags[0], train_bags[1]), axis=0) #there must be a better way... but not looking into it now

train_labels = np.array([1]*trainingBagsTp.shape[0] + [-1]*trainingBagsTn.shape[0])

test_bags = np.array([testBagsTp, testBagsTn])
test_bags = np.concatenate((test_bags[0], test_bags[1]), axis=0)
test_labels = np.array([1]*testBagsTp.shape[0] + [-1]*testBagsTn.shape[0])

#1. Obtain centers in each bag. Here I will use random points
print "computing centers"
centers = []
for bag in train_bags:
	
	#Number of instances in the bag
	
	instanceNum = bag.shape[0]
	randomInstanceInd = np.random.choice(range(0, instanceNum), 1)[0]
	
	randomInstance = bag[randomInstanceInd,:]
	centers.append(randomInstance)
	

#2. Compute the distance from each cluster center to other points and take the minimum

#Here we first define the distance functions that we will use per feature.
#Features in order: noOfGenes pLI	RVIS	overlappingTadBoundaries	hiCDegree	hiCBetweenness
#For noOfGenes and overlappingTadBoundaries we can compute the absolute distance (sum of distances)
#For the hiCDegree and hiCBetweenness, we could in principle do the same. If there is a very high degree there, the score will be high. Even if one degree is low, the total score will remain high.
#The same can be used for the pLI and RVIS. (these are not yet in the code right now)
distanceFunctions = ["absoluteDistance", "absoluteDistance", "absoluteDistance", "absoluteDistance"]
distanceFunctions = ["absoluteDistance", "absoluteDistance"]
#distanceFunctions = ["absoluteDistance", "absoluteDistance", "listAbsoluteDistance", "listAbsoluteDistance"]

def absoluteDistance(center, instances): #compute distance simply based on integer values
	distances = np.abs(center - instances)
	#Then take the sum across all differences, and report one column
	return np.sum(distances, axis=1)
	
def listAbsoluteDistance(center, instances): #Compute distance given a list of entries, for HiC degree and HiC betweenness
	#First sum all elements (so all lists) of the same features
	#This sum is to avoid computing the difference between lists that are of different sizes (e.g. different number of degrees in window of variant)
	summedInstances = np.zeros(instances.shape)
	for row in range(0, instances.shape[0]):
		for col in range(0, instances.shape[1]):
			summedInstances[row,col] = np.sum(instances[row,col])
	
	#Do the same summing for the center as well
	
	centerSum = center
	centerSum[0] = np.sum(center[0])
	centerSum[1] = np.sum(center[1])
	
	#Then take the difference between the lists
	difference = np.abs(centerSum - summedInstances)
	
	#Return the sum of the difference
	return np.sum(difference, axis=1)
	
	

#3. Generate a similarity matrix

###THIS PART IS REALLY, REALLY SLOW, BUT THE DATASET IS ALSO LARGE SO I CONSIDER IT FOR NOW SUFFICIENT

print "generating similarity matrix"
similarityMatrix = np.zeros([len(train_bags), len(centers)])

print "size of the similarity matrix: "
print similarityMatrix.shape

print "size of the training bags: "
print train_bags.shape

print train_bags[0]
#The goal here is:
#- we want to have a matrix where the rows are the instances, and the columns are the bags.
#- To compute a distance from the bag to the instances, we want a bag center to all other instances, and take the smallest distance. 
#- Can we compute the distances quicker than now? At the moment, we first take a center for each bag, and then compute the distance to all other elements.
#- Currently, we have 350 bags. Then there will be 350x350 distances (bag to smallest instance in each bag).
#- Can we compute the distance to all instances at once and then take the minimum?

#Get the number of different disance functions we need to call, we can do it once for all features at once that we treat in the same fashion, should speed up the code
uniqueDistanceFunctions = np.unique(np.array(distanceFunctions))

for centerInd in range(0, len(centers)):
	print "center: ", centerInd
	for bagInd in range(0, len(train_bags)):
		
		#Skip this bag if the center is in the current bag
		if centerInd == bagInd:
			continue
		
		#Compute the distance from center to all instances for each set of features with the same distance function at once
		allDistances = np.zeros([train_bags[bagInd].shape[0], len(uniqueDistanceFunctions)])
		distanceFunctionInd = 0
		for distanceFunction in uniqueDistanceFunctions:
			#Check which feature indices match with the current distance function
			matchingFeatureIndices = np.where(np.array(distanceFunctions) == distanceFunction)[0]
			
			#Extract the feature columns
			#Provide to this function the center and the feature columns for all instances in the current bag
			currentDistances = locals()[distanceFunction](centers[centerInd][matchingFeatureIndices], train_bags[bagInd][:,matchingFeatureIndices])
			
			allDistances[:,distanceFunctionInd] = currentDistances
			distanceFunctionInd += 1
		
		#Here we should find the row (or the instance) that is closest to this center
		
		#Sum across the rows
		#Identify the row (instance) that has the smallest total sum
		summedDistances = np.sum(allDistances, axis=1)

		smallestDistance = np.min(summedDistances)
		
		similarityMatrix[bagInd, centerInd] = smallestDistance
	
#Also get training labels per instance (not used in the classification, only to test the performance and also to plot)
trainLabels = []

#Fix the labels per instance
labelCount = 0
for bag in train_bags:
		
	for instance in bag:
		trainLabels.append(train_labels[labelCount])
	labelCount += 1

print "similarity matrix shape: "
print similarityMatrix.shape
print "test labels shape: "
print len(train_labels)


#4. Train a classifier on the similarity space
print "training the classifier in similarity space"
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
rfClassifier.fit(similarityMatrix, train_labels) #Use the bag labels, not the instance labels

#5. Test the classifier on a new point (should this also be in the similarity space? Yes, right?)


#Convert test bags to a format that we can use
testInstances = np.vstack(test_bags)

#Do we now compute the distance to the training data to get the same number of features?


testSimilarityMatrix = np.zeros([len(test_bags), len(centers)])
testLabels = []

#Fix the labels per instance
labelCount = 0
for bag in test_bags:
		
	for instance in bag:
		testLabels.append(test_labels[labelCount])
	labelCount += 1

#Computing the distance from every instance to all other instances is computationally expensive
#we should compute the distance from a random center
#But should the distance be computed to the center points of the training data? or to centers of the test data? training, right? otherwise the space is not the same.

#Here we should not use the distance from center to bags, but from center to each instance. Otherwise we cannot do instance-based classification!

print "computing test data similarity matrix"

for centerInd in range(0, len(centers)):
	print "center: ", centerInd
	
	for bagInd in range(0, len(test_bags)):
		
		#Skip this bag if the center is in the current bag
		if centerInd == bagInd:
			continue
		
		allDistances = np.zeros([test_bags[bagInd].shape[0], len(uniqueDistanceFunctions)])
		distanceFunctionInd = 0
		for distanceFunction in uniqueDistanceFunctions:
			#Check which feature indices match with the current distance function
			matchingFeatureIndices = np.where(np.array(distanceFunctions) == distanceFunction)[0]
			
			#Extract the feature columns
			#Provide to this function the center and the feature columns for all instances in the current bag
			currentDistances = locals()[distanceFunction](centers[centerInd][matchingFeatureIndices], test_bags[bagInd][:,matchingFeatureIndices])
			
			allDistances[:,distanceFunctionInd] = currentDistances
			distanceFunctionInd += 1
		
		#Here we should find the row (or the instance) that is closest to this center
		
		#Sum across the rows
		#Identify the row (instance) that has the smallest total sum
		summedDistances = np.sum(allDistances, axis=1)

		smallestDistance = np.min(summedDistances)
		
		testSimilarityMatrix[bagInd, centerInd] = smallestDistance

print "similarity matrix shape: "
print testSimilarityMatrix.shape
print "test labels shape: "
print len(test_labels)
		
print "scoring classifier"
score = rfClassifier.score(testSimilarityMatrix, test_labels) #First try with classifying bags, later we can do instance-based classification.
print score

print rfClassifier.predict(testSimilarityMatrix)
# 
# #Make a plot of the classifier
# from matplotlib.colors import ListedColormap
# import matplotlib.pyplot as plt
# def plotClassificationResult(allData, dataSubset, labels, clf):
# 
#     cmap_light = ListedColormap(['#FFAAAA', '#AAFFAA', '#AAAAFF'])
#     cmap_bold = ListedColormap(['#FF0000', '#00FF00', '#0000FF'])
#     X = allData
#     h = 0.1 #fine size of mesh
#     x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
#     y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
#     xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
#                              np.arange(y_min, y_max, h))
#     Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
# 
#     # Put the result into a color plot
#     Z = Z.reshape(xx.shape)
#     plt.figure()
#     plt.pcolormesh(xx, yy, Z, cmap=cmap_light)
# 
#     #Plot the training data with decision boundary
#     plt.scatter(dataSubset[:,0], dataSubset[:,1], c=labels, cmap=cmap_bold)
#     plt.show()
# 	
# #We cannot plot the bags using this function, only the instance labels.
# #There is a way to do this probably, but I'll leave it for now. 
# trainingSubset = np.vstack(train_bags)
# 
# plotClassificationResult(trainingSubset, trainingSubset, trainLabels, rfClassifier)

