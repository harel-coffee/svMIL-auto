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
			splitLine = line.split('\t')
			if lineCount < 1:
				header = splitLine
				lineCount += 1
				continue
		
			#Obtain required data by header information
			chr1Index = header.index("chr1")
			s1Index = header.index("s1")
			e1Index = header.index("e1")
			chr2Index = header.index("chr2")
			s2Index = header.index("s2")
			e2Index = header.index("e2")
			noOfGenesInWindowIndex = header.index("noOfGenesInWindow")
			overlappingTadBoundariesIndex = header.index("overlappingTadBoundaries")
			
			#For these features we can maybe take the maximum value.
			hiCDegreeIndex = header.index("hiCDegree")
			hiCBetweennessIndex = header.index("hiCBetweenness")
			
			noOfGenesInWindow = splitLine[noOfGenesInWindowIndex]
			overlappingTadBoundaries = splitLine[overlappingTadBoundariesIndex]
			
			if len(splitLine[hiCDegreeIndex]) > 0:
				hiCDegree = max(splitLine[hiCDegreeIndex]) #take max for now
			else:
				hiCDegree = 0
			if len(splitLine[hiCBetweennessIndex]) > 0:
				hiCBetweenness = max(splitLine[hiCBetweennessIndex])
			else:
				hiCBetweenness = 0
			currentAnnotations = [int(noOfGenesInWindow), int(overlappingTadBoundaries), int(hiCDegree), int(hiCBetweenness)]
			
			annotations.append(currentAnnotations)	

	return np.matrix(annotations, dtype="int")

# Define data
patient1P = np.array([[50,150],[60, 200]])
patient2P = np.array([[60, 200],[100, 50]])
patient3P = np.array([[100, 75],[80, 80]])
patient4P = np.array([[120, 50],[90, 100]])

#Negative bags. 
patient1N = np.array([[1, 4],[0, 9]])
patient2N = np.array([[10, 5],[14, 3]])
#Add more patients later for balance

bags = np.array([patient1P, patient2P, patient1N, patient3P, patient4P, patient2N])
labels = np.array([1,1,-1, 1, 1, -1])

train_bags = bags[3:]
print "original train_bags: "
print train_bags

train_labels = labels[3:]
test_bags = bags[:3]
test_labels = labels[:3]

###### Real data test

#Read some test and training data from the real datasets (annotated, the output of main.py)
truePositives = sys.argv[1]
trueNegatives = sys.argv[2]

#Read the annotation files into a list of features
### there is still a problem with the files, the output does not look as expected for the TN

tpAnnotations = readAnnotationData(truePositives)
tnAnnotations = readAnnotationData(trueNegatives)

#Subset the data for training and testing (arbitrary)
#We need to make some bags
#There is patient information in a original file, but this is not parsed right now. In the future, each bag can be a patient, and each instance is a variant of the patient. 
#For now, I will just do that arbitrarily and see what comes out

#Lots of conversions between numpy going on here, it could probably be better

allTpBags = []
allTnBags = []

#Make 10 bags with 50 instances
totalBagNum = 10
instanceNum = 50

currentInstanceNum = 0
for bagNum in range(0, totalBagNum):
	
	tpBag = tpAnnotations[currentInstanceNum:instanceNum, :]
	allTpBags.append(tpBag)
	
	tnBag = tnAnnotations[currentInstanceNum:instanceNum, :]
	allTnBags.append(tnBag)
	
	currentInstanceNum += instanceNum
	instanceNum += instanceNum


allTpBags = np.array(allTpBags)

allTnBags = np.array(allTnBags)

#Determine training/test labels

#Is this still the good format?? It does not really look like the other one anymore, there's now types which is a bit weird
#Also I should change the naming conventions in the rest of the code to look the same as in the rest of the program

train_bags = np.array([allTpBags[5:,:], allTnBags[5:,:]])
train_bags = np.concatenate((train_bags[0], train_bags[1]), axis=0) #there must be a better way... but not looking into it now
train_labels = np.array([1,1,1,1,1,-1,-1,-1,-1,-1]) #some random labels for now
test_bags = np.array([allTpBags[:5,:], allTnBags[:5,:]])
test_bags = np.concatenate((test_bags[0], test_bags[1]), axis=0)
test_labels = train_labels

print "training data:"
print train_bags

#1. Obtain centers in each bag. Here I will use random points
centers = []
for bag in train_bags:
	
	#Number of instances in the bag
	
	instanceNum = bag.shape[0]
	randomInstanceInd = np.random.choice(range(0, instanceNum), 1)[0]
	
	randomInstance = bag[randomInstanceInd,:]
	centers.append(randomInstance)
	

#2. Compute the distance from each cluster center to other points and take the minimum
#3. Generate a similarity matrix
similarityMatrix = np.zeros([len(centers), len(train_bags)])

for centerInd in range(0, len(centers)):
	
	for bagInd in range(0, len(train_bags)):
		
		#Skip this bag if the center is in the current bag
		if centerInd == bagInd:
			continue
		
		smallestDistance = float("inf")
		for instance in train_bags[bagInd]:
			distance = np.sum(np.abs(instance - centers[centerInd]))
			
			if distance < smallestDistance:
				smallestDistance = distance
		
		similarityMatrix[centerInd, bagInd] = smallestDistance
		
				
print similarityMatrix
		
#4. Train a classifier on the similarity space


rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=10)
rfClassifier.fit(similarityMatrix, train_labels)

#5. Test the classifier on a new point (should this also be in the similarity space? Yes, right?)
print test_bags

#Convert test bags to a format that we can use
testInstances = np.vstack(test_bags)

#Do we now compute the distance to the training data to get the same number of features?


testSimilarityMatrix = np.zeros([testInstances.shape[0], len(test_bags)])
testLabels = []

#Fix the labels per instance
labelCount = 0
for bag in test_bags:
		
	for instance in bag:
		testLabels.append(test_labels[labelCount])
	labelCount += 1



for centerInd in range(0, testInstances.shape[0]):
	
	
	for bagInd in range(0, len(train_bags)):
		
		#Skip this bag if the center is in the current bag
		if centerInd == bagInd:
			continue
		
		smallestDistance = float("inf")
		for instance in test_bags[bagInd]:
			distance = np.sum(np.abs(instance - testInstances[centerInd]))
			
			if distance < smallestDistance:
				smallestDistance = distance
		
		testSimilarityMatrix[centerInd, bagInd] = smallestDistance
		

score = rfClassifier.score(testSimilarityMatrix, testLabels)
print score

print rfClassifier.predict(testSimilarityMatrix)
