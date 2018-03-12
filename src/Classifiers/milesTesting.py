

""" Basic test for a MILES-based classifier


#Rather than following all steps in the paper, we can try a slightly different implementation

#Steps:

#1. Obtain cluster centers in each bag

#2. Compute distance from cluster centers to the instances of other bags and take the minimum

#3. Generate a similarity space

#4. Train a classifier on the similarity space (using instance labels)

#5. Classify each instance in a test bag as positive or negative

"""

import numpy as np
from sklearn.ensemble import RandomForestClassifier


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

#Read some test and training data from the real datasets




# Spilt dataset arbitrarily to train/test sets
train_bags = bags[3:]
train_labels = labels[3:]
test_bags = bags[:3]
test_labels = labels[:3]

#1. Obtain centers in each bag. Here I will use random points
centers = []
for bag in train_bags:
	
	#Number of instances in the bag
	
	instanceNum = bag.shape[0]
	randomInstanceInd = np.random.choice(range(0, instanceNum), 1)[0]
	
	randomInstance = bag[randomInstanceInd,:]
	centers.append(randomInstance)
	
	
print centers
	
#2. Compute the distance from each cluster center to other points and take the minimum
#3. Generate a similarity matrix
similarityMatrix = np.zeros([len(centers), len(bags)])

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


testSimilarityMatrix = np.zeros([testInstances.shape[0], len(bags)])
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






