import sys
import numpy as np


somaticScores = np.loadtxt(sys.argv[1], dtype="object")
germlineScores = np.loadtxt(sys.argv[2], dtype="object")

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

print "computing centers"
centers = []
for bag in bags:
	
	#Number of instances in the bag
	bag = np.array(bag)
	instanceNum = bag.shape[0]
	randomInstanceInd = np.random.choice(range(0, instanceNum), 1)[0]
	
	randomInstance = bag[randomInstanceInd,:]
	centers.append(randomInstance)

print len(centers)

print "generating similarity matrix"
similarityMatrix = np.zeros([len(bags), len(centers)])

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
for bagInd in range(0, len(bags)):
	for gene in bags[bagInd]:
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


#4. Train a classifier on the similarity space
from sklearn.ensemble import RandomForestClassifier
print "training the classifier in similarity space"
rfClassifier = RandomForestClassifier(max_depth=5, n_estimators=2)
rfClassifier.fit(similarityMatrix, trainLabels) #Use the bag labels, not the instance labels

predictions = rfClassifier.predict(testSimilarityMatrix)
print predictions
print np.average(labels == np.sign(predictions))
positiveInd = np.where(np.array(predictions) == 1)[0]
print positiveInd


