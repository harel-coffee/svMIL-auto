"""
	1. Read the gene scores
	2. Read the DEGs based on patients with and patients without SVs
	3. Try a random forest classifier
"""

import sys
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score

scoreFile = sys.argv[1]
degFile = sys.argv[2]

geneScores = np.loadtxt(scoreFile, dtype="object")

#Filter the gene scores to a feature matrix
featureMatrix = []
for gene in geneScores:
	
	samples = gene[len(gene)-1]
	splitSamples = samples.split(",")
	sampleNum = len(splitSamples)
	
	geneFeatures = list(gene[4:len(gene)-3]) #everything minus the samples, and the total score
	#geneFeatures.append(sampleNum)
	featureMatrix.append(geneFeatures)
	
featureMatrix = np.array(featureMatrix, dtype="float")		
print featureMatrix

#Get the labels based on the differential expression
#degs = np.loadtxt(degFile, dtype="object")

#Get the cosmic genes
cosmicGenes = []
with open(degFile, 'rb') as f:
	lineCount = 0
	
	for line in f:
		line = line.strip()
		if lineCount == 0:
			lineCount += 1
			continue
		
		splitLine = line.split("\t")
		
		geneName = splitLine[0]
		cosmicGenes.append(geneName)
print cosmicGenes
# 
# geneLabels = []
# for gene in geneScores:
# 	if gene[0] in degs[:,0]:
# 		geneLabels.append(1)
# 	else:
# 		geneLabels.append(0)

#For COSMIC
geneLabels = []
for gene in geneScores:
	if gene[0] in cosmicGenes:
		geneLabels.append(1)
	else:
		geneLabels.append(0)

		
#print geneLabels		


#Try a random forest
X_train, X_test, y_train, y_test = train_test_split(featureMatrix, geneLabels, test_size=.4, random_state=42) #60%/40%
print X_train.shape
print X_test.shape

print np.where(np.array(y_train) == 1)[0].shape
print np.where(np.array(y_test) == 1)[0].shape



clf = RandomForestClassifier()

clf.fit(X_train, y_train)
importances = clf.feature_importances_
std = np.std([tree.feature_importances_ for tree in clf.estimators_],
             axis=0)
indices = np.argsort(importances)[::-1]

# Print the feature ranking
import matplotlib.pyplot as plt
print("Feature ranking:")

for f in range(X_train.shape[1]):
    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

# Plot the feature importances of the forest
plt.figure()
plt.title("Feature importances")
plt.bar(range(X_train.shape[1]), importances[indices],
       color="r", yerr=std[indices], align="center")
plt.xticks(range(X_train.shape[1]), indices)
plt.xlim([-1, X_train.shape[1]])
plt.show()
exit()
score = clf.score(X_test, y_test)
print("Classification score: ", score)

#With cross validation 10 fold
clf = RandomForestClassifier()
folds = 10 #change this value to use a different number of folds
scores = cross_val_score(clf, featureMatrix, geneLabels, cv=folds)
print('classification score per fold: ', scores) #show accuracy of each fold
print('total score: ', sum(scores) / len(scores))


