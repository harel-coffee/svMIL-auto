"""
	Playing around with regression.
	We have the rule-based pairs with gains and losses, and TAD-based pairs, for which we have the z-scores of z-scores.
	The question is now: Can we build a model that identifies the most interesting rule-based pairs, where interesting is defined as having a more divergent z-score from 0?
	What are then the features that the model selects? Can we say anything about what using those features means for our data? 

"""

import numpy as np

#1. Load the rule-based pairs

rulePairs = np.loadtxt('Output/RankedGenes/0/BRCA/nonCoding_geneSVPairs.txt_', dtype='object')

#2. Load the tad-based pairs and z-scores of z-scores
tadPairs = np.loadtxt('tadBasedRanks.rnk', dtype='object')
print(tadPairs)

features = []
#get the same order for both data
labels = []
for pair in rulePairs:
	
	splitPair = pair[0].split('_')
	gene = splitPair[0]
	patient = splitPair[7]
	newPair = patient + '_' + gene
	
	if newPair not in tadPairs[:,0]:
		continue
	
	#get the rank
	rank = tadPairs[tadPairs[:,0] == newPair][0][1]
	labels.append(float(rank))
	features.append([float(i) for i in pair[1:]])
	
print(labels)

#3. Train a regressor model

from sklearn import linear_model
from sklearn.model_selection import train_test_split

#split dataset
X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.33, random_state=42)

clf= linear_model.LassoLars()

clf.fit(X_train, y_train)
print(clf.predict(X_test),'\n')

