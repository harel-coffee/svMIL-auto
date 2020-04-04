### test alternative to MIL

import sys
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn import model_selection
from sklearn.metrics import plot_roc_curve, auc, average_precision_score
import matplotlib.pyplot as plt
from scipy import interp
import random

positivePairs = np.loadtxt(sys.argv[1], dtype='object')
negativePairs = np.loadtxt(sys.argv[2], dtype='object')

#try simple classification on this.

#divide into chromosomes.

positivePairsPerChromosome = dict()
for pair in positivePairs:

	splitPair = pair[0].split('_')

	if splitPair[12] != 'DEL':
		continue

	if splitPair[1] not in positivePairsPerChromosome:
		positivePairsPerChromosome[splitPair[1]] = []

	positivePairsPerChromosome[splitPair[1]].append(pair[1:])

negativePairsPerChromosome = dict()
for pair in negativePairs:

	splitPair = pair[0].split('_')

	if splitPair[12] != 'DEL':
		continue

	if splitPair[1] not in negativePairsPerChromosome:
		negativePairsPerChromosome[splitPair[1]] = []

	negativePairsPerChromosome[splitPair[1]].append(pair[1:])

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
			   'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
			   'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
			   'chr20', 'chr21', 'chr22']


from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis


names = ["Nearest Neighbors", "Linear SVM", "RBF SVM", "Gaussian Process",
         "Decision Tree", "Random Forest", "Neural Net", "AdaBoost",
         "Naive Bayes", "QDA"]

classifiers = [
    KNeighborsClassifier(3),
    SVC(kernel="linear", C=0.025),
    SVC(gamma=2, C=1),
    GaussianProcessClassifier(1.0 * RBF(1.0)),
    DecisionTreeClassifier(max_depth=5),
    RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
    MLPClassifier(alpha=1, max_iter=1000),
    AdaBoostClassifier(),
    GaussianNB(),
    QuadraticDiscriminantAnalysis()]

clfInd = -1
for classifier in classifiers:
	clfInd += 1
	print(names[clfInd])
	
	#go through the chromosomes and make the positive/negative sets
	#subsample the negative set to the same size.
	aucs = []
	for chromosome in chromosomes:
	
		testPositivePairs = np.array(positivePairsPerChromosome[chromosome])
		testNegativePairs = np.array(negativePairsPerChromosome[chromosome])

		#randomly subsample
		randInd = random.sample(range(0, testNegativePairs.shape[0]), testPositivePairs.shape[0])
		testNegativeSubset = testNegativePairs[randInd]

		testSet = np.concatenate((testPositivePairs, testNegativeSubset))
		testLabels = [1]*testPositivePairs.shape[0] + [0]*testNegativeSubset.shape[0]

		trainingSet = []
		trainingLabels = []
		for chromosome2 in positivePairsPerChromosome:

			if chromosome2 == chromosome:
				continue

			chrPositivePairs = np.array(positivePairsPerChromosome[chromosome2])
			chrNegativePairs = np.array(negativePairsPerChromosome[chromosome2])

			#randomly subsample
			randInd = random.sample(range(0, chrNegativePairs.shape[0]), chrPositivePairs.shape[0])
			chrNegativeSubset = chrNegativePairs[randInd]

			#chrTrainingSet = np.concatenate((chrPositivePairs, chrNegativeSubset))

			#trainingSet.append(list(chrTrainingSet))

			for pair in chrPositivePairs:
				trainingSet.append(pair)
			for pair in chrNegativeSubset:
				trainingSet.append(pair)

			chrLabels = [1]*chrPositivePairs.shape[0] + [0]*chrNegativeSubset.shape[0]
			trainingLabels += chrLabels

		trainingSet = np.array(trainingSet)


		#test classifier for this fold
		# print(chromosome)
		#
		# print(trainingSet.shape)
		# print(len(trainingLabels))
		# print(testSet.shape)
		# print(len(testLabels))

		clf = RandomForestClassifier(n_estimators= 200, random_state=42)
		clf.fit(trainingSet, trainingLabels)
		# print('train: ', clf.score(trainingSet, trainingLabels))
		# print('test: ', clf.score(testSet, testLabels))
		#
		fig, ax = plt.subplots()
		viz = plot_roc_curve(clf, testSet, testLabels,
							 name='roc',
							 alpha=0.3, lw=1, ax=ax)

		#print('auc: ', np.mean(viz.roc_auc))
		aucs.append(np.mean(viz.roc_auc))
		plt.close()

	print(np.mean(aucs))


