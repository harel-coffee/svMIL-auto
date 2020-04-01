"""

	Here the MIL classifier performance will be tested on the input similarity matrices
	This can work with a single similarity matrix, but can also be run multiple times
	in the feature elimination way.

"""

import sys
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn import model_selection
from sklearn.metrics import plot_roc_curve, auc
import matplotlib.pyplot as plt
from scipy import interp
from random import shuffle
import os
import os.path
import glob
import re

import matplotlib
matplotlib.use('Agg')

featureElimination = sys.argv[2]
leaveOnePatientOut = sys.argv[3]
svTypes = ['DEL', 'DUP', 'INV', 'ITX']
svTypes = ['INV', 'ITX']
outDir = sys.argv[1]

def cvClassification(similarityMatrix, bagLabels, clf, svType, title, plot, plotOutputFile, outputFile):

	#get the kfold model
	kfold = model_selection.StratifiedKFold(n_splits=10, shuffle=True, random_state=10)

	#make the ROC curve and compute the AUC
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)

	fig, ax = plt.subplots()
	importances = []
	for i, (train, test) in enumerate(kfold.split(similarityMatrix, bagLabels)):
		clf.fit(similarityMatrix[train], bagLabels[train])
		viz = plot_roc_curve(clf, similarityMatrix[test], bagLabels[test],
							 name='ROC fold {}'.format(i),
							 alpha=0.3, lw=1, ax=ax)
		interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(viz.roc_auc)
		importances.append(clf.feature_importances_)

	#write this to a file rather than to screen
	score = str(np.mean(aucs)) + '\t' + str(np.std(aucs)) + '\n'

	with open(outputFile, 'a') as outF:
		outF.write(score)

	if plot == True:

		ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
				label='Chance', alpha=.8)

		mean_tpr = np.mean(tprs, axis=0)
		mean_tpr[-1] = 1.0
		mean_auc = auc(mean_fpr, mean_tpr)
		std_auc = np.std(aucs)

		ax.plot(mean_fpr, mean_tpr, color='b',
				label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (np.mean(aucs), np.std(aucs)),
				lw=2, alpha=.8)

		std_tpr = np.std(tprs, axis=0)
		tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
		tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
		ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
						label=r'$\pm$ 1 std. dev.')

		ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
			   title="Receiver operating characteristic: " + title)
		ax.legend(loc="lower right")
		plt.tight_layout()
		plt.savefig(plotOutputFile)

#function to classify on a case where 1 patient is left out at a time
#the similarity matrices are pre-made, with each time 1 patient being left out.
#so we can just load these in, and then train our standard classifier on that.
#report on the performance of each patient.
def leaveOnePatientOutCV(leaveOneOutDataFolder, classifier):

	#first get the names of all patients
	allFiles = glob.glob(leaveOneOutDataFolder + '*_' + svType + '.npy')

	patientFiles = dict()
	for dataFile in allFiles:

		#get the patient ID
		splitFileId = dataFile.split('_')
		patientId = splitFileId[len(splitFileId)-2]

		if patientId not in patientFiles:
			patientFiles[patientId] = []
		patientFiles[patientId].append(dataFile)


	#for each patient, get the train/test combination, and run the classifier
	performances = []
	aucs = []
	for patient in patientFiles:
		print(patient)

		for dataFile in patientFiles[patient]:

			if re.search('similarityMatrixTrain', dataFile):
				similarityMatrixTrain = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('similarityMatrixTest', dataFile):
				similarityMatrixTest = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('bagLabelsTrain', dataFile):
				bagLabelsTrain = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('bagLabelsTest', dataFile):
				bagLabelsTest = np.load(dataFile, encoding='latin1', allow_pickle=True)

		print(bagLabelsTest)

		#then train the classifier
		classifier.fit(similarityMatrixTrain, bagLabelsTrain)
		print('train: ', classifier.score(similarityMatrixTrain, bagLabelsTrain))
		print('test: ', classifier.score(similarityMatrixTest, bagLabelsTest))
		performances.append(classifier.score(similarityMatrixTest, bagLabelsTest))

		fig, ax = plt.subplots()
		viz = plot_roc_curve(classifier, similarityMatrixTest, bagLabelsTest,
							 name='roc',
							 alpha=0.3, lw=1, ax=ax)
		aucs.append(np.mean(viz.roc_auc))
		print('auc: ', np.mean(viz.roc_auc))
		
	print(aucs)
	print(np.mean(aucs))
	print(performances)
	print(np.mean(performances))

	#output these values to a file
	finalOutDir = outDir + '/multipleInstanceLearning/leaveOnePatientOutCV/'
	if not os.path.exists(finalOutDir):
		os.makedirs(finalOutDir)

	outFile = finalOutDir + '/leaveOnePatientOutCV_' + svType + '.txt'

	strAucs = [str(i) for i in aucs]
	strAuc = '\t'.join(strAucs)
	strAccs = [str(i) for i in performances]
	strAcc = '\t'.join(strAccs)
	with open(outFile, 'a') as outF:
		outF.write(strAuc)
		outF.write('\n')
		outF.write(str(np.mean(aucs)))
		outF.write('\n')
		outF.write(strAcc)
		outF.write('\n')
		outF.write(str(np.mean(performances)))
		outF.write('\n')

	return 0

for svType in svTypes:

	if svType == 'DEL':
		classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'deletions'
	elif svType == 'DUP':
		#classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=2, min_samples_leaf=2, max_features='sqrt', max_depth=110, bootstrap=False)
		classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		#classifier = RandomForestClassifier(n_estimators= 1200, min_samples_split=2, min_samples_leaf=4, max_features='sqrt', max_depth=70, bootstrap=False)
		title = 'duplications'
	elif svType == 'INV':
		classifier = RandomForestClassifier(n_estimators= 200, min_samples_split=5, min_samples_leaf=4, max_features='auto', max_depth=10, bootstrap=True)
		title = 'inversions'
	elif svType == 'ITX':
		classifier = RandomForestClassifier(n_estimators= 1000, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'translocations'
	else:
		classifier = RandomForestClassifier(n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'All SV types'

	#obtain the right similarity matrix and bag labels for this SV type
	dataPath = outDir + '/multipleInstanceLearning/similarityMatrices/'

	if os.path.isfile(dataPath + '/similarityMatrix_' + svType + '.npy') == False:
		continue

	similarityMatrix = np.load(dataPath + '/similarityMatrix_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	bagLabels = np.load(dataPath + '/bagLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)

	plot = False
	#don't make the plots for each feature to eliminate
	if featureElimination == "False" and leaveOnePatientOut == 'False':
		plot = True
	
		plotOutputPath = outDir + '/multipleInstanceLearning/rocCurves/'
		if not os.path.exists(plotOutputPath):
			os.makedirs(plotOutputPath)

		plotOutputFile = plotOutputPath + '/rocCurve_' + svType + '.svg'
		outputFile = outDir + '/multipleInstanceLearning/performance_' + svType + '.txt'
		cvClassification(similarityMatrix, bagLabels, classifier, svType, title, plot, plotOutputFile, outputFile)

		#repeat, but then with random labels.
		#mind here, if multiple iterations, the baglabels are permanently shuffled!

		shuffle(bagLabels)

		plotOutputFile = plotOutputPath + '/rocCurve_' + svType + '_randomLabels.svg'
		outputFile = outDir + '/multipleInstanceLearning/performance_' + svType + '_randomBagLabels.txt'
		cvClassification(similarityMatrix, bagLabels, classifier, svType, title, plot, plotOutputFile, outputFile)
	elif featureElimination == 'True' and leaveOnePatientOut == 'False':
		
		#do the feature elimination here. Go through each similarity matrix with eliminated feature,
		#and obtain the classification scores.
		#the plot should not be made, and the AUCs should be written to a file
		
		#get each similarity matrix
		
		featureEliminationDataFolder = outDir + '/multipleInstanceLearning/similarityMatrices/featureSelection/'
		featureEliminationFiles = glob.glob(featureEliminationDataFolder + '*_' + svType + '_*')

		#file to write to aucs to	
		outputFile = outDir + '/multipleInstanceLearning/featureEliminationResults_' + svType + '.txt'
		for fileInd in range(0, len(featureEliminationFiles)):

			#get the right similarity matrix for this shuffled feature
			similarityMatrix = np.load(featureEliminationDataFolder + '/similarityMatrix_' + svType + '_' + str(fileInd) + '.npy')
			cvClassification(similarityMatrix, bagLabels, classifier, svType, title, plot, '', outputFile)
	elif featureElimination == 'False' and leaveOnePatientOut == 'True':
		
		leaveOneOutDataFolder = outDir + '/multipleInstanceLearning/similarityMatrices/leaveOnePatientOut/'
		
		leaveOnePatientOutCV(leaveOneOutDataFolder, classifier)

		
	else:
		
		print('Combination of options not implemented')
		exit(1)
		
		
