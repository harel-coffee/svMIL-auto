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
leaveOnePatientOut = sys.argv[3] #1 patient at a time in the test set
leaveOneChromosomeOut = sys.argv[4] #1 chromosome at a time in the test set
leaveBagsOut = sys.argv[5] #random bags in each CV fold

svTypes = ['DEL', 'DUP', 'INV', 'ITX']
outDir = sys.argv[1]

def leaveOneChromosomeOutCV(leaveChromosomeOutDataFolder, classifier, svType, plotOutputFile, title):

	#first get the names of all patients
	allFiles = glob.glob(leaveChromosomeOutDataFolder + '*_' + svType + '.npy')

	chromosomeFiles = dict()
	for dataFile in allFiles:

		#get the chromosome name
		splitFileId = dataFile.split('_')
		chrId = splitFileId[len(splitFileId)-2]

		if chrId not in chromosomeFiles:
			chromosomeFiles[chrId] = []
		chromosomeFiles[chrId].append(dataFile)


	#for each patient, get the train/test combination, and run the classifier
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)
	fig, ax = plt.subplots()
	ind = 0
	for chromosome in chromosomeFiles:

		for dataFile in chromosomeFiles[chromosome]:

			if re.search('similarityMatrixTrain', dataFile):
				similarityMatrixTrain = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('similarityMatrixTest', dataFile):
				similarityMatrixTest = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('bagLabelsTrain', dataFile):
				bagLabelsTrain = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('bagLabelsTest', dataFile):
				bagLabelsTest = np.load(dataFile, encoding='latin1', allow_pickle=True)

		#then train the classifier
		classifier.fit(similarityMatrixTrain, bagLabelsTrain)

		#output the predictions to a file for the COSMIC analysis
		preds = classifier.predict(similarityMatrixTest)
	
		viz = plot_roc_curve(classifier, similarityMatrixTest, bagLabelsTest,
							 name='roc',
							 alpha=0.3, lw=1, ax=ax)
		interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(np.mean(viz.roc_auc))
		print('auc: ', np.mean(viz.roc_auc))

	print(np.mean(aucs))
	
	#make the CV plot, just plot the mean
	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
				label='Chance', alpha=.8)

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)

	ax.plot(mean_fpr, mean_tpr, color='b',
			label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (np.mean(aucs), np.std(aucs)),
			lw=2, alpha=.8)

	# std_tpr = np.std(tprs, axis=0)
	# tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	# tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	# ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
	# 				label=r'$\pm$ 1 std. dev.')

	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
		   title="Leave-one-patient-out CV: " + title)
	ax.legend(loc="lower right")
	plt.tight_layout()
	plt.savefig(plotOutputFile)


def bagsCVClassification(dataPath, clf, svType, title, plot, plotOutputFile, outputFile, shuffleLabels):

	#load the similarity matrices of each fold
	folds = 10
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)
	fig, ax = plt.subplots()
	for fold in range(0, folds):

		#match files on fold
		allFiles = glob.glob(dataPath + '*_' + svType + '_' + str(fold) + '.npy')

		for dataFile in allFiles:

			if re.search('similarityMatrixTrain', dataFile):
				similarityMatrixTrain = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('similarityMatrixTest', dataFile):
				similarityMatrixTest = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('bagLabelsTrain', dataFile):
				bagLabelsTrain = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('bagLabelsTest', dataFile):
				bagLabelsTest = np.load(dataFile, encoding='latin1', allow_pickle=True)

		if len(bagLabelsTest) < 1 or len(bagLabelsTrain) < 1:
			continue #for cases with not enough data to fill 10 folds.

		if shuffleLabels == True:
			shuffle(bagLabelsTrain)
			shuffle(bagLabelsTest)

		clf.fit(similarityMatrixTrain, bagLabelsTrain)

		viz = plot_roc_curve(clf, similarityMatrixTest, bagLabelsTest,
							 name='ROC fold {}'.format(fold),
							 alpha=0.3, lw=1, ax=ax)

		interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(viz.roc_auc)

	print(np.mean(aucs))

	if plot == True:

		fig, ax = plt.subplots()
		ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
				label='Chance', alpha=.8)

		mean_tpr = np.mean(tprs, axis=0)
		mean_tpr[-1] = 1.0
		mean_auc = auc(mean_fpr, mean_tpr)
		std_auc = np.std(aucs)

		ax.plot(mean_fpr, mean_tpr, color='b',
				label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (np.mean(aucs), np.std(aucs)),
				lw=2, alpha=.8)

		# std_tpr = np.std(tprs, axis=0)
		# tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
		# tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
		# ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
		# 				label=r'$\pm$ 1 std. dev.')

		ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
			   title="Leave-random-bags-out CV: " + title)
		ax.legend(loc="lower right")
		plt.tight_layout()
		plt.savefig(plotOutputFile)

#function to classify on a case where 1 patient is left out at a time
#the similarity matrices are pre-made, with each time 1 patient being left out.
#so we can just load these in, and then train our standard classifier on that.
#report on the performance of each patient.
def leaveOnePatientOutCV(leaveOneOutDataFolder, classifier, svType, plotOutputFile, title):

	#first get the names of all patients
	allFiles = glob.glob(leaveOneOutDataFolder + '*_[0-9]*' + svType + '.npy')

	patientFiles = dict()
	for dataFile in allFiles:

		#get the patient ID
		splitFileId = dataFile.split('_')
		patientId = splitFileId[len(splitFileId)-2]

		if patientId not in patientFiles:
			patientFiles[patientId] = []
		patientFiles[patientId].append(dataFile)


	#for each patient, get the train/test combination, and run the classifier
	tprs = []
	aucs = []
	mean_fpr = np.linspace(0, 1, 100)
	fig, ax = plt.subplots()
	ind = 0
	predictions = dict()
	for patient in patientFiles:

		for dataFile in patientFiles[patient]:

			if re.search('similarityMatrixTrain', dataFile):
				similarityMatrixTrain = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('similarityMatrixTest', dataFile):
				similarityMatrixTest = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('bagLabelsTrain', dataFile):
				bagLabelsTrain = np.load(dataFile, encoding='latin1', allow_pickle=True)
			if re.search('bagLabelsTest', dataFile):
				bagLabelsTest = np.load(dataFile, encoding='latin1', allow_pickle=True)

		#then train the classifier
		classifier.fit(similarityMatrixTrain, bagLabelsTrain)

		#output the predictions to a file for the COSMIC analysis
		preds = classifier.predict(similarityMatrixTest)
		predictions[patient] = preds
		

		viz = plot_roc_curve(classifier, similarityMatrixTest, bagLabelsTest,
							 name='roc',
							 alpha=0.3, lw=1, ax=ax)
		interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(np.mean(viz.roc_auc))
		print('auc: ', np.mean(viz.roc_auc))

	print(np.mean(aucs))
	
	#make the CV plot, just plot the mean
	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
				label='Chance', alpha=.8)

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)

	ax.plot(mean_fpr, mean_tpr, color='b',
			label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (np.mean(aucs), np.std(aucs)),
			lw=2, alpha=.8)

	# std_tpr = np.std(tprs, axis=0)
	# tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	# tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	# ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
	# 				label=r'$\pm$ 1 std. dev.')

	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
		   title="Leave-one-patient-out CV: " + title)
	ax.legend(loc="lower right")
	plt.tight_layout()
	plt.savefig(plotOutputFile)
	

	#output these values to a file
	finalOutDir = outDir + '/multipleInstanceLearning/leaveOnePatientOutCV/'
	if not os.path.exists(finalOutDir):
		os.makedirs(finalOutDir)

	outFile = finalOutDir + '/leaveOnePatientOutCV_' + svType + '.txt'

	strAucs = [str(i) for i in aucs]
	strAuc = '\t'.join(strAucs)

	with open(outFile, 'a') as outF:

		for patient in predictions:
			outF.write(patient + '\t' + '\t'.join([str(i) for i in predictions[patient]]) + "\n")


for svType in svTypes:

	if svType == 'DEL':
		classifier = RandomForestClassifier(random_state=785, n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'deletions'
	elif svType == 'DUP':
		classifier = RandomForestClassifier(random_state=785, n_estimators= 600, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'duplications'
	elif svType == 'INV':
		classifier = RandomForestClassifier(random_state=785, n_estimators= 200, min_samples_split=5, min_samples_leaf=4, max_features='auto', max_depth=10, bootstrap=True)
		title = 'inversions'
	elif svType == 'ITX':
		classifier = RandomForestClassifier(random_state=785, n_estimators= 1000, min_samples_split=5, min_samples_leaf=1, max_features='auto', max_depth=80, bootstrap=True)
		title = 'translocations'
	else:
		print('Please provide a SV type')
		exit(1)

	#obtain the right similarity matrix and bag labels for this SV type


	#if os.path.isfile(dataPath + '/similarityMatrix_' + svType + '.npy') == False:
	#	continue

	#similarityMatrix = np.load(dataPath + '/similarityMatrix_' + svType + '.npy', encoding='latin1', allow_pickle=True)
	#bagLabels = np.load(dataPath + '/bagLabels_' + svType + '.npy', encoding='latin1', allow_pickle=True)

	plot = False
	#don't make the plots for each feature to eliminate
	if featureElimination == "False" and leaveBagsOut == 'True':
		plot = True
	
		plotOutputPath = outDir + '/multipleInstanceLearning/rocCurves/'
		if not os.path.exists(plotOutputPath):
			os.makedirs(plotOutputPath)

		plotOutputFile = plotOutputPath + '/rocCurve_' + svType + '_leaveBagsOut.svg'
		outputFile = outDir + '/multipleInstanceLearning/performance_' + svType + '.txt'
		dataPath = outDir + '/multipleInstanceLearning/similarityMatrices/leaveBagsOut/'
		
		bagsCVClassification(dataPath, classifier, svType, title, plot, plotOutputFile, outputFile, False)
		
		#repeat, but then with random labels.
		#mind here, if multiple iterations, the baglabels are permanently shuffled!

		#ths is still broken
		#plotOutputFile = plotOutputPath + '/rocCurve_' + svType + '_randomLabels.svg'
		#outputFile = outDir + '/multipleInstanceLearning/performance_' + svType + '_randomBagLabels.txt'
		#cvClassification(dataPath, classifier, svType, title, plot, plotOutputFile, outputFile, True)
	elif featureElimination == 'True' and leaveOnePatientOut == 'False':
		
		#do the feature elimination here. Go through each similarity matrix with eliminated feature,
		#and obtain the classification scores.
		#the plot should not be made, and the AUCs should be written to a file
		
		#get each similarity matrix
		
		featureEliminationDataFolder = outDir + '/multipleInstanceLearning/similarityMatrices/featureSelection/'
		featureEliminationFiles = glob.glob(featureEliminationDataFolder + '*_' + svType + '_*')

		#

		#file to write to aucs to	
		outputFile = outDir + '/multipleInstanceLearning/featureEliminationResults_' + svType + '.txt'
		for fileInd in range(0, len(featureEliminationFiles)):



			#get the right similarity matrix for this shuffled feature
			#similarityMatrix = np.load(featureEliminationDataFolder + '/similarityMatrix_' + svType + '_' + str(fileInd) + '.npy')
			cvClassification(similarityMatrix, bagLabels, classifier, svType, title, plot, '', outputFile)
	elif featureElimination == 'False' and leaveOnePatientOut == 'True':
		
		leaveOneOutDataFolder = outDir + '/multipleInstanceLearning/similarityMatrices/leaveOnePatientOut/'
		plotOutputFile = outDir + '/multipleInstanceLearning/rocCurves/rocCurve_' + svType + '_leaveOnePatientOut.svg'
		leaveOnePatientOutCV(leaveOneOutDataFolder, classifier, svType, plotOutputFile, title)

	elif featureElimination == 'False' and leaveOneChromosomeOut == 'True':
		
		leaveOneChromosomeOutDataFolder = outDir + '/multipleInstanceLearning/similarityMatrices/leaveOneChromosomeOut/'
		plotOutputFile = outDir + '/multipleInstanceLearning/rocCurves/rocCurve_' + svType + '_leaveOneChromosomeOut.svg'
		leaveOneChromosomeOutCV(leaveOneChromosomeOutDataFolder, classifier, svType, plotOutputFile, title)
	
	else:
		
		print('Combination of options not implemented')
		exit(1)
		
		
