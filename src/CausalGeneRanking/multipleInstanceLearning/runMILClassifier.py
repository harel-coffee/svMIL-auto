"""

	Here the MIL classifier performance will be tested on the input similarity matrices
	There are 3 CV possibilities:
	- leave-one-patient-out CV (also with random labels)
	- leave-one-chromosome-out CV (also possible with feature elimination)
	- leave-bags-out-CV

"""

import sys
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn import model_selection
from sklearn.metrics import plot_roc_curve, auc
import matplotlib.pyplot as plt
from scipy import interp
import os
import os.path
import glob
import re
from random import shuffle

import matplotlib
matplotlib.use('Agg')

featureElimination = sys.argv[2]
leaveOnePatientOut = sys.argv[3] #1 patient at a time in the test set
leaveOneChromosomeOut = sys.argv[4] #1 chromosome at a time in the test set
leaveBagsOut = sys.argv[5] #random bags in each CV fold
randomLabels = sys.argv[6] #running CV with randomized labels, only implemented for lopoCV

svTypes = ['DEL', 'DUP', 'INV', 'ITX']

outDir = sys.argv[1]

finalOutDir = outDir + '/rocCurves_fig3_S2/'

if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

def leaveOneChromosomeOutCV(leaveChromosomeOutDataFolder, classifier, svType, plotOutputFile, title, featureInd):
	"""
		Test classifier using leave-one-chromosome-out CV.

		Get the test sim matrix for each test chromosome from a file, and also get the training sim matrix for all other chromosomes.
		Also read the train/test labels.

		leaveChromosomeOutDataFolder (str): folder where the leave-one-chromosome-out data can be found
		classifier (object): classifier that we use to test performance
		svType (str): which SV type are we testing for?
		plotOutputFile (str): file where we output the roc curve to
		title (str): title of the plot
		featureInd (str): index of the feature to test for, specific for when we do feature elimination
	"""

	#first get the names of all patients
	if featureInd == 'False': #set to false to use regular chrCV
		allFiles = glob.glob(leaveChromosomeOutDataFolder + '*_' + svType + '.npy')
	else:
		allFiles = glob.glob(leaveChromosomeOutDataFolder + '*_' + svType + '_*_' + str(featureInd) + '.npy')

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
	totalTP = 0
	totalFP = 0
	totalTN = 0
	totalFN = 0
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
		predictions = classifier.predict(similarityMatrixTest)
		for labelInd in range(0, len(bagLabelsTest)):

			if bagLabelsTest[labelInd] == 1 and predictions[labelInd] == 1:
				totalTP += 1
			elif bagLabelsTest[labelInd] == 0 and predictions[labelInd] == 1:
				totalFP += 1
			elif bagLabelsTest[labelInd] == 1 and predictions[labelInd] == 0:
				totalFN += 1
			else:
				totalTN += 1

	
		viz = plot_roc_curve(classifier, similarityMatrixTest, bagLabelsTest,
							 name='roc',
							 alpha=0.3, lw=1, ax=ax)
		interp_tpr = interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(np.mean(viz.roc_auc))
		print('auc: ', np.mean(viz.roc_auc))

	print(np.mean(aucs))
	tpr = totalTP / (totalTP + totalFN)
	fpr = totalFP / (totalTN + totalFP)
	print('tpr', tpr)
	print('fpr', fpr)

	#write to the output file in case of feature elimination
	if featureInd != 'False':
		finalOutDir = outDir + '/multipleInstanceLearning/featureElimination/'
		if not os.path.exists(finalOutDir):
			os.makedirs(finalOutDir)

		outFile = finalOutDir + '/featureElimination_' + svType + '.txt'

		with open(outFile, 'a') as outF:

			outF.write(str(featureInd) + '\t' + str(np.mean(aucs)) + '\n')

	else:

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
			   title="Leave-one-chromosome-out CV: " + title)
		ax.legend(loc="lower right")
		plt.tight_layout()
		plt.savefig(plotOutputFile)


def bagsCVClassification(dataPath, clf, svType, title, plot, plotOutputFile, outputFile, shuffleLabels):
	"""
		Test classifier using leave-bags-out CV.

		Get the test sim matrix for each test bag from a file, and also get the training sim matrix for all other bags.
		Also read the train/test labels.

		dataPath (str): folder where the bag data can be found
		clf (object): classifier that we use to test performance
		svType (str): which SV type are we testing for?
		plot (boolean): if true we output a plot, otherwise not. 
		plotOutputFile (str): file where we output the roc curve to
		outputFile (str): obsolete
		title (str): title of the plot
		shuffleLabels (boolean): shuffle the labels when computing performance? 
	"""


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

		
		ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
			   title="Leave-random-bags-out CV: " + title)
		ax.legend(loc="lower right")
		plt.tight_layout()
		plt.savefig(plotOutputFile)


def leaveOnePatientOutCV(leaveOneOutDataFolder, classifier, svType, plotOutputFile, title, shuffleLabels):
	"""
		function to classify on a case where 1 patient is left out at a time the similarity matrices are pre-made, with each time 1 patient being left out.
		so we can just load these in, and then train our standard classifier on that. report on the performance of each patient.

		leaveOneOutDataFolder (str): folder where the per-patient data can be found
		classifier (object): classifier that we use to test performance
		svType (str): which SV type are we testing for?
		plotOutputFile (str): file where we output the roc curve to
		title (str): title of the plot
		shuffleLabels (boolean): shuffle the labels when computing performance?
	"""

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

		if shuffleLabels == True:
			shuffle(bagLabelsTrain)
			shuffle(bagLabelsTest)

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

	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
		   title="Leave-one-patient-out CV: " + title)
	ax.legend(loc="lower right")
	plt.tight_layout()
	plt.savefig(plotOutputFile)
	

	#output these values to a file
	if shuffleLabels == False:
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

	#Classifiers are the result of classifier optimization. Use the right one per SV type.
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


	plot = False #in feature elimination don't plot curves. 
	if featureElimination == "False" and leaveBagsOut == 'True': ### leave-bags-out CV
		plot = True


		plotOutputFile = finalOutDir + '/rocCurve_' + svType + '_leaveBagsOut.svg'
		outputFile = outDir + '/multipleInstanceLearning/performance_' + svType + '.txt'
		dataPath = outDir + '/multipleInstanceLearning/similarityMatrices/leaveBagsOut/'

		bagsCVClassification(dataPath, classifier, svType, title, plot, plotOutputFile, outputFile, False)
		
	elif featureElimination == 'True' and leaveOnePatientOut == 'False': #### Feature elimination scenario
		
		#do the feature elimination here. Go through each similarity matrix with eliminated feature,
		#and obtain the classification scores.
		
		featureEliminationDataFolder = outDir + '/multipleInstanceLearning/similarityMatrices/featureSelection/'
		#check the folder to get the max feature to check
		featureEliminationFiles = glob.glob(featureEliminationDataFolder + '/bagLabelsTest_' + svType + '*')

		#get the feature indices
		featureIndices = []
		for feFile in featureEliminationFiles:
			splitName = feFile.split('_')
			ind = splitName[len(splitName)-1]
			splitInd = ind.split('.')
			featureInd = splitInd[0]
			featureIndices.append(int(featureInd))

		#file to write to aucs to	
		outputFile = outDir + '/multipleInstanceLearning/featureEliminationResults_' + svType + '.txt'
		for featureInd in range(0, max(featureIndices)):
			#each time provide the feature index
			leaveOneChromosomeOutCV(featureEliminationDataFolder, classifier, svType, '', title, featureInd)
	elif featureElimination == 'False' and leaveOnePatientOut == 'True': ### Leave-one-patient-out CV
		
		shuffleLabels = False
		if randomLabels == 'True':
			shuffleLabels = True
			plotOutputFile = finalOutDir + '/rocCurve_' + svType + '_leaveOnePatientOut_random.svg'
		else:
			plotOutputFile = finalOutDir + '/rocCurve_' + svType + '_leaveOnePatientOut.svg'

		leaveOneOutDataFolder = outDir + '/multipleInstanceLearning/similarityMatrices/leaveOnePatientOut/'
		leaveOnePatientOutCV(leaveOneOutDataFolder, classifier, svType, plotOutputFile, title, shuffleLabels)

	elif featureElimination == 'False' and leaveOneChromosomeOut == 'True': ### leave-one-chromosome-out CV

		plotOutputFile = finalOutDir + '/rocCurve_' + svType + '_leaveOneChromosomeOut.svg'

		leaveOneChromosomeOutDataFolder = outDir + '/multipleInstanceLearning/similarityMatrices/leaveOneChromosomeOut/'
		leaveOneChromosomeOutCV(leaveOneChromosomeOutDataFolder, classifier, svType, plotOutputFile, title, 'False')
	
	else:
		
		print('Combination of options not implemented')
		exit(1)
		
		
