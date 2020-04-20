"""

	For the predicted labels in the lopoCV, check which predictions were correct
	How many of these are COSMIC?


"""

import sys
import glob
import re
import numpy as np
import random
from scipy import stats
path = sys.argv[2]
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVsGenes/')

import settings
from inputParser import InputParser

outDir = sys.argv[1]

svTypes = ['DEL', 'DUP', 'INV', 'ITX']

#get the cosmic genes
cosmicGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
cosmicGeneNames = []
for gene in cosmicGenes:
	cosmicGeneNames.append(gene[3].name)

nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], cosmicGenes) #In the same format as the causal genes.

#Combine the genes into one set.
allGenes = np.concatenate((cosmicGenes, nonCausalGenes), axis=0)

bcGeneNames = []
bcGenesFile = '../../data/genes/breastCancerCausalGenes.txt' #make setting
with open(bcGenesFile, 'r') as inF:

	for line in inF:
		line = line.strip()
		bcGeneNames.append(line)

#get the original labels for each patient
leaveOneOutDataFolder = outDir + '/multipleInstanceLearning/similarityMatrices/leaveOnePatientOut/'

def getCancerGeneEnrichment(outDir, leaveOneOutDataFolder, allGenes, cosmicGeneNames, svTypes):
	patientsWithCorrectCosmic = dict()
	patientsWithCorrectCosmicRandom = dict() #store per shuffled iteration.
	cosmicPairs = []
	totalCosmicGenes = 0
	for svType in svTypes:


		#get the predictions
		predOutFile = outDir + '/multipleInstanceLearning/leaveOnePatientOutCV/leaveOnePatientOutCV_' + svType + '.txt'

		perPatientPredictions = dict()
		with open(predOutFile, 'r') as inF:

			for line in inF:
				line = line.strip()
				splitLine = line.split('\t')

				patient = splitLine[0]
				predictions = splitLine[1:]
				perPatientPredictions[patient] = [float(i) for i in predictions]

		#get the original labels
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
		totalTP = 0
		totalFP = 0
		totalTN = 0
		totalFN = 0

		predictions = dict()
		for patient in patientFiles:



			for dataFile in patientFiles[patient]:

				if re.search('bagLabelsTest', dataFile):
					bagLabelsTest = np.load(dataFile, encoding='latin1', allow_pickle=True)
				if re.search('bagPairLabels', dataFile):
					bagPairLabels = np.load(dataFile, encoding='latin1', allow_pickle=True)

			for labelInd in range(0, len(bagPairLabels)):
				pairLabel = bagPairLabels[labelInd]
				splitLabel = pairLabel.split('_')

				if bagLabelsTest[labelInd] == 1:
					if splitLabel[0] in cosmicGeneNames:
						totalCosmicGenes += 1

				if bagLabelsTest[labelInd] == 1 and perPatientPredictions[patient][labelInd] == 1:
					pairLabel = bagPairLabels[labelInd]
					splitLabel = pairLabel.split('_')

					if splitLabel[0] in cosmicGeneNames:
						print(splitLabel[0])
						if splitLabel[7] not in patientsWithCorrectCosmic:
							patientsWithCorrectCosmic[splitLabel[7]] = 0
						patientsWithCorrectCosmic[splitLabel[7]] += 1
						cosmicPairs.append(splitLabel[0] + '_' + splitLabel[7] + '_' + svType)

				if bagLabelsTest[labelInd] == 1 and perPatientPredictions[patient][labelInd] == 1:
					totalTP += 1
				elif bagLabelsTest[labelInd] == 0 and perPatientPredictions[patient][labelInd] == 1:
					totalFP += 1
				elif bagLabelsTest[labelInd] == 1 and perPatientPredictions[patient][labelInd] == 0:
					totalFN += 1
				else:
					totalTN += 1

			#randomly sample as many genes as there are bags. Use these to check significance.
			for randIteration in range(0, 100):

				if randIteration not in patientsWithCorrectCosmicRandom:
					patientsWithCorrectCosmicRandom[randIteration] = dict()

				#sample as many random genes as there are bags
				randomInd = random.sample(range(0, allGenes.shape[0]), len(bagLabelsTest))
				randomGenes = allGenes[randomInd]

				#now do the same thing. For the indices that are true, how many genes are COSMIC?
				for labelInd in range(0, len(bagPairLabels)):

					if bagLabelsTest[labelInd] == 1 and perPatientPredictions[patient][labelInd] == 1:
						pairLabel = bagPairLabels[labelInd]
						splitLabel = pairLabel.split('_')
						geneName = randomGenes[labelInd, 3].name

						if geneName in cosmicGeneNames:
							if splitLabel[7] not in patientsWithCorrectCosmicRandom[randIteration]:
								patientsWithCorrectCosmicRandom[randIteration][splitLabel[7]] = 0
							patientsWithCorrectCosmicRandom[randIteration][splitLabel[7]] += 1


		#do this check here, because it is not implemented in the lopoCV properly and takes time to re-run.
		tpr = totalTP / (totalTP + totalFN)
		fpr = totalFP / (totalTN + totalFP)
		print(svType)
		print('tpr', tpr)
		print('fpr', fpr)

	print(totalCosmicGenes) #how many cosmic genes were linked to positive SV-gene pairs in total?


	patientCount = len(patientsWithCorrectCosmic)
	patientCountNegative = []
	for iteration in patientsWithCorrectCosmicRandom:
		negativeCount = len(patientsWithCorrectCosmicRandom[iteration])
		patientCountNegative.append(negativeCount)

	print(patientCount)
	print(patientCountNegative)


	patientCount = np.sum(list(patientsWithCorrectCosmic.values()))
	patientCountNegative = []
	for iteration in patientsWithCorrectCosmicRandom:
		negativeCount = np.sum(list(patientsWithCorrectCosmicRandom[iteration].values()))
		patientCountNegative.append(negativeCount)

	print(patientCount)
	print(patientCountNegative)

	z = patientCount - np.mean(patientCountNegative) / np.std(patientCountNegative)
	pValue = stats.norm.sf(abs(z))*2

	print(np.mean(patientCountNegative))
	print(np.std(patientCountNegative))
	print(z)
	print(pValue)

	return cosmicPairs

#first check for all cosmic genes
cosmicPairs = getCancerGeneEnrichment(outDir, leaveOneOutDataFolder, allGenes, cosmicGeneNames, svTypes)
#then check specific for breast cancer genes
bcPairs = getCancerGeneEnrichment(outDir, leaveOneOutDataFolder, allGenes, bcGeneNames, svTypes)

#check mutations. Are these genes frequently mutated in other patients?
mutDir = outDir + '/patientGeneMutationPairs/'
snvPatients = np.load(mutDir + 'snvPatients.npy', allow_pickle=True, encoding='latin1').item()

svPatientsDel = np.load(mutDir + 'svPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
svPatientsDup = np.load(mutDir + 'svPatientsDup.npy', allow_pickle=True, encoding='latin1').item()
svPatientsInv = np.load(mutDir + 'svPatientsInv.npy', allow_pickle=True, encoding='latin1').item()
svPatientsItx = np.load(mutDir + 'svPatientsItx.npy', allow_pickle=True, encoding='latin1').item()

cnvPatientsDel = np.load(mutDir + 'cnvPatientsDel.npy', allow_pickle=True, encoding='latin1').item()
cnvPatientsAmp = np.load(mutDir + 'cnvPatientsAmp.npy', allow_pickle=True, encoding='latin1').item()


noMutPairs = []
for pair in cosmicPairs:
	
	print(pair)
	
	splitPair = pair.split('_')
	svType = splitPair[2]
	
	if splitPair[0] in snvPatients[splitPair[1]]:
		print('snv')
	if splitPair[0] in cnvPatientsAmp[splitPair[1]]:
		print('cnv amp')
	if splitPair[0] in cnvPatientsDel[splitPair[1]]:
		print('cnv del')
	if splitPair[0] in svPatientsDel[splitPair[1]]:
		print('sv del')
	if splitPair[0] in svPatientsDup[splitPair[1]]:
		print('sv dup')
	if splitPair[0] in svPatientsInv[splitPair[1]]:
		print('sv inv')
	if splitPair[0] in svPatientsItx[splitPair[1]]:
		print('sv itx')
		

