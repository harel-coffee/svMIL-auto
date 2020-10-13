"""
	Make figure 1 showing the most frequently affected cosmic genes in different cancer types

"""

import sys
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os.path
import pandas as pd
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
path = sys.argv[1]
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVsGenes/')

import settings
from inputParser import InputParser

class DriverPlotter:

	#the cancer types to show the results for
	cancerTypes = ['BRCA', 'PCAWG_OV', 'LUAD', 'LIHC', 'COAD']
	#use names to search for the tissue in the cosmic cancer type list
	cancerTypeNames = {'BRCA': 'breast', 'PCAWG_OV': 'ovarian', 'LUAD': 'lung',
					   'LIHC': 'liver', 'COAD': 'colorectal'}
	outDirPrefix = 'output/'
	svTypes = ['DEL', 'DUP', 'INV', 'ITX']

	def plotCosmicFrequencyScatter(self):

		correctPairsPerCancerType = dict()
		for cancerType in self.cancerTypes:
			cosmicGeneNames, cosmicGeneCancerTypes = self.getCosmicGenes()
			correctCosmicPairs = self.getCorrectlyPredictedCosmicPairs(cancerType, cosmicGeneNames)
			correctPairsPerCancerType[cancerType] = correctCosmicPairs

		#make the plot across all cancer types
		self.generateFrequencyScatterPlot(correctPairsPerCancerType, cosmicGeneCancerTypes)

	def plotAUC(self):

		#we should get this from a file automatically, this is just for the figure demo.
		auc = {'BRCA': [0.89, 0.89, 0.9, 0.81],
			   'PCAWG_OV': [0.76, 0.87, 0.87, 0.67],
			   'LUAD': [0.85, 0.63, 0.8, 0],
			   'LIHC': [0.65, 0.72, 0.76, 0.79],
			   'COAD': [0.72, 0.64, 0.8, 0.72]}

		#show the auc as dots in a scatterplot
		cancerTypeInd = 0
		plottedCancerTypes = []
		svTypeColors = ['#b5ffb9', '#f9bc86', '#a3acff', '#FF6B6C']
		#svTypeColors = [0, 1, 2, 3]
		jitter = [-0.03, -0.01, 0.01, 0.03]
		plotData = []
		for cancerType in auc:
			plottedCancerTypes.append(cancerType)

			for svTypeInd in range(0, len(self.svTypes)):
				plotData.append([cancerTypeInd+jitter[svTypeInd], auc[cancerType][svTypeInd], svTypeColors[svTypeInd]])
				#plt.scatter(cancerTypeInd + jitter[svTypeInd], auc[cancerType][svTypeInd], color=svTypeColors[svTypeInd])

			cancerTypeInd += 1

		#plotData = np.array(plotData)
		data = pd.DataFrame(plotData)
		data.columns = ['cancer type', 'AUC', 'color']
		data = data.drop_duplicates()
		print(data)
		#exit()
		sns.scatterplot(data=data, x='cancer type', y='AUC', hue=data.color,
						palette=sns.color_palette("Set1", data.color.nunique()), legend=False)

		plt.ylim([0.5,1])
		plt.xticks(np.arange(0, len(auc)), list(auc.keys()), rotation='vertical')
		plt.tight_layout()

		plt.show()

	def plotPathogenicSVFrequency(self):

		#read the line count of the pathogenic SV files.
		pathogenicSVCounts = dict()
		svTypeDistribution = dict()
		plotData = []
		for cancerType in self.cancerTypes:
			pathogenicSVCounts[cancerType] = 0
			svTypeDistribution[cancerType] = dict()
			svTypeDistribution[cancerType]['DEL'] = 0
			svTypeDistribution[cancerType]['DUP'] = 0
			svTypeDistribution[cancerType]['INV'] = 0
			svTypeDistribution[cancerType]['ITX'] = 0
			countsPerSample = dict()

			pathogenicSVFile = 'output/' + cancerType + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt'

			with open(pathogenicSVFile, 'r') as inF:
				for line in inF:
					pathogenicSVCounts[cancerType] += 1

					splitLine = line.split('\t')
					pair = splitLine[0]
					splitPair = pair.split('_')
					sample = splitPair[7]
					if sample not in countsPerSample:
						countsPerSample[sample] = 0
					countsPerSample[sample] += 1
					
					svType = splitPair[12]
					svTypeDistribution[cancerType][svType] += 1
					
			print(cancerType)
			print(len(countsPerSample))
			for sample in countsPerSample:
				plotData.append([cancerType, sample, countsPerSample[sample], pathogenicSVCounts[cancerType]])

			#plotData.append([cancerType, pathogenicSVCounts[cancerType]])
		exit()
		#plotData = np.array(plotData)
		data = pd.DataFrame(plotData)
		data.columns = ['cancerType', 'sample', 'sampleCount', 'svCount']
	
		#sns.scatterplot(data=data, x='cancerType', y='svCount', legend=False)
		v = sns.violinplot(data=data, x='cancerType', y='sampleCount', legend=False)
		
		
		# add n = X to show total count.
		
		cancerTypesWithCounts = []
		for cancerType in self.cancerTypes:
			cancerTypesWithCounts.append(cancerType + ' (N = ' + str(pathogenicSVCounts[cancerType]) + ')')

		plt.xticks(np.arange(0, len(self.cancerTypes)), cancerTypesWithCounts, rotation='vertical')
		plt.ylim([0, 750])
		plt.tight_layout()
		plt.show()
		
		#make the SV type bar chart
		
		#make a dictionary per SV type, where each array is then the cancer type.
		plotData = dict()
		for svType in ['DEL', 'DUP', 'INV', 'ITX']:
			plotData[svType] = []
			for cancerType in svTypeDistribution:
				plotData[svType].append(svTypeDistribution[cancerType][svType])
				
		df = pd.DataFrame(plotData)
		
		# From raw value to percentage
		totals = [i+j+k+l for i,j,k,l in zip(df['DEL'], df['DUP'], df['INV'], df['ITX'])]
		delBars = [i / j * 100 for i,j in zip(df['DEL'], totals)]
		dupBars = [i / j * 100 for i,j in zip(df['DUP'], totals)]
		invBars = [i / j * 100 for i,j in zip(df['INV'], totals)]
		itxBars = [i / j * 100 for i,j in zip(df['ITX'], totals)]
		 
		# plot
		barWidth = 0.85
		r = np.arange(0, len(self.cancerTypes))
		names = self.cancerTypes
		# Create green Bars
		plt.bar(r, delBars, color='#e41a1c', edgecolor='white', width=barWidth, label = 'DEL')
		# Create orange Bars
		plt.bar(r, dupBars, bottom=delBars, color='#377eb8', edgecolor='white', width=barWidth, label = 'DUP')
		# Create blue Bars
		plt.bar(r, invBars, bottom=[i+j for i,j in zip(delBars, dupBars)], color='#4daf4a', edgecolor='white', width=barWidth, label = 'INV')
		plt.bar(r, itxBars, bottom=[i+j+k for i,j,k in zip(delBars, dupBars, invBars)], color='#984ea3', edgecolor='white', width=barWidth, label = 'ITX')
		 
		# Custom x axis
		plt.xticks(r, names)
		#plt.xlabel("Cancer type")
		
		plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
		plt.tight_layout()
		# Show graphic
		plt.show()

			
			
		
		
		#plot here the total SV count
		# plotData = []
		# for cancerType in self.cancerTypes:
		# 	
		# 	if settings.general['source'] == 'PCAWG':
		# 		print("Reading SV data PCAWG")
		# 		svDir = settings.files['svDir']
		# 
		# 		svData = InputParser().getSVsFromFile_pcawg(svDir, cancerType)
		# 		plotData.append([cancerType, svData.shape[0]])
		# 
		# data = pd.DataFrame(plotData)
		# data.columns = ['cancerType', 'svCount']
		# 
		# sns.barplot(x="cancerType", y="svCount", data=data)
		# 
		# plt.xticks(np.arange(0, len(self.cancerTypes)), self.cancerTypes, rotation='vertical')
		# 
		# plt.show()

	def getCodingFrequency(self, gene, mutationPairs, codingFrequency):
		
		for patient in mutationPairs:
			if gene in mutationPairs[patient]:
				codingFrequency[patient] = 0
				
		return codingFrequency


	def generateFrequencyScatterPlot(self, allCosmicPairs, cosmicGeneCancerTypes):

		#Create an order for the genes and cancer types
		cancerTypesIndex = dict()
		cosmicGenesIndex = dict()
		geneFrequencies = dict()
		geneInd = 0
		for cancerTypeInd in range(0, len(allCosmicPairs)):
			cancerType = list(allCosmicPairs.keys())[cancerTypeInd]
			cancerTypesIndex[cancerType] = cancerTypeInd

			if cancerType not in geneFrequencies:
				geneFrequencies[cancerType] = dict()

			for pair in allCosmicPairs[cancerType]:
				splitPair = pair.split('_')
				gene = splitPair[0]
				if gene not in cosmicGenesIndex:
					cosmicGenesIndex[gene] = geneInd
					geneInd += 1
				if gene not in geneFrequencies[cancerType]:
					geneFrequencies[cancerType][gene] = 0
				geneFrequencies[cancerType][gene] += 1

		#instead of frequency by non-coding SVs, use number of coding events as size
		print('Calculating coding events...')
		codingFrequency = dict()
		normalizedCodingFrequency = dict()
		patientCounts = dict()
		
		#aside from normal codng events, also sample random genes to compare to
		iterationCount = 1
		#get all genes to sample from
		causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
		nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
		
		#Combine the genes into one set.
		allGenes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
		allGenes = nonCausalGenes
		randomCodingFrequency = dict()
		
		for cancerTypeInd in range(0, len(allCosmicPairs)):
			cancerType = self.cancerTypes[cancerTypeInd]
			
			snvPatients = np.load(self.outDirPrefix + '/' + cancerType + '/patientGeneMutationPairs/snvPatients.npy', encoding='latin1', allow_pickle=True).item()
			cnvPatientsAmp = np.load(self.outDirPrefix + '/' + cancerType + '/patientGeneMutationPairs/cnvPatientsAmp.npy', encoding='latin1', allow_pickle=True).item()
			cnvPatientsDel = np.load(self.outDirPrefix + '/' + cancerType + '/patientGeneMutationPairs/cnvPatientsDel.npy', encoding='latin1', allow_pickle=True).item()
			svPatientsDel = np.load(self.outDirPrefix + '/' + cancerType + '/patientGeneMutationPairs/svPatientsDel.npy', encoding='latin1', allow_pickle=True).item()
			svPatientsDup = np.load(self.outDirPrefix + '/' + cancerType + '/patientGeneMutationPairs/svPatientsDup.npy', encoding='latin1', allow_pickle=True).item()
			svPatientsInv = np.load(self.outDirPrefix + '/' + cancerType + '/patientGeneMutationPairs/svPatientsInv.npy', encoding='latin1', allow_pickle=True).item()
			svPatientsItx = np.load(self.outDirPrefix + '/' + cancerType + '/patientGeneMutationPairs/svPatientsItx.npy', encoding='latin1', allow_pickle=True).item()

			#check how many events are there.
			codingFrequency[cancerType] = dict()
			patientCounts[cancerType] = dict()
			normalizedCodingFrequency[cancerType] = dict()
			randomCodingFrequency[cancerType] = []
			for pair in allCosmicPairs[cancerType]:
				
				splitPair = pair.split('_')
				gene = splitPair[0]
				patient = splitPair[1]
				
				patientCounts[cancerType][patient] = 0
				
				if gene not in codingFrequency:
					codingFrequency[cancerType][gene] = 0
	
				codingPatients = dict()
				codingPatients = self.getCodingFrequency(gene, snvPatients, codingPatients)
				# codingPatients = self.getCodingFrequency(gene, cnvPatientsAmp, codingPatients)
				# codingPatients = self.getCodingFrequency(gene, cnvPatientsDel, codingPatients)
				# codingPatients = self.getCodingFrequency(gene, svPatientsDel, codingPatients)
				# codingPatients = self.getCodingFrequency(gene, svPatientsDup, codingPatients)
				# codingPatients = self.getCodingFrequency(gene, svPatientsInv, codingPatients)
				# codingPatients = self.getCodingFrequency(gene, svPatientsItx, codingPatients)
				
				codingFrequency[cancerType][gene] = len(codingPatients)
				
			#normalize the coding frequencies by the sample count with pathogenic SVs
			#for gene in codingFrequency[cancerType]:
			#	normalizedCodingFrequency[cancerType][gene] = len(codingFrequency[cancerType][gene]) / len(patientCounts[cancerType])
				#codingFrequency[cancerType][gene] = codingFrequency[cancerType][gene] / len(patientCounts[cancerType])
				
			for i in range(0, iterationCount):
				
				driverCount = len(codingFrequency[cancerType])
				
				#randomly sample driverCount genes
				randomGenes = np.random.choice(allGenes[:,3], driverCount)

				for gene in randomGenes:
					
					print(gene.name)
				
					
					geneCodingFrequency = dict()
					
					
					
					#get the coding events for these genes
					geneCodingFrequency = self.getCodingFrequency(gene.name, snvPatients, geneCodingFrequency)
					# geneCodingFrequency = self.getCodingFrequency(gene.name, cnvPatientsAmp, geneCodingFrequency)
					# geneCodingFrequency = self.getCodingFrequency(gene.name, cnvPatientsDel, geneCodingFrequency)
					# geneCodingFrequency = self.getCodingFrequency(gene.name, svPatientsDel, geneCodingFrequency)
					# geneCodingFrequency = self.getCodingFrequency(gene.name, svPatientsDup, geneCodingFrequency)
					# geneCodingFrequency = self.getCodingFrequency(gene.name, svPatientsInv, geneCodingFrequency)
					# geneCodingFrequency = self.getCodingFrequency(gene.name, svPatientsItx, geneCodingFrequency)
					
					if len(geneCodingFrequency) == 0:
						print('missing gene: ', gene.name)
						continue #skip for now because these mess up the statistics
					
					randomCodingFrequency[cancerType].append(len(geneCodingFrequency))
		
			break
		print(codingFrequency)
		print(randomCodingFrequency)
		print(np.mean(randomCodingFrequency['BRCA']))
		
		pValues = dict()
		for cancerType in codingFrequency:
			pValues[cancerType] = []
			
			uncorrectedPValues = []
			for gene in codingFrequency[cancerType]:
				print(gene)
				print(codingFrequency[cancerType][gene], np.mean(randomCodingFrequency[cancerType]))
				z = (codingFrequency[cancerType][gene] - np.mean(randomCodingFrequency[cancerType])) / np.std(randomCodingFrequency[cancerType])
				print(z)
				pValue = stats.norm.sf(abs(z))
				print(pValue)
				uncorrectedPValues.append([gene, z, pValue])
				
			#do mtc
			uncorrectedPValues = np.array(uncorrectedPValues, dtype = 'object')

			reject, pAdjusted, _, _ = multipletests(uncorrectedPValues[:,2], method='bonferroni') #fdr_bh or bonferroni
			
			print(reject)
			print(pAdjusted)
			exit()

			signPatients = []
			for pValueInd in range(0, len(uncorrectedPValues[:,2])):
				
				if reject[pValueInd] == True and uncorrectedPValues[pValueInd, 1] > 0:
			
					signPatients.append([uncorrectedPValues[pValueInd][0], uncorrectedPValues[pValueInd][1], pAdjusted[pValueInd]])
					
			signPatients = np.array(signPatients, dtype='object')
			pValues[cancerType] = signPatients
		
		print(pValues)
		exit()
		
		print(cancerTypesIndex)
		print(cosmicGenesIndex)
		
		
	

		#create the scatter plot in this order, use the frequency as point size
		genePlotIndices = dict()
		currentGenePlotIndex = 0
		plotData = []
		plotFrequencies = []
		pointColors = []
		for cancerType in allCosmicPairs:
			cancerTypeIndex = cancerTypesIndex[cancerType]
			cancerTypeName = self.cancerTypeNames[cancerType]
			for pair in allCosmicPairs[cancerType]:
				splitPair = pair.split('_')
				gene = splitPair[0]

				#get frequency of this gene
				#geneFrequency = geneFrequencies[cancerType][gene]
				#use frequency of coding events
				geneFrequency = normalizedCodingFrequency[cancerType][gene]
				if geneFrequency > 0:
					if gene not in genePlotIndices:
						genePlotIndices[gene] = currentGenePlotIndex
						currentGenePlotIndex += 1

					#determine the color based on if this gene is cancer-type specific
					edgecolors = 1
					facecolors = 'black'

					if re.search(cancerTypeName, cosmicGeneCancerTypes[gene], re.IGNORECASE):
						print('match', cancerType, gene)
						edgecolors = 2
						facecolors = 'red'

			
					plotData.append([genePlotIndices[gene], cancerTypeIndex, edgecolors, geneFrequency])
					plotFrequencies.append(geneFrequency)
					#plt.scatter(cancerTypeIndex, genePlotIndices[gene], s = geneFrequency*10)
					plt.scatter(genePlotIndices[gene], cancerTypeIndex, color=facecolors, s = geneFrequency*5)
		# plotData = np.array(plotData)
		# data = pd.DataFrame(plotData)
		# data.columns = ['gene', 'cancer type', 'color', 'frequency']
		# data = data.drop_duplicates()
		# print(data)
		# #exit()
		# sns.scatterplot(data=data, x='gene', y='cancer type', size='frequency', hue=data.color,
		# 				palette=sns.color_palette("Set2", data.color.nunique()), legend=False)

		#plt.yticks(np.arange(0, len(genePlotIndices)), list(genePlotIndices.keys()))
		#plt.xticks(np.arange(0, len(cancerTypesIndex)), list(cancerTypesIndex.keys()), rotation = 'vertical')

		plt.xticks(np.arange(0, len(genePlotIndices)), list(genePlotIndices.keys()), rotation = 'vertical')
		plt.yticks(np.arange(0, len(cancerTypesIndex)), list(cancerTypesIndex.keys()))

		plt.tight_layout()
		plt.savefig('frequency_scatter.svg')
		plt.show()


		return 0

	def getCosmicGenes(self):
		cosmicGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
		cosmicGeneNames = []
		cancerTypes = dict()
		for gene in cosmicGenes:
			cosmicGeneNames.append(gene[3].name)
			cancerTypes[gene[3].name] = gene[4]

		return cosmicGeneNames, cancerTypes


	#1. get the cosmic genes predicted correctly in the lopoCV
	def getCorrectlyPredictedCosmicPairs(self, cancerType, cosmicGeneNames):
		"""
			Iterate over the results of the lopoCV for each cancer type, and check if the
			predicted labels are correct. If the SV-gene pair is pathogenic AND the gene is
			also a known cosmic gene, report it.

			Parameters:
			- cancerType: cancer type to identify correctly predicted cosmic genes for

			Return:
			- correctlyPredictedCosmicGenes: correctly predicted cosmic genes
			- correctlyPredictedSpecificCosmicGenes: correctly predicted cosmic genes, specific for this cancer type

		"""
		
		leaveOneOutDataFolder = self.outDirPrefix + '/' + cancerType + '/multipleInstanceLearning/similarityMatrices/leaveOnePatientOut/'
		
		patientsWithCorrectCosmic = dict()
		cosmicPairs = []
		totalCosmicGenes = 0
		for svType in self.svTypes:

			#get the predictions
			predOutFile = self.outDirPrefix + '/' + cancerType + '/multipleInstanceLearning/leaveOnePatientOutCV/leaveOnePatientOutCV_' + svType + '.txt'

			#check if file exists, skip otherwise if the run failed for this sv type

			if os.path.isfile(predOutFile) is False:
				continue

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
							cosmicPairs.append(splitLabel[0] + '_' + splitLabel[7] + '_' + svType)





		return cosmicPairs
	
	
	
#2. Make the plot
#DriverPlotter().plotPathogenicSVFrequency()
#DriverPlotter().plotAUC()
DriverPlotter().plotCosmicFrequencyScatter()