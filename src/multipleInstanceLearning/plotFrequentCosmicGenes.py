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
import gzip
path = sys.argv[1]
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVsGenes/')

import settings
from inputParser import InputParser

class DriverPlotter:

	#the cancer types to show the results for
	cancerTypes = ['HMF_Breast', 'HMF_Ovary', 'HMF_Liver', 'HMF_Lung', 'HMF_Colorectal',
				   'HMF_UrinaryTract', 'HMF_Prostate', 'HMF_Esophagus', 'HMF_Skin',
				   'HMF_Pancreas', 'HMF_Uterus', 'HMF_Kidney', 'HMF_NervousSystem']
	#because the cancer types in the metadata have slashes and spaces, we cannot use them, so use those
	#converted names here to read in the data.
	cancerTypeMetadataNames = {'HMF_Breast': 'Breast', 'HMF_Ovary': 'Ovary', 'HMF_Lung': 'Lung',
					   'HMF_Liver': 'Liver', 'HMF_Colorectal': 'Colon/Rectum',
					   'HMF_UrinaryTract': 'Urinary tract', 'HMF_Prostate': 'Prostate',
					   'HMF_Esophagus': 'Esophagus',
					   'HMF_Skin': 'Skin', 'HMF_Pancreas': 'Pancreas',
					   'HMF_Uterus': 'Uterus', 'HMF_Kidney': 'Kidney',
					   'HMF_NervousSystem': 'Nervous system'}

	#cancerTypes = ['HMF_Colorectal']

	#use names to search for the tissue in the cosmic cancer type list
	#this is using regex, so partial matches are ok if necessary
	cancerTypeNames = {'HMF_Breast': ['breast'], 'HMF_Ovary': ['ovarian'], 'HMF_Lung': ['lung'],
					   'HMF_Liver': ['liver'], 'HMF_Colorectal': ['colorectal'],
					   'HMF_UrinaryTract': ['bladder'], 'HMF_Prostate': ['prostate'],
					   'HMF_Esophagus': ['esophagus', 'esophageal'],
					   'HMF_Skin': ['skin', 'melanoma'], 'HMF_Pancreas': ['pancreas', 'pancreatic'],
					   'HMF_Uterus': ['uter'], 'HMF_Kidney': ['kidney', 'renal'],
					   'HMF_NervousSystem': ['brain', 'nervous']}
	outDirPrefix = 'output/'
	svTypes = ['DEL', 'DUP', 'INV', 'ITX']

	def plotCosmicFrequencyScatter(self):

		correctPairsPerCancerType = dict()
		pathogenicSNVCounts = dict()
		for cancerType in self.cancerTypes:
			cosmicGeneNames, cosmicGeneCancerTypes = self.getCosmicGenes()
			correctCosmicPairs = self.getCorrectlyPredictedCosmicPairs(cancerType, cosmicGeneNames)
			correctPairsPerCancerType[cancerType] = correctCosmicPairs

			#only the part after 'HMF' is in the metadata
			splitCancerType = cancerType.split('_')
			#pathogenicSNVCounts[cancerType] = DriverPlotter().getPathogenicSNVsPerGene(splitCancerType[1])

		#save the pathogenic snv counts for later
		#np.save('pathogenicSNVCounts.npy', pathogenicSNVCounts)

		pathogenicSNVCounts = np.load('pathogenicSNVCounts.npy', encoding='latin1', allow_pickle=True).item()

		#make the plot across all cancer types
		self.generateFrequencyScatterPlot(correctPairsPerCancerType, cosmicGeneCancerTypes, pathogenicSNVCounts)

	def plotAUC(self):

		#we should get this from a file automatically, this is just for the figure demo.
		# auc = {'BRCA': [0.89, 0.89, 0.9, 0.81],
		# 	   'PCAWG_OV': [0.76, 0.87, 0.87, 0.67],
		# 	   'LUAD': [0.85, 0.63, 0.8, 0],
		# 	   'LIHC': [0.65, 0.72, 0.76, 0.79],
		# 	   'COAD': [0.72, 0.64, 0.8, 0.72]}

		auc = {'HMF_Breast': [0.95, 0.87, 0.89, 0.89],
			   'HMF_Ovary': [0.82, 0.74, 0.73, 0.72],
			   'HMF_Liver': [0.89, 0, 0.43, 0.93],
			   'HMF_Lung': [0.85, 0.86, 0.79, 0.79],
			   'HMF_Colorectal': [0.92, 0.84, 0.73, 0.74],

			   'HMF_UrinaryTract': [0.95, 0.92, 0.92, 0.91],
			   'HMF_Prostate': [0.73, 0.86, 0.76, 0.66],
			   'HMF_Esophagus': [0.91, 0.65, 0.76, 0.69],
			   'HMF_Skin': [0.88, 0.86, 0.87, 0.80],
			   'HMF_Pancreas': [0.79, 0.84, 0.82, 0.82],
			   'HMF_Uterus': [0.5, 0.83, 0.56, 0.62],
			   'HMF_Kidney': [0.7, 0.8, 0.83, 0.96],
			   'HMF_NervousSystem': [0.76, 0.69, 0.83, 0.75],
			   }

		#show the auc as dots in a scatterplot
		cancerTypeInd = 0
		plottedCancerTypes = []
		svTypeColors = ['#b5ffb9', '#f9bc86', '#a3acff', '#FF6B6C']
		#svTypeColors = [0, 1, 2, 3]
		jitter = [-0.15, -0.05, 0.05, 0.15]
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
		fig, ax = plt.subplots(1,1)
		sns.scatterplot(data=data, x='cancer type', y='AUC', hue=data.color,
						palette=sns.color_palette("Set1", data.color.nunique()), legend=False,
						s = 60, edgecolor = 'k')

		#set separators

		ax.set_xticks(np.arange(0, len(auc)-1)+0.5, minor=True)
		ax.grid(b=True, which='minor', linewidth=0.5, linestyle='--')

		plt.ylim([0.4,1])
		plt.xticks(np.arange(0, len(auc)), list(auc.keys()), rotation='vertical')
		plt.tight_layout()

		plt.show()

	def plotPathogenicSVFrequency(self):

		# plotData = dict()
		# plotData['pathogenicSVs'] = []
		# plotData['totalSVs'] = []
		plotData = []
		for cancerType in self.cancerTypes:

			#count how many pathogenic SVs we have
			pathogenicSVFile = 'output/' + cancerType + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt'

			pathogenicSVCount = 0
			with open(pathogenicSVFile, 'r') as inF:
				for line in inF:
					pathogenicSVCount += 1

			plotData.append([cancerType, pathogenicSVCount])

			#plotData['pathogenicSVs'].append(pathogenicSVCount)

			#count the total number of SVs
			# svDir = settings.files['svDir']
			# svData = InputParser().getSVs_hmf(svDir, self.cancerTypeMetadataNames[cancerType])
			# #plotData['totalSVs'].append(svData.shape[0])
			#
			#
			# plotData.append([cancerType, svData.shape[0], 'SV'])

		data = pd.DataFrame(plotData)
		data.columns = ['cancerType', 'svCount']

		#make bar plot
		ax = sns.barplot(x="cancerType", y="svCount", data=data, color='#a2d5f2')
		plt.xticks(np.arange(0, len(self.cancerTypes)), self.cancerTypes, rotation='vertical')

		plt.tight_layout()
		# Show graphic
		plt.show()
		exit()

		plotData = []
		samplePlotData = []
		for cancerType in self.cancerTypes:

			#count the total number of SVs
			svDir = settings.files['svDir']
			svData = InputParser().getSVs_hmf(svDir, self.cancerTypeMetadataNames[cancerType])
			#plotData['totalSVs'].append(svData.shape[0])


			plotData.append([cancerType, svData.shape[0]])
			samplePlotData.append([cancerType, len(np.unique(svData[:,7]))])

		data = pd.DataFrame(plotData)
		data.columns = ['cancerType', 'svCount']

		#make bar plot
		ax = sns.barplot(x="cancerType", y="svCount", data=data, color='#07689f')
		plt.xticks(np.arange(0, len(self.cancerTypes)), self.cancerTypes, rotation='vertical')

		plt.tight_layout()
		# Show graphic
		plt.show()

		data = pd.DataFrame(samplePlotData)
		data.columns = ['cancerType', 'sampleCount']

		#make bar plot
		ax = sns.barplot(x="cancerType", y="sampleCount", data=data, color='#ff7e67')
		plt.xticks(np.arange(0, len(self.cancerTypes)), self.cancerTypes, rotation='vertical')

		plt.tight_layout()
		# Show graphic
		plt.show()
		exit()

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
		#exit()
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
		plt.ylim([0, 200])
		plt.tight_layout()
		plt.show()


		###make a plot showing how many pathogenic SVs vs total SVs



		




		
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
		plt.xticks(r, names, rotation='vertical')
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


	def generateFrequencyScatterPlot(self, allCosmicPairs, cosmicGeneCancerTypes, pathogenicSNVCounts):

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

		#check distribution of genes/cosmic etc
		uniqueGenes = dict()
		uniqueCosmicGenes = dict()
		uniqueSpecificGenes = dict()
		for cancerTypeInd in range(0, len(allCosmicPairs)):
			cancerType = list(allCosmicPairs.keys())[cancerTypeInd]
			cancerTypeNames = self.cancerTypeNames[cancerType]

			for pair in allCosmicPairs[cancerType]:
				splitPair = pair.split('_')
				gene = splitPair[0]

				uniqueGenes[gene] = 0
				if gene in cosmicGeneCancerTypes:
					uniqueCosmicGenes[gene] = 0
					for keyword in cancerTypeNames:
						if re.search(keyword, cosmicGeneCancerTypes[gene], re.IGNORECASE):
							uniqueSpecificGenes[gene] = 0

		print('total drivers: ', len(uniqueGenes))
		print('total known drivers: ', len(uniqueCosmicGenes))
		print('total specific drivers: ', len(uniqueSpecificGenes))
		exit()
				






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
		
		allGeneNames = []
		for gene in allGenes:
			allGeneNames.append(gene[3].name)
		cosmicGeneNames = []
		for gene in causalGenes:
			cosmicGeneNames.append(gene[3].name)
		#allGenes = nonCausalGenes

		np.random.seed(42)
		randomGenes = np.random.choice(allGeneNames, 100)
		
		geneFrequencies = dict()
		nonCodingOnlyGenes = dict()
		allPValues = []
		for cancerType in self.cancerTypes:
			nonCodingOnlyGenes[cancerType] = dict()
			geneFrequencies[cancerType] = dict()

			randomDistribution = []
			for gene in randomGenes:

				if gene in pathogenicSNVCounts[cancerType]:
					randomDistribution.append(pathogenicSNVCounts[cancerType][gene])
				else:
					randomDistribution.append(0)
			# print(cancerType)
			# print(randomDistribution)
			# print(np.mean(randomDistribution))

			pValues = []
			for pair in allCosmicPairs[cancerType]:

				splitPair = pair.split('_')
				gene = splitPair[0]

				score = 0
				if gene in pathogenicSNVCounts[cancerType]:
					#print(gene, ': ', pathogenicSNVCounts[cancerType][gene])
					score = pathogenicSNVCounts[cancerType][gene]
				else:
					#print(gene, ' not pathogenic')
					#don't count duplicates, that would be more than 1 per patient
					nonCodingOnlyGenes[cancerType][gene] = 0

				z = (score - np.mean(randomDistribution)) / np.std(randomDistribution)

				pValue = stats.norm.sf(abs(z))
				pValues.append([gene, z, pValue])
				allPValues.append([gene, cancerType, z, pValue])

			#uncorrectedPValues = np.array(pValues, dtype = 'object')

		#adjust across cancer types
		
		#print(uncorrectedPValues)
		
		uncorrectedPValues = np.array(allPValues, dtype='object')

		reject, pAdjusted, _, _ = multipletests(uncorrectedPValues[:,3], method='bonferroni') #fdr_bh or bonferroni

		signPatients = []
		for pValueInd in range(0, len(uncorrectedPValues[:,3])):
			
			gene = uncorrectedPValues[pValueInd, 0]
			cancerType = uncorrectedPValues[pValueInd, 1]



			if reject[pValueInd] == True and uncorrectedPValues[pValueInd, 2] > 0:

				geneFrequencies[cancerType][gene] = uncorrectedPValues[pValueInd, 2]

				signPatients.append([uncorrectedPValues[pValueInd][0], uncorrectedPValues[pValueInd][2], pAdjusted[pValueInd]])

		signPatients = np.array(signPatients, dtype='object')

		print(signPatients)

	

		#create the scatter plot in this order, use the frequency as point size
		genePlotIndices = dict()
		currentGenePlotIndex = 0
		plotData = []
		plotFrequencies = []
		pointColors = []
		for cancerType in allCosmicPairs:
			cancerTypeIndex = cancerTypesIndex[cancerType]
			cancerTypeNames = self.cancerTypeNames[cancerType]
			for pair in allCosmicPairs[cancerType]:
				splitPair = pair.split('_')
				gene = splitPair[0]

				#get frequency of this gene
				if gene in geneFrequencies[cancerType]:
					geneFrequency = geneFrequencies[cancerType][gene]
					#use frequency of coding events
					#geneFrequency = normalizedCodingFrequency[cancerType][gene]
					#3.5
				
					if gene not in genePlotIndices:
						genePlotIndices[gene] = currentGenePlotIndex
						currentGenePlotIndex += 1

					#determine the color based on if this gene is cancer-type specific
					edgecolors = 1
					facecolors = 'black'

					if gene in cosmicGeneCancerTypes:
						facecolors = 'green'
						edgecolors = 3
						for keyword in cancerTypeNames:
							if re.search(keyword, cosmicGeneCancerTypes[gene], re.IGNORECASE):
								print('match', cancerType, gene)
								edgecolors = 2
								facecolors = 'red'


					plotData.append([genePlotIndices[gene], cancerTypeIndex, edgecolors, geneFrequency*500])
					
					#plt.scatter(genePlotIndices[gene], cancerTypeIndex, color=facecolors, s = geneFrequency*5)
		plotData = np.array(plotData)
		data = pd.DataFrame(plotData)
		data.columns = ['gene', 'cancerType', 'color', 'frequency']
		data = data.drop_duplicates()
		print(data)
		#exit()
		sns.scatterplot(data=data, x='gene', y='cancerType', size=data.frequency, hue=data.cancerType,
						legend=False, style=data.color, edgecolor = 'k', sizes=(20, 300),
						palette=sns.color_palette("hls", data.cancerType.nunique()))


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

					#if bagLabelsTest[labelInd] == 1 and perPatientPredictions[patient][labelInd] == 1:
					if perPatientPredictions[patient][labelInd] == 1:
						pairLabel = bagPairLabels[labelInd]
						splitLabel = pairLabel.split('_')

						#if splitLabel[0] in cosmicGeneNames:
						
						cosmicPairs.append(splitLabel[0] + '_' + splitLabel[7] + '_' + svType)




		return cosmicPairs

	#check if there are more high/moderate impact snvs in the patients than random genes
	def getPathogenicSNVsPerGene(self, cancerType):
		
		#doesn't work to use that as the folder name... 
		if cancerType == 'Colorectal':
			cancerType = 'Colon/Rectum'
		if cancerType == 'UrinaryTract':
			cancerType = 'Urinary tract'
		if cancerType == 'NervousSystem':
			cancerType = 'Nervous system'

		pathogenicSNVCounts = dict()

		#get the samples of this cancer type
		metadataFile = settings.files['metadataHMF']
		snvDir = settings.files['snvDir']

		#save the IDs of the patients with this cancer type
		cancerTypeIds = dict()
		with open(metadataFile, 'rb') as inF:

			for line in inF:
				line = line.decode('ISO-8859-1')

				splitLine = line.split('\t')
				if splitLine[6] == cancerType:
					sampleId = splitLine[1]
					patientId = splitLine[0]

					cancerTypeIds[sampleId] = patientId

		#get the snv files
		for sampleId in cancerTypeIds:
			patientId = cancerTypeIds[sampleId]

			matchedFile = glob.glob(snvDir + '/*_' + patientId + '/' + sampleId + '.purple.somatic.vcf.gz')[0]

			#open the .gz file
			with gzip.open(matchedFile, 'rb') as inF:

				for line in inF:
					line = line.strip().decode('ISO-8859-1')

					if re.search('^#', line): #skip header
						continue

					#skip the SV if it did not pass.
					splitLine = line.split("\t")
					filterInfo = splitLine[6]
					if filterInfo != 'PASS':
						continue

					#check if there is a pathogenic snv
					if re.search('HIGH', line) or re.search('MODERATE', line):

						#Check if this SNV has any affiliation with a gene. This means that in the info field, a gene is mentioned somewhere. That is, there is an ENSG identifier.
						infoField = splitLine[7]

						geneSearch = re.search('(ENSG\d+)', infoField)
						if geneSearch:
							#the gene name is always directly before the ENSG identifier
							geneMatch = re.search('.+[\|\=\(](.+)?\|ENSG\d+', infoField).group(1)

							if geneMatch not in pathogenicSNVCounts:
								pathogenicSNVCounts[geneMatch] = 0
							pathogenicSNVCounts[geneMatch] += 1


		
		return pathogenicSNVCounts

	
	
	
#2. Make the plot
DriverPlotter().plotPathogenicSVFrequency()
#DriverPlotter().plotAUC()
#DriverPlotter().plotCosmicFrequencyScatter()