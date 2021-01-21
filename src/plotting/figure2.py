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

#We need the settings file to determine the paths to where SVs are stored.
#so any file with svDir is ok to pass
#path = sys.argv[1]
path = 'settings/settings_HMF_Breast_hmec/'
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVsGenes/')
from inputParser import InputParser
import settings

class Figure2:
	"""
		Class for generating all panels of figure 2.

		Panels:
		A: Predicted driver genes in each cancer type with significant driver potential
		B: Contribution of non-coding SV drivers vs SNV drivers across cancer types

	"""

	#because the cancer types in the metadata have slashes and spaces, we cannot use them, so use those
	#converted names here to read in the data.
	cancerTypeMetadataNames = {'HMF_Breast_hmec': 'Breast', 'HMF_Ovary_ov': 'Ovary', 'HMF_Lung_luad': 'Lung',
					   'HMF_Colorectal_coad': 'Colon/Rectum',
					   'HMF_UrinaryTract_urinaryTract': 'Urinary tract',
					   'HMF_Prostate_prostate': 'Prostate',
					   'HMF_Esophagus_esophagus': 'Esophagus',
					   'HMF_Skin_skin': 'Skin', 'HMF_Pancreas_pancreas': 'Pancreas',
					   'HMF_Uterus_uterus': 'Uterus', 'HMF_Kidney_kidney': 'Kidney',
					   'HMF_NervousSystem_nervousSystem': 'Nervous system',
					   'HMF_Breast_CTCF': 'Breast', 'HMF_Colorectal_CTCF': 'Colon/Rectum',
					   'HMF_Lung_CTCF': 'Lung'}

	#use names to search for the tissue in the cosmic cancer type list
	#this is using regex, so partial matches are ok if necessary
	cancerTypeNames = {'HMF_Breast_hmec': ['breast'], 'HMF_Ovary_ov': ['ovarian'],
					   'HMF_Lung_luad': ['lung'],'HMF_Colorectal_coad': ['colorectal'],
					   'HMF_UrinaryTract_urinaryTract': ['bladder'], 'HMF_Prostate_prostate': ['prostate'],
					   'HMF_Esophagus_esophagus': ['esophagus', 'esophageal'],
					   'HMF_Skin_skin': ['skin', 'melanoma'], 'HMF_Pancreas_pancreas': ['pancreas', 'pancreatic'],
					   'HMF_Uterus_uterus': ['uter'], 'HMF_Kidney_kidney': ['kidney', 'renal'],
					   'HMF_NervousSystem_nervousSystem': ['brain', 'nervous'],
					   'HMF_Breast_CTCF': ['breast'], 'HMF_Colorectal_CTCF': ['colorectal'],
					   'HMF_Lung_CTCF': ['lung']}

	def generatePathogenicSNVCounts(self, cancerTypes):
		"""
			For each provided cancer type, get the pathogenic SNV counts.

			Parameters:
			- cancerTypes: cancer types to get the high-impact SNVs for. In format: HMF_cancerType.*

			Return:
			- pathogenicSNVCounts: dictionary with the cancer type as key, and as value
			a dictionary with each gene as key, and as value the count of high-impact
			SNVs affecting that gene in this cancer type.
		"""
		pathogenicSNVCounts = dict()
		for cancerType in cancerTypes:

			#only the part after 'HMF' is in the metadata
			splitCancerType = cancerType.split('_')
			pathogenicSNVCounts[cancerType] = self.getPathogenicSNVsPerGene(splitCancerType[1])


		#Save the counts to use for later to save time
		np.save('pathogenicSNVCounts.npy', pathogenicSNVCounts)

		return pathogenicSNVCounts

	def getCodingFrequency(self, gene, mutationPairs, codingFrequency):
		"""
			Helper function to quickly re-format a dictionary with all patients as
			key. Used by @generateFrequencyScatterPlot
		"""

		for patient in mutationPairs:
			if gene in mutationPairs[patient]:
				codingFrequency[patient] = 0

		return codingFrequency


	def generateFrequencyScatterPlot(self, cancerTypes, pathogenicSNVCounts):
		"""
			Create figures S3A+B and figure 2A.

			Parameters:
			- cancerTypes: cancer types to include in the plot
			- pathogenicSNVCounts: dictionary with the cancer type as key, and as value
			a dictionary with each gene as key, and as value the count of high-impact
			SNVs affecting that gene in this cancer type.

		"""

		#Get the predicted positive SV-gene pairs (no longer cosmic-specific)
		allCosmicPairs = dict()
		for cancerType in cancerTypes:
			cosmicGeneNames, cosmicGeneCancerTypes = self.getCosmicGenes()
			correctCosmicPairs = self.getCorrectlyPredictedCosmicPairs(cancerType, cosmicGeneNames)
			allCosmicPairs[cancerType] = correctCosmicPairs


		#Create an order for the genes and cancer types
		cancerTypesIndex = dict()
		cosmicGenesIndex = dict()
		geneFrequencies = dict()
		geneInd = 0
		cancerTypePlotNames = []
		for cancerTypeInd in range(0, len(allCosmicPairs)):
			cancerType = list(allCosmicPairs.keys())[cancerTypeInd]
			cancerTypesIndex[cancerType] = cancerTypeInd

			splitCancerType = cancerType.split('_')
			cancerType2 = '_'.join(splitCancerType[1:2])
			cancerTypePlotNames.append(cancerType2)

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

		plotData = []
		plotDataAllGenes = []
		for cancerTypeInd in range(0, len(allCosmicPairs)):
			cancerType = list(allCosmicPairs.keys())[cancerTypeInd]
			cancerTypeNames = self.cancerTypeNames[cancerType]

			uniqueGenesC = dict()
			uniqueCosmicGenesC = dict()
			uniqueSpecificGenesC = dict()

			uniquePatients = dict()
			genesPerPatient = dict()

			for pair in allCosmicPairs[cancerType]:
				splitPair = pair.split('_')
				gene = splitPair[0]
				uniquePatients[splitPair[1]] = 0



				uniqueGenes[gene] = 0
				uniqueGenesC[gene] = 0
				geneType = 'Predicted driver gene'
				if gene in cosmicGeneCancerTypes:
					geneType = 'CGC gene'
					uniqueCosmicGenes[gene] = 0
					uniqueCosmicGenesC[gene] = 0
					for keyword in cancerTypeNames:
						if re.search(keyword, cosmicGeneCancerTypes[gene], re.IGNORECASE):
							geneType = 'Cancer-type specific CGC gene'
							uniqueSpecificGenes[gene] = 0
							uniqueSpecificGenesC[gene] = 0



				if splitPair[1] not in genesPerPatient:
					genesPerPatient[splitPair[1]] = []
				genesPerPatient[splitPair[1]].append(gene)

			print('cancer type: ', cancerType)
			print('genes: ', len(uniqueGenesC))
			print('cosmic genes: ', len(uniqueCosmicGenesC))
			print('specific genes: ', len(uniqueSpecificGenesC))
			print(uniqueSpecificGenesC)
			print('number of patients: ', len(uniquePatients))
			print('genes per patient: ', len(uniqueGenesC)/len(uniquePatients))

			perPatientGeneDistribution = []
			perPatientCosmicGeneDistribution = []
			perPatientSCosmicGeneDistribution = []
			for patient in genesPerPatient:

				geneCount = 0
				cosmicGeneCount = 0
				sCosmicGeneCount = 0
				for gene in genesPerPatient[patient]:
					geneCount += 1
					if gene in cosmicGeneCancerTypes:
						cosmicGeneCount += 1
						for keyword in cancerTypeNames:
							if re.search(keyword, cosmicGeneCancerTypes[gene], re.IGNORECASE):
								sCosmicGeneCount += 1

				perPatientGeneDistribution.append(geneCount)
				perPatientCosmicGeneDistribution.append(cosmicGeneCount)
				perPatientSCosmicGeneDistribution.append(sCosmicGeneCount)

				plotDataAllGenes.append([cancerType, 'Predicted driver genes', geneCount, patient])
				plotData.append([cancerType, 'CGC genes', cosmicGeneCount, patient])
				plotData.append([cancerType, 'Cancer type-specific CGC genes', sCosmicGeneCount, patient])
			
			
		print('total drivers: ', len(uniqueGenes))
		print('total known drivers: ', len(uniqueCosmicGenes))
		print('total specific drivers: ', len(uniqueSpecificGenes))

		#plot Fig S3A and S3B
		data = pd.DataFrame(plotData)
		data.columns = ['Cancer type', 'Gene type', 'Gene count per patient', 'Patient']

		v = sns.boxplot(y='Gene count per patient', x='Cancer type', data=data, hue='Gene type',
						palette=['#57db5f', '#5f57db'])

		plt.xticks(np.arange(0, len(cancerTypes)), cancerTypePlotNames, rotation='vertical')
		plt.tight_layout()

		plt.savefig('output/figures/figureS3A.svg')

		data = pd.DataFrame(plotDataAllGenes)
		data.columns = ['Cancer type', 'Gene type', 'Gene count per patient', 'Patient']

		v = sns.boxplot(y='Gene count per patient', x='Cancer type', data=data, hue='Gene type',
						palette=['#db5f57'])

		plt.xticks(np.arange(0, len(cancerTypes)), cancerTypePlotNames, rotation='vertical')
		plt.tight_layout()

		plt.savefig('output/figures/figureS3B.svg')
		
		
		####Then use the same information to output figure 2A
		
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
		
		#instead, sample 10.000 times X genes of the same set size
		#take the average of that set.

		np.random.seed(1)
		randomSampleIterations = 100

		geneFrequencies = dict()
		nonCodingOnlyGenes = dict()
		allPValues = []
		for cancerType in cancerTypes:

			#if checking results for CTCF, make sure that we can find the results
			#in the pathogenic SNV pairs data. 
			if cancerType == 'HMF_Breast_CTCF':
				cancerType2 = 'HMF_Breast'
			elif cancerType == 'HMF_Colorectal_CTCF':
				cancerType2 = 'HMF_Colorectal'
			elif cancerType == 'HMF_Lung_CTCF':
				cancerType2 = 'HMF_Lung'
			else:
				splitCancerType = cancerType.split('_')
				cancerType2 = '_'.join(splitCancerType[0:2])

			nonCodingOnlyGenes[cancerType] = dict()
			geneFrequencies[cancerType] = dict()

			trueGenes = dict()
			for pair in allCosmicPairs[cancerType]:

				splitPair = pair.split('_')
				gene = splitPair[0]
				trueGenes[gene] = 0

			randomDistribution = []
			for iteration in range(0, randomSampleIterations):

				#sample random genes of the same size.
				randomGenes = np.random.choice(allGeneNames, len(trueGenes))
				for gene in randomGenes:
					if gene in pathogenicSNVCounts[cancerType2]:
						randomDistribution.append(pathogenicSNVCounts[cancerType2][gene])
					else:
						randomDistribution.append(0)

			randomMean = np.mean(randomDistribution)
			randomStd = np.std(randomDistribution)
			pValues = []
			#allPValues = []
			for pair in allCosmicPairs[cancerType]:

				splitPair = pair.split('_')
				gene = splitPair[0]

				score = 0
				if gene in pathogenicSNVCounts[cancerType2]:
					score = pathogenicSNVCounts[cancerType2][gene]
				else:
					#don't count duplicates, that would be more than 1 per patient
					nonCodingOnlyGenes[cancerType][gene] = 0

				z = (score - randomMean) / randomStd

				pValue = stats.norm.sf(abs(z))
				pValues.append([gene, z, pValue])
				allPValues.append([gene, cancerType, z, pValue, score])



			if len(allPValues) < 1:
				continue
		uncorrectedPValues = np.array(allPValues, dtype='object')

		#sort by most significant first
		uncorrectedPValues = uncorrectedPValues[np.argsort(uncorrectedPValues[:,3])]

		reject, pAdjusted, _, _ = multipletests(uncorrectedPValues[:,3], method='fdr_bh', alpha=0.1) #fdr_bh or bonferroni

		signPatients = []
		for pValueInd in range(0, len(uncorrectedPValues[:,3])):

			gene = uncorrectedPValues[pValueInd, 0]
			cancerType = uncorrectedPValues[pValueInd, 1]

			if reject[pValueInd] == True and uncorrectedPValues[pValueInd, 2] > 0:
			
				geneFrequencies[cancerType][gene] = uncorrectedPValues[pValueInd, 2]

				signPatients.append([uncorrectedPValues[pValueInd][0], uncorrectedPValues[pValueInd][2], pAdjusted[pValueInd], uncorrectedPValues[pValueInd][3], uncorrectedPValues[pValueInd][4]])

		signPatients = np.array(signPatients, dtype='object')
		print(signPatients)


		#create the scatter plot in this order, use the frequency as point size
		genePlotIndices = dict()
		currentGenePlotIndex = 0
		plotData = []
		plotFrequencies = []
		pointColors = []
		cancerTypePlotNames = []
		for cancerType in allCosmicPairs:
			splitCancerType = cancerType.split('_')
			cancerType2 = '_'.join(splitCancerType[1:2])
			cancerTypePlotNames.append(cancerType2)
			cancerTypeIndex = cancerTypesIndex[cancerType]
			cancerTypeNames = self.cancerTypeNames[cancerType]
			for pair in allCosmicPairs[cancerType]:
				splitPair = pair.split('_')
				gene = splitPair[0]

				#get frequency of this gene
				if gene in geneFrequencies[cancerType] and gene in signPatients[:,0]:
					geneFrequency = geneFrequencies[cancerType][gene]
					
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
					
		plotData = np.array(plotData)
		data = pd.DataFrame(plotData)
		data.columns = ['Gene', 'Cancer type', 'color', 'frequency']
		data = data.drop_duplicates()

		#make sure to use the same colors as in the other plots, and not skip colors because
		#not all cacner types have significant genes.
		customPalette = sns.color_palette("hls", len(cancerTypes))
		finalPalette = []
		for colorInd in range(0, len(customPalette)):
			if colorInd not in plotData[:,1]:
				continue
			else:
				finalPalette.append(customPalette[colorInd])


		plt.figure(figsize=(10, 6))
		sns.scatterplot(data=data, x='Gene', y='Cancer type', size=data.frequency, hue=data['Cancer type'],
						legend=False, style=data.color, edgecolor = 'k', sizes=(20, 300),
						palette=finalPalette)
		
		plt.xticks(np.arange(0, len(genePlotIndices)), list(genePlotIndices.keys()), rotation = 'vertical')
		plt.yticks(np.arange(0, len(cancerTypesIndex)), cancerTypePlotNames)

		plt.tight_layout()
		plt.savefig('output/figures/figure2A.svg')
		plt.clf()


	def plotSVContributionPanel(self, cancerTypes, pathogenicSNVCounts):
		"""
			Figure 2B: For every cancer type, make a plot showing the count of predicted
			driver non-coding SVs versus the high-impact of SNVs in that cancer type.

			Parameters:
			- cancerTypes: list of cancer types to include in the plot
			- pathogenicSNVCounts: dictionary with cancer types as keys, and as
			values a dictionary with each gene and the count of high-impact SNVs
			in that gene in this cancer type.

		"""

		#Get the predicted non-coding SV drivers
		correctPairsPerCancerType = dict()
		for cancerType in cancerTypes:
			cosmicGeneNames, cosmicGeneCancerTypes = self.getCosmicGenes()
			correctCosmicPairs = self.getCorrectlyPredictedCosmicPairs(cancerType, cosmicGeneNames)
			correctPairsPerCancerType[cancerType] = correctCosmicPairs

		#check how many SNVs there are compared to ncSV drivers.
		plotData = []
		for cancerType in cancerTypes:
			splitCancerType = cancerType.split('_')
			cancerType2 = '_'.join(splitCancerType[0:2])
			geneSNVCounts = 0
			snvGenes = dict()
			for gene in pathogenicSNVCounts[cancerType2]:
				geneSNVCounts += pathogenicSNVCounts[cancerType2][gene]
				snvGenes[gene] = 0

			geneNcSVCounts = len(correctPairsPerCancerType[cancerType])
			ncSVGenes = dict()
			for pair in correctPairsPerCancerType[cancerType]:
				splitPair = pair.split('_')
				ncSVGenes[splitPair[0]] = 0

			#check which are in common, and which are unique.
			uniqueSNVs = np.setdiff1d(list(snvGenes.keys()), list(ncSVGenes.keys()))
			uniqueSVs = np.setdiff1d(list(ncSVGenes.keys()), list(snvGenes.keys()))
			intersect = np.intersect1d(list(snvGenes.keys()), list(ncSVGenes.keys()))

			plotData.append(['_'.join(splitCancerType[1:2]), len(uniqueSNVs), len(uniqueSVs)])

		#make the plot
		plt.figure(figsize=(6, 6))
		data = pd.DataFrame(plotData)
		data.columns = ['Cancer type', 'Genes affected by driver SNVs', 'Genes affected by driver non-coding SVs']
		filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
		sns.scatterplot(data=data, x='Genes affected by driver non-coding SVs',
						y='Genes affected by driver SNVs', hue=data['Cancer type'],
						palette=sns.color_palette("hls", len(cancerTypes)), legend='brief',
						s = 60, edgecolor = 'k', markers=filled_markers, style=data['Cancer type'])

		plt.tight_layout()
		plt.savefig('output/figures/figure2B.svg')

	def getCosmicGenes(self):
		"""
			Read all names of COSMIC genes into a dictionary. Keys are the names,
			associated cancer type as value.
		"""
		cosmicGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
		cosmicGeneNames = []
		cancerTypes = dict()
		for gene in cosmicGenes:
			cosmicGeneNames.append(gene[3].name)
			cancerTypes[gene[3].name] = gene[4]

		return cosmicGeneNames, cancerTypes


	#1. get the cosmic genes predicted correctly in the lopoCV
	def getCorrectlyPredictedCosmicPairs(self, cancerType, cosmicGeneNames, svTypes = ['DEL', 'DUP', 'INV', 'ITX']):
		"""
			Iterate over the results of the lopoCV for each cancer type, and check the
			predicted labels. If the SV-gene pair is pathogenic, report it.

			Used to look at correctly predicted (e.g. matching labels), but
			now only looks at all predicted positives. Also used to be
			specific for COSMIC genes, but now just reports all genes. 

			Parameters:
			- cancerType: cancer type to identify predicted genes for

			Return:
			- correctlyPredictedCosmicGenes: predicted genes
			- correctlyPredictedSpecificCosmicGenes: predicted genes, specific for this cancer type

		"""

		leaveOneOutDataFolder = 'output/' + cancerType + '/multipleInstanceLearning/similarityMatrices/leaveOnePatientOut/'

		patientsWithCorrectCosmic = dict()
		cosmicPairs = []
		totalCosmicGenes = 0
		for svType in svTypes:


			#get the predictions
			predOutFile = 'output/' + cancerType + '/multipleInstanceLearning/leaveOnePatientOutCV/leaveOnePatientOutCV_' + svType + '.txt'

			#check if file exists, skip otherwise if the run failed for this sv type

			if os.path.isfile(predOutFile) is False:
				continue


			perPatientPredictions = dict()
			with open(predOutFile, 'r') as inF:

				for line in inF:
					line = line.strip()
					splitLine = line.split('\t')

					pair = splitLine[0]

					splitPairLabel = pair.split("_")
					trueLabel = splitLine[1]
					prediction = splitLine[2]

					if prediction == "1":

						cosmicPairs.append(splitPairLabel[0] + '_' + splitPairLabel[7] + '_' + svType)

		return cosmicPairs

	#check if there are more high/moderate impact snvs in the patients than random genes
	def getPathogenicSNVsPerGene(self, cancerType):
		"""
			For the provided cancer type, go through the SNV files and return the
			genes that have SNVs with either HIGH or MODERATE impact.

			Parameters:
			- cancerType: list of cancer types to check for, should correspond to the names
			in the HMF metadata file.

			Return:
			- pathogenicSNVCounts: dictionary with gene names as keys and counts of
			high-impact SNVs as values.
		"""

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

#cancer types to create the plots for
cancerTypes = ['HMF_Breast_hmec', 'HMF_Ovary_ov', 'HMF_Lung_luad', 'HMF_Colorectal_coad',
			   'HMF_UrinaryTract_urinaryTract', 'HMF_Prostate_prostate', 'HMF_Esophagus_esophagus', 'HMF_Skin_skin',
			   'HMF_Pancreas_pancreas', 'HMF_Uterus_uterus', 'HMF_Kidney_kidney', 'HMF_NervousSystem_nervousSystem']


#1. Generate pathogenic SNV counts
#pathogenicSNVCounts = Figure2().generatePathogenicSNVCounts(cancerTypes)
## OR load them to save time
pathogenicSNVCounts = np.load('pathogenicSNVCounts.npy', encoding='latin1', allow_pickle=True).item()

#2. Generate figure 2A
#Figure2().generateFrequencyScatterPlot(cancerTypes, pathogenicSNVCounts)

#3. Generate figure 2B
Figure2().plotSVContributionPanel(cancerTypes, pathogenicSNVCounts)
