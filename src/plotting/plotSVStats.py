"""
	Plot stats related to SVs.
	Fig S2A: Number of samples per cancer type
	Fig S2B: Number of SVs per cancer type
	Fig S2C: Number of pathogenic SVs per cancer type
	Fig S2D: Percentage of pathogenic SVs compared to total SVs per cancer type

	Fig 4B: Number of pathogenic SVs and % when running with CTCF loops

	Fig S4: Ratio of pathogenic SV types per cancer type

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os.path
import pandas as pd

#We need the settings file to determine the paths to where SVs are stored.
#so any file with svDir is ok to pass
path = sys.argv[1]
sys.path.insert(1, path)
sys.path.insert(1, 'linkSVsGenes/')
from inputParser import InputParser
import settings

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

def plotSVStatsPanels(cancerTypes, loopType):
	"""
		Plot panels A, B, C and D of Figure S2. Also C and D can be pt 1 and 2
		of figure 4B if CTCF data is provided.

		Parameters:
		- cancerTypes: cancerTypes to output for. Should equal output folder names
		- loopType: used to determine output figure name. TAD for Figure S3, CTCF for Figure 4B.
	"""

	#1. plot pathogenic SV count
	plotData = []
	cancerTypePlotNames = []
	for cancerType in cancerTypes:
		splitCancerType = cancerType.split('_')
		cancerType2 = '_'.join(splitCancerType[1:2])
		if loopType == 'TAD':
			cancerTypePlotNames.append(cancerType2)
		else:
			cancerTypePlotNames.append(cancerType2 + '_CTCF')

		#count how many pathogenic SVs we have
		pathogenicSVFile = 'output/' + cancerType + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt'

		pathogenicSVCount = 0
		with open(pathogenicSVFile, 'r') as inF:
			for line in inF:
				pathogenicSVCount += 1

		plotData.append([cancerType, pathogenicSVCount])

	pathogenicSVCounts = pd.DataFrame(plotData)
	pathogenicSVCounts.columns = ['Cancer type', 'Number of pathogenic SVs']

	#make bar plot
	ax = sns.barplot(data=pathogenicSVCounts, x="Cancer type", y="Number of pathogenic SVs", color='#a2d5f2')
	plt.xticks(np.arange(0, len(cancerTypes)), cancerTypePlotNames, rotation='vertical')

	plt.tight_layout()
	#Plot and save based on CTCF/TAD input
	if loopType == 'TAD':
		plt.savefig('output/figures/figureS2C.svg')
	else:
		plt.savefig('output/figures/figure4A_A.svg')
	plt.clf()

	#2. Make plot of the total SV counts
	plotData = []
	samplePlotData = []
	for cancerType in cancerTypes:

		#count the total number of SVs
		svDir = settings.files['svDir']
		svData = InputParser().getSVs_hmf(svDir, cancerTypeMetadataNames[cancerType])

		plotData.append([cancerType, svData.shape[0]])
		samplePlotData.append([cancerType, len(np.unique(svData[:,7]))])

	totalSVCounts = pd.DataFrame(plotData)
	totalSVCounts.columns = ['Cancer type', 'Number of SVs']

	#make bar plot
	ax = sns.barplot(x="Cancer type", y="Number of SVs", data=totalSVCounts, color='#07689f')
	plt.xticks(np.arange(0, len(cancerTypes)), cancerTypePlotNames, rotation='vertical')

	plt.tight_layout()
	if loopType == 'TAD':
		plt.savefig('output/figures/figureS2B.svg')
	plt.clf()


	#3. Make plot of the sample counts
	sampleCounts = pd.DataFrame(samplePlotData)
	sampleCounts.columns = ['Cancer type', 'Number of samples']

	#make bar plot
	ax = sns.barplot(x="Cancer type", y="Number of samples", data=sampleCounts, color='#ff7e67')
	plt.xticks(np.arange(0, len(cancerTypes)), cancerTypePlotNames, rotation='vertical')

	plt.tight_layout()
	if loopType == 'TAD':
		plt.savefig('output/figures/figureS2A.svg')
	plt.clf()

	#4. Show the relative % of pathogenic compared to total SVs.
	plotData = []
	for index, row in totalSVCounts.iterrows():

		cancerType = row['Cancer type']
		svGenePairFile = 'output/' + cancerType + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_'

		svGenePairCount = 0
		with open(svGenePairFile, 'r') as inF:
			for line in inF:
				svGenePairCount += 1

		pathogenicSVCount = pathogenicSVCounts[pathogenicSVCounts['Cancer type'] == row['Cancer type']]['Number of pathogenic SVs']
		relativeFrequency = (pathogenicSVCount / float(svGenePairCount)) * 100

		plotData.append([row['Cancer type'], relativeFrequency])

	sampleCounts = pd.DataFrame(plotData)
	sampleCounts.columns = ['Cancer type', 'Relative pathogenic SV frequency']

	#make bar plot
	ax = sns.barplot(x="Cancer type", y="Relative pathogenic SV frequency", data=sampleCounts, color='black')
	plt.xticks(np.arange(0, len(cancerTypes)), cancerTypePlotNames, rotation='vertical')

	plt.tight_layout()
	if loopType == 'TAD':
		plt.savefig('output/figures/figureS2D.svg')
	else:
		plt.savefig('output/figures/figure4A_B.svg')
		
	plt.clf()


def plotSVTypeDistribution(cancerTypes):
	"""
		Plot Figure S4, an overview of the ratio of each pathogenic SV type compared to
		all pathogenic SVs across cancer types.

		Parameters:
		- cancerTypes: cancer types to include in the plot. Equals output folder names.
	"""

	#make a dictionary per SV type, where each array is then the cancer type.

	#read the line count of the pathogenic SV files.
	pathogenicSVCounts = dict()
	svTypeDistribution = dict()
	plotData = []
	cancerTypePlotNames = []
	for cancerType in cancerTypes:
		#split names for clarity in plots
		splitCancerType = cancerType.split('_')
		cancerType2 = '_'.join(splitCancerType[1:2])
		cancerTypePlotNames.append(cancerType2)

		pathogenicSVCounts[cancerType] = 0
		svTypeDistribution[cancerType] = dict()
		svTypeDistribution[cancerType]['DEL'] = 0
		svTypeDistribution[cancerType]['DUP'] = 0
		svTypeDistribution[cancerType]['INV'] = 0
		svTypeDistribution[cancerType]['ITX'] = 0

		pathogenicSVFile = 'output/' + cancerType + '/linkedSVGenePairs/nonCoding_geneSVPairs.txt_pathogenicPairsFeatures.txt'

		with open(pathogenicSVFile, 'r') as inF:
			for line in inF:
				pathogenicSVCounts[cancerType] += 1

				splitLine = line.split('\t')
				pair = splitLine[0]
				splitPair = pair.split('_')

				svType = splitPair[12]
				svTypeDistribution[cancerType][svType] += 1

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
	r = np.arange(0, len(cancerTypes))
	names = cancerTypePlotNames
	# Create green Bars
	plt.bar(r, delBars, color='#e41a1c', edgecolor='white', width=barWidth, label = 'DEL')
	# Create orange Bars
	plt.bar(r, dupBars, bottom=delBars, color='#377eb8', edgecolor='white', width=barWidth, label = 'DUP')
	# Create blue Bars
	plt.bar(r, invBars, bottom=[i+j for i,j in zip(delBars, dupBars)], color='#ff6000ff', edgecolor='white', width=barWidth, label = 'INV')
	plt.bar(r, itxBars, bottom=[i+j+k for i,j,k in zip(delBars, dupBars, invBars)], color='#984ea3', edgecolor='white', width=barWidth, label = 'ITX')

	# Custom x axis
	plt.xticks(r, names, rotation='vertical')
	#plt.xlabel("Cancer type")

	plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
	plt.tight_layout()
	# Show graphic
	plt.savefig('output/figures/figureS4.svg')

cancerTypes = ['HMF_Breast_hmec', 'HMF_Ovary_ov', 'HMF_Lung_luad', 'HMF_Colorectal_coad',
			   'HMF_UrinaryTract_urinaryTract', 'HMF_Prostate_prostate', 'HMF_Esophagus_esophagus', 'HMF_Skin_skin',
			   'HMF_Pancreas_pancreas', 'HMF_Uterus_uterus', 'HMF_Kidney_kidney', 'HMF_NervousSystem_nervousSystem']

#1. Make the panels of Fig S2
#plotSVStatsPanels(cancerTypes, 'TAD')

#2. Make the panels of Fig 4B
cancerTypesCTCF = ['HMF_Breast_CTCF', 'HMF_Colorectal_CTCF', 'HMF_Lung_CTCF']
#plotSVStatsPanels(cancerTypesCTCF, 'CTCF')

#3. Make Fig S4
plotSVTypeDistribution(cancerTypes)