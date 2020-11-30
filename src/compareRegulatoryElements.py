# load regulatory elements of one cancer type, and compare it to others

import sys
import numpy as np
import glob
import re
import os
path = sys.argv[1]
sys.path.insert(1, path) #path to the settings file
import settings
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

#this code depends on the input parser from the linking part. This is quick and dirty, is there a better solution?
sys.path.insert(1, 'linkSVsGenes/')

from inputParser import InputParser
from neighborhoodDefiner import NeighborhoodDefiner

finalOutDir = 'output/featureComparisons/'
if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

def readBedFile(elementFile):

	elementData = []
	with open(elementFile, "r") as f:
		lineCount = 0
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")

			chrom = splitLine[0]
			if chrom == 'chr23':
				chrom = 'chrX'
			elif chrom == 'chr24':
				chrom = 'chrY'
			elif chrom == 'chr25':
				chrom = 'chrM'

			#chr, start, end
			elementData.append([chrom, int(splitLine[1]), int(splitLine[2])])

	elementData = np.array(elementData, dtype='object')

	return elementData

def getEnhancersFromFile(enhancerFile):

	enhancers = []
	with open(enhancerFile, 'r') as f:

		lineCount = 0
		for line in f:
			if lineCount < 1:
				lineCount += 1
				continue

			line = line.strip()
			splitLine = line.split("\t")

			interaction = splitLine[0]
			splitInteraction = interaction.split(",") #first part is the enhancer, 2nd part the gene

			enhancerInfo = splitInteraction[0]
			splitEnhancerInfo = enhancerInfo.split(":")
			#Add the chr notation for uniformity.
			chrMatch = re.search("chr", splitEnhancerInfo[0], re.IGNORECASE)
			chrName = ""
			if chrMatch is None:
				chrName = "chr" + splitEnhancerInfo[0]
			else:
				chrName = splitEnhancerInfo[0]
			splitPosInfo = splitEnhancerInfo[1].split("-") #get the positions
			start = int(splitPosInfo[0])
			end = int(splitPosInfo[1])


			element = [chrName, start, end]
			enhancers.append(element)

	return np.array(enhancers, dtype='object')

causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.

#Combine the genes into one set.
genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)

#all cancer types for easy loading
cancerTypes = dict()
cancerTypes['hmec'] = 0
cancerTypes['ov'] = 1
cancerTypes['coad'] = 2
cancerTypes['esophagus'] = 3
cancerTypes['kidney'] = 4
cancerTypes['luad'] = 5
cancerTypes['nervousSystem'] = 6
cancerTypes['pancreas'] = 7
cancerTypes['prostate'] = 8
cancerTypes['skin'] = 9
cancerTypes['urinaryTract'] = 10
cancerTypes['uterus'] = 11
cancerTypes['gm12878'] = 12

regulatoryElements = ['tads', 'enhancers', 'superEnhancers', 'h3k4me1', 'h3k27ac',
					  'h3k27me3', 'h3k4me3', 'ctcf', 'dnase', 'rnapol']

#have a mal that we fill up for each chromosome depending on the size
#get the sizes of each chromosome
#this can be done better without re-reading the file but no time
def getChromosomeMal():

	coordinates = dict()
	with open("../data/chromosomes/hg19Coordinates.txt", 'r') as inF:
		lineCount = 0
		for line in inF:
			line = line.strip()
			if lineCount < 1:
				lineCount += 1
				continue
			splitLine = line.split('\t')
			chrEnd = splitLine[3]
			chrId = splitLine[1]
			coordinates["chr" + chrId] = int(chrEnd)

	malDict = dict()
	for chrom in coordinates:
		#mal = np.zeros(coordinates[chrom])
		malDict[chrom] = dict()
	return malDict



#Option 1: make vectors and correlate, or compute absolute distance
for regulatoryElement in regulatoryElements:
	print('Comparing: ', regulatoryElement)

	differenceMatrix = np.zeros([len(cancerTypes), len(cancerTypes)])
	for cancerType1 in cancerTypes:
		print('cancer type: ', cancerType1)

		#load the regulatory elements
		superEnhancerFile = glob.glob("../data/" + regulatoryElement + "/" + cancerType1 + "/*")
		#select for the clustered one
		if len(superEnhancerFile) > 1:
			superEnhancerFile = glob.glob("../data/" + regulatoryElement + "/" + cancerType1 + "/*_clustered.bed")

		if len(superEnhancerFile) < 1:
			continue

		if regulatoryElement == 'enhancers':
			superEnhancerData = getEnhancersFromFile(superEnhancerFile[0])
		else:
			superEnhancerData = readBedFile(superEnhancerFile[0])

		malDictTissue1 = getChromosomeMal()
		for row in superEnhancerData:
			chrom = row[0]
			malDictTissue1[chrom][int(row[1])] = 1


		for cancerType2 in cancerTypes:
			if cancerType1 == cancerType2:
				continue

			superEnhancerFile2 = glob.glob("../data/" + regulatoryElement + "/" + cancerType2 + "/*")
			#select for the clustered one
			if len(superEnhancerFile2) > 1:
				superEnhancerFile2 = glob.glob("../data/" + regulatoryElement + "/" + cancerType2 + "/*_clustered.bed")

			if len(superEnhancerFile2) < 1:
				continue
			if regulatoryElement == 'enhancers':
				superEnhancerData2 = getEnhancersFromFile(superEnhancerFile2[0])
			else:
				superEnhancerData2 = readBedFile(superEnhancerFile2[0])

			malDictTissue2 = getChromosomeMal()
			for row in superEnhancerData2:
				chrom = row[0]
				
				malDictTissue2[chrom][int(row[1])] = 1

			#compute the differences.
			for chrom in malDictTissue1:
				
				difference = 0
				for pos in malDictTissue1[chrom]:
					if pos not in malDictTissue2[chrom]:
						difference += 1

				for pos in malDictTissue2[chrom]:
					if pos not in malDictTissue1[chrom]:
						difference += 1
				
				#difference = np.abs(malDictTissue1[chrom] - malDictTissue2[chrom])

				#differenceMatrix[cancerTypes[cancerType1]][cancerTypes[cancerType2]] += np.sum(difference)
				differenceMatrix[cancerTypes[cancerType1]][cancerTypes[cancerType2]] += difference
	
	print(differenceMatrix)
					
	#plot heatmap to show results visually
	fig =plt.figure(figsize=(15,10))
	
	data = pd.DataFrame(differenceMatrix)
	g=sns.heatmap(data,annot=False,square=True, linewidths=0.5,
				  yticklabels=list(cancerTypes.keys()), xticklabels=list(cancerTypes.keys()),
				  cmap="OrRd")
	
	
	plt.tight_layout()
	
	plt.savefig(finalOutDir + '/featureCorrelation_' + regulatoryElement + '.svg')
	#plt.show()


