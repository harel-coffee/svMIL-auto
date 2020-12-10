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

import matplotlib
matplotlib.use('Agg')

#this code depends on the input parser from the linking part. This is quick and dirty, is there a better solution?
sys.path.insert(1, 'linkSVsGenes/')

from inputParser import InputParser
from neighborhoodDefiner import NeighborhoodDefiner

finalOutDir = 'output/featureComparisons/'
if not os.path.exists(finalOutDir):
	os.makedirs(finalOutDir)

def readBedFile(elementFile):

	possibleChromosomes, coordinates = getChromosomeMal()

	elementData = []
	with open(elementFile, "r") as f:
		lineCount = 0
		for line in f:
			line = line.strip()
			splitLine = line.split("\t")

			chrom = splitLine[0]

			if chrom[0] != 'c':
				chrom = 'chr' + chrom
				
			if chrom == 'chr23':
				chrom = 'chrX'
			elif chrom == 'chr24':
				chrom = 'chrY'
			elif chrom not in possibleChromosomes:
				continue

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

regulatoryElements = ['eQTLs', 'tads', 'enhancers', 'superEnhancers', 'h3k4me1', 'h3k27ac',
					  'h3k27me3', 'h3k4me3', 'ctcf', 'dnase', 'rnapol']
#regulatoryElements = ['tads']
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
	return malDict, coordinates



#Option 1: make vectors and correlate, or compute absolute distance
optimalDistances = dict()
regulatoryElementsMap = dict()
ind = 0
for regulatoryElement in regulatoryElements:
	regulatoryElementsMap[regulatoryElement] = ind
	ind += 1

for regulatoryElement in regulatoryElements:
	print('Comparing: ', regulatoryElement)

	differenceMatrix = np.zeros([len(cancerTypes), len(cancerTypes)])
	for cancerType1 in cancerTypes:
		print('cancer type: ', cancerType1)

		if cancerType1 not in optimalDistances:
			optimalDistances[cancerType1] = dict()
		if regulatoryElement not in optimalDistances[cancerType1]:
			#optimalDistances[cancerType1][regulatoryElement] = [float("inf"), None]
			optimalDistances[cancerType1][regulatoryElement] = []

		#load the regulatory elements
		superEnhancerFile = glob.glob("../data/" + regulatoryElement + "/" + cancerType1 + "/*")

		#select for the clustered one
		if len(superEnhancerFile) > 1:
			superEnhancerFile = glob.glob("../data/" + regulatoryElement + "/" + cancerType1 + "/*_" + regulatoryElement + ".bed")

		if len(superEnhancerFile) < 1:
			continue

		if regulatoryElement == 'enhancers':
			superEnhancerData = getEnhancersFromFile(superEnhancerFile[0])
		else:
			superEnhancerData = readBedFile(superEnhancerFile[0])


		malDictTissue1, coordinates = getChromosomeMal()
		for row in superEnhancerData:
			chrom = row[0]

			chromCoordinates = coordinates[chrom]


			#for pos in range(int(row[1]), int(row[2])):
			#	malDictTissue1[chrom][pos] = 1
			malDictTissue1[chrom][int(row[1])] = 1
			malDictTissue1[chrom][int(row[2])] = 1

		for cancerType2 in cancerTypes:
			if cancerType1 == cancerType2:
				continue

			superEnhancerFile2 = glob.glob("../data/" + regulatoryElement + "/" + cancerType2 + "/*")

			#select for the clustered one
			if len(superEnhancerFile2) > 1:
				superEnhancerFile2 = glob.glob("../data/" + regulatoryElement + "/" + cancerType2 + "/*_" + regulatoryElement + ".bed")

			if len(superEnhancerFile2) < 1:
				continue
			if regulatoryElement == 'enhancers':
				superEnhancerData2 = getEnhancersFromFile(superEnhancerFile2[0])
			else:
				superEnhancerData2 = readBedFile(superEnhancerFile2[0])

			malDictTissue2, coordinates = getChromosomeMal()
			for row in superEnhancerData2:
				chrom = row[0]
				
				#for pos in range(int(row[1]), int(row[2])):
				#	malDictTissue2[chrom][pos] = 1

				malDictTissue2[chrom][int(row[1])] = 1
				malDictTissue2[chrom][int(row[2])] = 1


				
			#for TADs, compute overlap
			if regulatoryElement not in regulatoryElements:
				
				similarityPerChrom = dict()
				for row in superEnhancerData:
					chrom = row[0]
					
					if chrom not in similarityPerChrom:
						similarityPerChrom[chrom] = 0
	
					chrMatches = superEnhancerData2[superEnhancerData2[:,0] == chrom]
	
					startMatches = row[1] <= chrMatches[:,2]
					endMatches = row[2] >= chrMatches[:,1]
	
					allMatches = chrMatches[startMatches * endMatches]
	
					for match in allMatches:
						if match[2] > row[2] and match[1] < row[1]:
							similarityPerChrom[chrom] += row[2] - row[1]
						elif row[2] > match[2] and row[1] < match[1]:
							similarityPerChrom[chrom] += match[2] - match[1]
						elif row[1] < match[1] and row[2] < match[2]:
							similarityPerChrom[chrom] += row[2] - match[1]
						elif row[2] > match[2] and row[1] > match[1]:
							similarityPerChrom[chrom] += match[2] - row[1]
							
				#subtract chr length from similarity
				difference = 0
				for chrom in coordinates:

					if chrom not in similarityPerChrom:
						continue

					end = coordinates[chrom]
					difference += end - similarityPerChrom[chrom]
				

				differenceMatrix[cancerTypes[cancerType1]][cancerTypes[cancerType2]] += difference
			else:
				#compute the differences.
				for chrom in malDictTissue1:

					difference = 0
					for pos in malDictTissue1[chrom]:
						if chrom not in malDictTissue2:
							difference += 1
						if pos not in malDictTissue2[chrom]:
							difference += 1

					for pos in malDictTissue2[chrom]:
						if chrom not in malDictTissue1:
							difference += 1
						if pos not in malDictTissue1[chrom]:
							difference += 1

					#difference = np.abs(malDictTissue1[chrom] - malDictTissue2[chrom])

					#differenceMatrix[cancerTypes[cancerType1]][cancerTypes[cancerType2]] += np.sum(difference)
					differenceMatrix[cancerTypes[cancerType1]][cancerTypes[cancerType2]] += difference

			#save the result if smaller than for a previous cancer type.
			#so for cancer type 1, and this regulatory element, is it smaller than
			#what we have saved so far?
			if differenceMatrix[cancerTypes[cancerType1]][cancerTypes[cancerType2]] == 0:
				continue #do not select missing data or itself
			#if differenceMatrix[cancerTypes[cancerType1]][cancerTypes[cancerType2]] < optimalDistances[cancerType1][regulatoryElement][0]:

			#	optimalDistances[cancerType1][regulatoryElement][0] = differenceMatrix[cancerTypes[cancerType1]][cancerTypes[cancerType2]]
			#	optimalDistances[cancerType1][regulatoryElement][1] = cancerType2
			optimalDistances[cancerType1][regulatoryElement].append([int(differenceMatrix[cancerTypes[cancerType1]][cancerTypes[cancerType2]]), cancerType2])
					
	
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


for cancerType in optimalDistances:
	print('cancerType: ', cancerType)
	
	
	optimalDistanceMatrix = np.zeros([len(cancerTypes), len(regulatoryElements)])
	#make sure that we do not select itself as minimum
	#optimalDistanceMatrix[cancerTypes[cancerType]][regulatoryElementsMap[regulatoryElement]] = len(cancerTypes)-1
	for regulatoryElement in optimalDistances[cancerType]:
		print(regulatoryElement, ": ", optimalDistances[cancerType][regulatoryElement])

		if len(optimalDistances[cancerType][regulatoryElement]) < 1:
			continue #if regulatory data is missing for this cancer type.

		#sort the distances for this regulatory element
		distanceArray = np.array(optimalDistances[cancerType][regulatoryElement], dtype='object')

		sortedDistanceArray = distanceArray[np.argsort(distanceArray[:,0])][::-1]
		for arrayInd in range(0, sortedDistanceArray.shape[0]):
			cancerType2 = sortedDistanceArray[arrayInd,1]

			optimalDistanceMatrix[cancerTypes[cancerType2]][regulatoryElementsMap[regulatoryElement]] = arrayInd


		# cancerType2 = optimalDistances[cancerType][regulatoryElement][1]
		# if cancerType2 is None:
		# 	continue
		# optimalDistanceMatrix[cancerTypes[cancerType2]][regulatoryElementsMap[regulatoryElement]] = 1

	#determine the overall best cancer type to use:


	cancerTypeAvg = np.mean(optimalDistanceMatrix, axis=1)
	bestInd = np.argmax(cancerTypeAvg)
	for arrayInd in range(0, len(cancerTypes)):
		cancerType2 = list(cancerTypes.keys())[arrayInd]
		print(cancerType2, ": ", cancerTypeAvg[arrayInd])

	print('Best match: ', list(cancerTypes.keys())[bestInd])

	#merge the ranking with the distance matrix for plotting

	#get the ranks
	temp = cancerTypeAvg.argsort()
	ranks = np.empty_like(temp)
	ranks[temp] = np.arange(len(cancerTypeAvg))

	optimalDistanceMatrixRanked = np.zeros([len(cancerTypes), len(regulatoryElements)+1])
	optimalDistanceMatrixRanked[0:optimalDistanceMatrix.shape[0], 0:optimalDistanceMatrix.shape[1]] = optimalDistanceMatrix
	optimalDistanceMatrixRanked[:,optimalDistanceMatrix.shape[1]] = ranks


	#plot heatmap per cancer type
	fig =plt.figure(figsize=(5,5))

	#data = pd.DataFrame(optimalDistanceMatrix)
	data = pd.DataFrame(optimalDistanceMatrixRanked)
	g=sns.heatmap(data,annot=False,square=True, linewidths=0.5,
				  xticklabels=list(regulatoryElementsMap.keys()) + ['Ranking'], yticklabels=list(cancerTypes.keys()),
				  cmap="OrRd")

	plt.title(cancerType + ', best match: ' + list(cancerTypes.keys())[bestInd])
	plt.tight_layout()

	plt.savefig(finalOutDir + '/optimum_' + cancerType + '.svg')
	#plt.show()






