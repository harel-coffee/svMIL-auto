"""
	Class intended to read all files that are provided to the program.
	There will be functions to:
	- Read SV files
	- Read SNV files
	- Read gene files
	- Read all the individual data types (e.g. eQTLs, enhancers, Hi-C data)
	
	Each of these input files will require a specific format with at least some required fields. 

"""
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from random import randint
import re

class InputParser:

	
	def processSVGenePairs(self, pairsFile, pairsDegFile):
			
		svGenePairs = np.loadtxt(pairsFile, dtype='object')
		svGeneDegPairs = np.load(pairsDegFile, allow_pickle=True, encoding='latin1')
		
		#Split the data into positive (all DEG pairs) and negative (pairs that are not in the DEG set)
		
		print("splitting pairs")
		#Split the data into positive (all DEG pairs) and negative (pairs that are not in the DEG set)
		positivePairs = svGeneDegPairs[:,0]
		#negativePairs = np.setdiff1d(svGenePairs, svGeneDegPairs[:,0])
		
		#setdiff is very slow in python 3 somehow, so switch to ultrafast dictionaries
		posDict = dict()
		for pair in positivePairs:
			posDict[pair] = 1
		for pair in svGenePairs:
			if pair not in posDict:
				posDict[pair] = 0
		
		negativePairs = []
		for pair in posDict:
			if posDict[pair] == 0:
				negativePairs.append(pair)
		
		negativePairs = np.array(negativePairs, dtype='object')

		#First focus only on intrachromosomal SVs.
		intraPairsPositive = []
		for pair in positivePairs:
			splitPair = pair.split("_")
			
			chr1 = splitPair[1]
			chr2 = splitPair[4]
			
			if chr1 == chr2:
				intraPairsPositive.append(pair)
		intraPairsNegative = []
		for pair in negativePairs:
			splitPair = pair.split("_")
			
			chr1 = splitPair[1]
			chr2 = splitPair[4]
			
			if chr1 == chr2:
				intraPairsNegative.append(pair)
		
		positivePairs = np.array(intraPairsPositive)
		negativePairs = np.array(intraPairsNegative)
		
		#Randomly subsample the SV-pairs to be balanced with the DEG pairs
		np.random.seed(0)
		negativePairsSubsampled = np.random.choice(negativePairs, positivePairs.shape[0])
		
		print(positivePairs.shape)
		print(negativePairsSubsampled.shape)
		
		allPairs = np.concatenate((positivePairs, negativePairsSubsampled))
		labels = np.array([1]*positivePairs.shape[0] + [0]*negativePairsSubsampled.shape[0])
		
		#split the pairs to be able to sort them
		splitPairs = []
		for pair in allPairs:
			splitPair = pair.split("_")
			splitPairs.append([splitPair[0], splitPair[1], int(splitPair[2]), int(splitPair[3]), splitPair[4], int(splitPair[5]), int(splitPair[6]), splitPair[7]])
		
		splitPairs = np.array(splitPairs, dtype='object')
		sortedPairs = splitPairs[splitPairs[:,1].argsort()]
		sortedLabels = labels[splitPairs[:,1].argsort()] #TO DO here make labels as well with the right sorting

		return sortedPairs, sortedLabels


	#Read TAD data
	def getTADsFromFile(self, tadFile):
		"""
			Read the TADs from the provided TAD file. 
			
			tadFile: (string) location of the TAD file on disk
			
			return:
			sortedTads: (numpy array) array with the TADs and their information. Sorted by chromosome & start position. chr, start, end, tadObject
			
		"""
		
		
		#Read the gene list data into a list
		tadData = []
		with open(tadFile, "r") as f:
			lineCount = 0
			for line in f:
				if lineCount < 2: #skip header
					lineCount += 1
					continue
				line = line.strip()
				splitLine = line.split("\t")
				
				
				#chr, start, end
				tadData.append([splitLine[0], int(splitLine[1]), int(splitLine[2])])
		
		#Also convert the other dataset to numpy
		tadData = np.array(tadData, dtype='object')
		
		#Make sure to sort the tads oer chromosome
		sortedTads = []
		chroms = np.unique(tadData[:,0])
		for chromosome in chroms:
			tadSubset = tadData[tadData[:,0] == chromosome]
			
			sortedSubset = tadSubset[tadSubset[:,1].argsort()]
			for tad in sortedSubset:
			
				sortedTads.append(tad)
			
		
		return np.array(sortedTads)
	
	#Reading eQTL file
	def getEQTLsFromFile(self, eQTLFile):
		"""
			Read the eQTLs from the file.
			
			eQTLFile: (string) Location of the eQTL file on disk.
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			neighborhoodDefiner: (object) neighborhoodDefiner class. Used to assign the elements to genes to determine losses later on. 
			
			return
			eQTLs: (numpy array) array with eQTL elements. chr, start, end, ElementObject
			
		"""
		
		eQTLs = []
		with open(eQTLFile, 'r') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				
				#Add the chr notation for uniformity. 		
				chrMatch = re.search("chr", splitLine[0], re.IGNORECASE)
				chrName = ""
				if chrMatch is None:
					chrName = "chr" + splitLine[0]
				else:
					chrName = splitLine[0]

		
				eQTL = [chrName, int(splitLine[1]), int(splitLine[2]), "eQTL", splitLine[3]]

				#eQTLs.append([chrName, int(splitLine[1]), int(splitLine[2]), eQTLObject, "eQTL"]) #Keep the eQTL information raw as well for quick overlapping.
				eQTLs.append(eQTL) #Keep the eQTL information raw as well for quick overlapping. 
		
		
		return np.array(eQTLs, dtype='object')
	
	def getEnhancersFromFile(self, enhancerFile):
		"""
			Read the enhancers from the file.
			
			enhancerFile: (string) Location of the enhancer file on disk.
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			neighborhoodDefiner: (object) neighborhoodDefiner class. Used to assign the elements to genes to determine losses later on. 
			
			return
			enhancers: (numpy array) array with enhancer elements. chr, start, end, ElementObject
		
		"""
		
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
				splitInteraction = interaction.split("_") #first part is the enhancer, 2nd part the gene
				
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
				
				#Get the gene name
				splitGeneInfo = splitInteraction[1].split("$")
				geneName = splitGeneInfo[1]
				
				element = [chrName, start, end, "enhancer", geneName]
				enhancers.append(element)
		
		
		return np.array(enhancers, dtype='object')
	
	def getPromotersFromFile(self, promoterFile, genes, neighborhoodDefiner):
		"""
			Read the promoters from the file.
			
			promoterFile: (string) Location of the promoter file on disk.
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			neighborhoodDefiner: (object) neighborhoodDefiner class. Used to assign the elements to genes to determine losses later on. 
			
			return
			promoters: (numpy array) array with promoter elements. chr, start, end, ElementObject
		"""
		
		geneDict = dict()
		
		for gene in genes:
			if gene not in geneDict:
				geneDict[gene.name] = gene
		
		
		promoters = []
		with open(promoterFile, 'r') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				#Format of file:
				#chr		start	end	geneName_\d
				
				#Add the chr notation for uniformity. 		
				chrMatch = re.search("chr", splitLine[0], re.IGNORECASE)
				chrName = ""
				if chrMatch is None:
					chrName = "chr" + splitLine[0]
				else:
					chrName = splitLine[0]
					
				start = int(splitLine[1])
				end = int(splitLine[2])
				geneName = splitLine[3]
				splitGeneName = geneName.split("_")
				finalGeneName = splitGeneName[0] #only get the part before the underscore
				
				if finalGeneName not in geneDict:
					continue
				
				
				# elementObject = Element(chrName, start, end)
				# elementObject.type = "promoter"
				
				#The mapping information is in the file, so we can already do it here
				
				promoter = [chrName, start, end, "promoter", finalGeneName]
				promoters.append(promoter) #Keep the eQTL information raw as well for quick overlapping. 
		
		return np.array(promoters, dtype='object')	

	def getCpgIslandsFromFile(self, cpgFile):
		"""
			Read the CpG islands from the file.
			
			cpgFile: (string) Location of the CpG file on disk.
			
			
			return
			cpgIslands: (numpy array) array with CpG elements. chr, start, end, ElementObject
		"""
		
		cpgIslands = []
		with open(cpgFile, 'r') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				#Add the chr notation for uniformity. 		
				chrMatch = re.search("chr", splitLine[1], re.IGNORECASE)
				chrName = ""
				if chrMatch is None:
					chrName = "chr" + splitLine[1]
				else:
					chrName = splitLine[1]
					
				start = int(splitLine[2])
				end = int(splitLine[3])
				
				# elementObject = Element(chrName, start, end)
				# elementObject.type = "cpg"
				
				cpgIsland = [chrName, start, end, "cpg", None] #None because it is not associated with a gene
				cpgIslands.append(cpgIsland) #Keep the eQTL information raw as well for quick overlapping. 
		
		return np.array(cpgIslands, dtype='object')	

	def getTranscriptionFactorsFromFile(self, tfFile):
		"""
			Read the transcription factors from the file.
			
			tfFile: (string) Location of the transcription factor file on disk.
			
			
			return
			tfs: (numpy array) array with TF elements. chr, start, end, ElementObject
		"""
		
		tfs = []
		with open(tfFile, 'r') as f:
			print(tfFile)
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")

				#Add the chr notation for uniformity. 		
				chrMatch = re.search("chr", splitLine[0], re.IGNORECASE)
				chrName = ""
				if chrMatch is None:
					chrName = "chr" + splitLine[0]
				else:
					chrName = splitLine[0]
					
				start = int(splitLine[1])
				end = int(splitLine[2])
				
				# elementObject = Element(chrName, start, end)
				# elementObject.type = "tf"
				# 

				tfs.append([chrName, start, end, "tf", None])
		
		return np.array(tfs, dtype='object')	

	def getHiCInteractionsFromFile(self, interactionsFile):
		"""
			Read the Hi-C interactions from the file. To speed this up, the interactions file should already be linked to TADs to skip an additional step in the tool.
			
			interactionsFile: (string) Location of the Hi-C interactions file on disk.

			return
			interactions: (dictionary) the TAD ID as provided in the file is the key, the start positions of the interactions are the values. 

		"""
		
		#Obtain the interaction indices per TAD
		interactions = dict()
		
		with open(interactionsFile, 'r') as inF:
			
			for line in inF:
				
				line = line.strip()
				splitLine = line.split("\t")
				
				tad = splitLine[0]
				
				interactionPositions = splitLine[1].split(",")
				interactions[tad] = interactionPositions

		return interactions	
	
	def getHistoneMarksFromFile(self, histoneFile, histoneType):
		"""
			Read the histone marks from the file. The histone marks are across multiple files in the same format, so we provide the type of histone as well
			to assign it to the Element object as type. 
			
			histoneFile: (string) Location of the histone marks file on disk.
			
			
			return
			histoneMarks: (numpy array) array with histone mark elements. chr, start, end, ElementObject
		"""
		histoneMarks = []
		with open(histoneFile, 'r') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")

				#Add the chr notation for uniformity. 		
				chrMatch = re.search("chr", splitLine[0], re.IGNORECASE)
				chrName = ""
				if chrMatch is None:
					chrName = "chr" + splitLine[0]
				else:
					chrName = splitLine[0]
					
				start = int(splitLine[1])
				end = int(splitLine[2])
				
			
				histoneMarks.append([chrName, start, end, histoneType, None])
		
		return np.array(histoneMarks, dtype='object')	
	
	def getDNAseIFromFile(self, dnaseIFile):
		"""
			Read the DNAse I hypersensitivity sites from the file. 
			
			dnaseIFile: (string) Location of the DNAse I file on disk.

			return
			dnaseISites: (numpy array) array with DNAse I sites elements. chr, start, end, ElementObject
		"""
		
		dnaseISites = []
		with open(dnaseIFile, 'r') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")

				#Add the chr notation for uniformity. 		
				chrMatch = re.search("chr", splitLine[0], re.IGNORECASE)
				chrName = ""
				if chrMatch is None:
					chrName = "chr" + splitLine[0]
				else:
					chrName = splitLine[0]
					
				start = int(splitLine[1])
				end = int(splitLine[2])
				
			
				dnaseISites.append([chrName, start, end, "dnaseI", None])
		
		return np.array(dnaseISites, dtype='object')	
	
	def readCausalGeneFile(self, causalGeneFile):
		"""
			Read the COSMIC genes from the file.
			
			causalGeneFile: (string) location of the file with COSMIC genes.
			
			return
			cosmicGenesSorted: (numpy array) array with the genes and their information, lexographically sorted by chromosome. chr, start, end, geneObject
		"""
		
			
		cosmicGenes = [] 
		with open(causalGeneFile, 'r') as geneFile:
			
			lineCount = 0
			header = []
			for line in geneFile:
				splitLine = line.split("\t")
				#First extract the header and store it in the dictionary to remove dependency on the order of columns in the file
				if lineCount < 1:
		
					header = splitLine
					lineCount += 1
					continue
					
				#Obtain the gene name and gene position
				
				geneSymbolInd = header.index('Gene Symbol')
				genePositionInd = header.index('Genome Location')
				
				geneSymbol = splitLine[geneSymbolInd]
				genePosition = splitLine[genePositionInd]
				
				#Split the gene position into chr, start and end
				
				colonSplitPosition = genePosition.split(":")
				dashSplitPosition = colonSplitPosition[1].split("-")
				
				chromosome = colonSplitPosition[0]
				start = dashSplitPosition[0].replace('"',"") #apparently there are some problems with the data, sometimes there are random quotes in there
				end = dashSplitPosition[1].replace('"', "")
				
				if start == '' or end == '':
					continue
				
				cosmicGenes.append(["chr" + chromosome, int(start), int(end), geneSymbol])
				
		#Sort the genes
		cosmicGenes = np.array(cosmicGenes, dtype='object')
		
		sortedInd = np.lexsort((cosmicGenes[:,1], cosmicGenes[:,0])) #first sort by chromosome and then by start position. 
		cosmicGenesSorted = cosmicGenes[sortedInd]
	
		return cosmicGenesSorted
		
	def readNonCausalGeneFile(self, nonCausalGeneFile, causalGenes):
		"""
			Read the non-causal genes. These are filtered for genes that are in COSMIC to make sure that these do not overlap. 
			
			
			nonCausalGeneFile: (string) location of the file with non-causal (non-COSMIC) genes. 
			causalGenes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			
			return:
			nonCausalGenes: (numpy array) array with the non-causal genes and their information. chr, start, end, geneObject
			
		"""
		causalGeneDict = dict() #for filtering out genes that are already present in the causal gene list
		for gene in causalGenes:
			causalGeneDict[gene[3]] = 1			
		
		nonCausalGeneList = []
		nonCausalGeneNameDict = dict() #dictionary to keep the names of genes that are already in our list and don't need to be sampled again. 
		
		with open(nonCausalGeneFile, 'r') as geneFile:
			
			lineCount = 0
			for line in geneFile:
				line = line.strip()
				splitLine = line.split("\t")
				
				
				
				#Obtain the name, chromosome and positions of the gene. 
				
				geneID = splitLine[3]
				
				
				chrom = splitLine[0]

				start = splitLine[1]
				end = splitLine[2]
				
				if geneID not in causalGeneDict:
				
					nonCausalGeneList.append([chrom, int(start), int(end), geneID])
				
		nonCausalGenes = np.array(nonCausalGeneList, dtype="object")
	
		
		return nonCausalGenes 
	