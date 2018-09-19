"""
	Class intended to read all files that are provided to the program.
	There will be functions to:
	- Read SV files
	- Read SNV files
	- Read gene files
	
	Each of these input files will require a specific format with at least some required fields. 

"""
from gene import Gene
from sv import SV
from snv import SNV
import numpy as np
from random import randint

class InputParser:
	
	def getSVsFromFile(self, svFile):
			
		variantsList = []
		
		with open(svFile, 'rb') as f:
			
			lineCount = 0
			header = []
			for line in f:
				line = line.strip()
				splitLine = line.split("\t")
				
				#First extract the header and store it in the dictionary to remove dependency on the order of columns in the file
				if lineCount < 1:
		
					header = splitLine
					lineCount += 1
					continue
				
				#Now extract the chromosome, start and end (there are multiple)
				chr1Index = header.index("chr1")
				s1Index = header.index("s1")
				e1Index = header.index("e1")
		
				chr2Index = header.index("chr2")
				s2Index = header.index("s2")
				e2Index = header.index("e2")
		
				cancerTypeIndex = header.index("cancer_type")
				sampleNameIndex = header.index("sample_name")
				
				cancerType = splitLine[cancerTypeIndex]
				sampleName = splitLine[sampleNameIndex]
				
				#Dirty workaround to make sure that the cancer type names are the same, we only focus on 1 type for this intial run
				if cancerType == "breast/gastric":
					cancerType = "breast"
				
				
				#if cancerType not in uniqueCancerTypes:
				#	uniqueCancerTypes.append(cancerType)
		
				#If the coordinates are missing on the second chromosome, we use the coordinates of the first chromosome unless chr 1 and chr 2 are different.
				if splitLine[chr1Index] == splitLine[chr2Index]:
					if splitLine[s2Index] == 'NaN':
						splitLine[s2Index] = int(splitLine[s1Index])
						
					if splitLine[e2Index] == 'NaN':
						splitLine[e2Index] = int(splitLine[e1Index])
				else:
					if splitLine[chr2Index] == 'NaN':
						continue # This line does not have correct chromosome 2 information (should we be skipping it?)
		
				s1 = int(splitLine[s1Index])
				e1 = int(splitLine[e1Index])
				s2 = int(splitLine[s2Index])
				e2 = int(splitLine[e2Index])
				chr2 = splitLine[chr2Index]
				
				chr1 = splitLine[chr1Index]
				
				svObject = SV(chr1, s1, e1, chr2, s2, e2, sampleName, cancerType)
				#chr 1, start, end, chr2, start2, end2
				variantsList.append([chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject])
		
		regions = np.array(variantsList, dtype='object')
		
		return regions
	

	#What do we want from the SNV file?
	#chromosome (combined with position in 'position' field)
	#position
	#cancer type (ID_tumour or Primary site) 
	#sample name (ID_SAMPLE)
	def getSNVsFromFile(self, snvFile):
		
		snvList = []
		
		with open(snvFile, 'rb') as f:
			
			lineCount = 0
			header = []
			for line in f:
				line = line.strip()
				splitLine = line.split("\t")
				
				#First extract the header and store it in the dictionary to remove dependency on the order of columns in the file
				if lineCount < 1:
		
					header = splitLine
					lineCount += 1
					continue
				
				#Now extract the chromosome, start and end (there are multiple)
				positionIndex = header.index("genome position")
				fullPosition = splitLine[positionIndex]
				
				#split the position to get the coordinates and the chromosome
				colonSplitPosition = fullPosition.split(":")
				dashSplitPosition = colonSplitPosition[1].split("-")
				
				chromosome = colonSplitPosition[0]
				start = dashSplitPosition[0]
				end = dashSplitPosition[1]
				
				cancerTypeIndex = header.index("Primary site") 
				sampleNameIndex = header.index("ID_SAMPLE")
				
				cancerType = splitLine[cancerTypeIndex]
				sampleName = splitLine[sampleNameIndex]
				
				#snvObject = SNV(chromosome, start, end, sampleName, cancerType)
				#snvList.append([chromosome, start, end, sampleName, cancerType, snvObject])
				snvList.append([chromosome, int(start), int(end), None, None, None, sampleName, cancerType]) #Add None because in the end we wish to treat variants the same way, but since objects take up memory it needs to be in a NP array but requires to be int he same column
				
		regions = np.array(snvList, dtype="object")
		
		return regions
	
	def readCausalGeneFile(self, causalGeneFile):
			
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
				
				gene = Gene(geneSymbol, chromosome, int(start), int(end)) #Keep in objects for easy access of properties related to the neighborhood of the gene
				
				cosmicGenes.append([chromosome, int(start), int(end), gene])
				
		#Sort the genes
		cosmicGenes = np.array(cosmicGenes, dtype='object')
		
		sortedInd = np.lexsort((cosmicGenes[:,1], cosmicGenes[:,0])) #first sort by chromosome and then by start position. 
		cosmicGenesSorted = cosmicGenes[sortedInd]
	
		return cosmicGenesSorted
		
	def readNonCausalGeneFile(self, nonCausalGeneFile, causalGenes):
		"""
			Read the non-causal genes. We currently take a random subset of these genes, and make sure that these are not overlapping with the COSMIC set.
			We first read all the genes, and then take a subset of the genes and make sure that these are not in cosmic. 
		"""
		
		nonCausalGeneList = []
		nonCausalGeneNameDict = dict() #dictionary to keep the names of genes that are already in our list and don't need to be sampled again. 
		
		with open(nonCausalGeneFile, 'r') as geneFile:
			
			lineCount = 0
			for line in geneFile:
				line = line.strip()
				splitLine = line.split("\t")
				if lineCount < 1:
					lineCount += 1
					continue
				
				#Obtain the name, chromosome and positions of the gene. 
				
				geneID = splitLine[4]
				chrom = splitLine[0]
				splitChrom = chrom.split("chr")
				finalChrom = splitChrom[1]

				start = splitLine[1]
				end = splitLine[2]
				
				geneObj = Gene(geneID, finalChrom, int(start), int(end))
				
				nonCausalGeneList.append([finalChrom, start, end, geneID, geneObj])
				
		nonCausalGenes = np.array(nonCausalGeneList)
		
		causalGeneNames = dict()
		for gene in causalGenes:
			causalGeneNames[gene[3].name] = 0
		
		
		randomNonCausalGenes = []
		#Then go through the list again, but take a random subset. Also make sure that the genes do not overlap with the COSMIC genes.
		noOfGenesToSample = 700 #there are not more than 700, apparently. 
		while len(nonCausalGeneNameDict) <= noOfGenesToSample:
		#while True:
			
			randomIndex = randint(0, nonCausalGenes.shape[0]-1)
			
			geneName = nonCausalGenes[randomIndex,3]
			if geneName not in nonCausalGeneNameDict:

				#Add the gene to our list if not in cosmic
				
				if geneName in causalGeneNames:
					
					randomNonCausalGenes.append([nonCausalGenes[randomIndex][0], nonCausalGenes[randomIndex][1], nonCausalGenes[randomIndex][2], nonCausalGenes[randomIndex][4]])
					nonCausalGeneNameDict[geneName] = 0
				
		
		randomNonCausalGenes = np.array(randomNonCausalGenes)
		
		return randomNonCausalGenes 
		
			
			
		
		
			
	