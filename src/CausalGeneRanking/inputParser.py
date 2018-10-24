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
import re

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
					
				#Skip anything that is not breast cancer for now. From here is the easiest way, saves time in processing as well
				if cancerType != "breast":
					continue
				
				svTypeIndex = header.index("sv_type")
				svType = splitLine[svTypeIndex]
				
				#Check if the SV type matches deletions
				match = re.search("deletion", svType, re.IGNORECASE)
				if match is None: #only focus on deletions for now
					continue
				
				
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
				
				#Make sure to switch the positions here as well
				#Some positions are swapped
				if int(e2) < int(e1):
					tmpE1 = e1
					e1 = e2
					e2 = tmpE1
					tmpS1 = s1
					s1 = s2
					s2 = tmpS1
				
				#Sometimes only the end is swapped.
				if int(e2) < int(s2):
					tmpS2 = s2
					s2 = e2
					e2 = tmpS2
					
				if int(e1) < int(s1):
					tmpS1 = s1
					s1 = e1
					e1 = tmpS1
				
				
				svObject = SV('chr' + chr1, s1, e1, 'chr' + chr2, s2, e2, sampleName, cancerType)
				#chr 1, start, end, chr2, start2, end2
				variantsList.append(['chr' + chr1, s1, e1, 'chr' + chr2, s2, e2, cancerType, sampleName, svObject])
		
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
				
				gene = Gene(geneSymbol, "chr" + chromosome, int(start), int(end)) #Keep in objects for easy access of properties related to the neighborhood of the gene
				
				cosmicGenes.append(["chr" + chromosome, int(start), int(end), gene])
				
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
				
				
				#Obtain the name, chromosome and positions of the gene. 
				
				geneID = splitLine[0]
				
				chrom = splitLine[1]

				start = splitLine[2]
				end = splitLine[3]
				
				geneObj = Gene(geneID, "chr" + chrom, int(start), int(end))
				
				nonCausalGeneList.append(["chr" + chrom, int(start), int(end), geneObj])
				
		nonCausalGenes = np.array(nonCausalGeneList, dtype="object")
	
		
		return nonCausalGenes 
		
			
			
		
		
			
	