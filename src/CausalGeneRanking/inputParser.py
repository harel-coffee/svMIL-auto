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
import numpy as np

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
	
	def readSNVFile(self, snvFile):
		
		1+1
	
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
		
			
	