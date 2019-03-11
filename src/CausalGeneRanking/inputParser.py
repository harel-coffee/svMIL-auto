"""
	Class intended to read all files that are provided to the program.
	There will be functions to:
	- Read SV files
	- Read SNV files
	- Read gene files
	- Read all the individual data types (e.g. eQTLs, enhancers, Hi-C data)
	
	Each of these input files will require a specific format with at least some required fields. 

"""
from gene import Gene
from sv import SV
from snv import SNV
from tad import TAD
from element import Element
import numpy as np
from random import randint
import re
import settings

class InputParser:
	
	def getSVsFromFile(self, svFile, typeFilter):
		"""
			Parse the SV data from the SV input file.
			
			svFile: (string) location of the SV file to read
			typeFilter: (string) meant to filter by which types of SVs to include in the model, but does not work. 
			
			return
			regions: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
		"""
		
			
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
				
				#Now extract the chromosome, start and end (there are multiple, 2 for an SV)
				chr1Index = header.index("chr1")
				s1Index = header.index("s1")
				e1Index = header.index("e1")
				o1Index = header.index("o1")
		
				chr2Index = header.index("chr2")
				s2Index = header.index("s2")
				e2Index = header.index("e2")
				o2Index = header.index("o2")
		
				cancerTypeIndex = header.index("cancer_type")
				sampleNameIndex = header.index("sample_name")
				
				cancerType = splitLine[cancerTypeIndex]
				sampleName = splitLine[sampleNameIndex]
				
				# #Dirty workaround to make sure that the cancer type names are the same, we only focus on 1 type for this intial run
				if cancerType == "breast/gastric":
					cancerType = "breast"
					
				#Skip anything that is not breast cancer for now. From here is the easiest way, saves time in processing as well
				if cancerType != settings.general['cancerType']:
					continue
				
				svTypeIndex = header.index("sv_type")
				svType = splitLine[svTypeIndex]
				
				# doubleTypeMatch = re.search("_", svType, re.IGNORECASE)
				# if doubleTypeMatch is not None:
				# 	continue
				#
				# 
				if svType != "invers":
					continue
				
				# #only keep the main 4 types
				# if svType != "del" and svType != "invers" and svType != "tandem_dup":
				# 	
				# 	interChrTypeMatch = re.search("chr", svType, re.IGNORECASE)
				# 	transTypeMatch = re.search("trans", svType, re.IGNORECASE)
				# 	rangeTypeMatch = re.search("range", svType, re.IGNORECASE)
				# 	if interChrTypeMatch is None and transTypeMatch is None and rangeTypeMatch is None:
				# 		continue
				# 	
					
				
				# if svType != "tandem_dup":
				# 	continue
				
				# interChrTypeMatch = re.search("chr", svType, re.IGNORECASE)
				# transTypeMatch = re.search("trans", svType, re.IGNORECASE)
				# rangeTypeMatch = re.search("range", svType, re.IGNORECASE)
				# if interChrTypeMatch is None and transTypeMatch is None and rangeTypeMatch is None:
				# 	continue
			
				
				# if typeFilter != "all":
				# 	#Check if the SV type matches deletions
				# 	match = re.search("deletion", svType, re.IGNORECASE)
				# 	if match is None: #only focus on deletions for now
				# 		continue
				# 
				
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
				o1 = splitLine[o1Index]
				o2 = splitLine[o2Index]
				
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
				
				
				svObject = SV('chr' + chr1, s1, e1, o1, 'chr' + chr2, s2, e2, o2, sampleName, cancerType, svType)
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
		"""
			TO DO:
			- Re-add SNVs to the model and fully document
		
		"""
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
				
				chromosome = "chr" + colonSplitPosition[0]
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
				
				gene = Gene(geneSymbol, "chr" + chromosome, int(start), int(end)) #Keep in objects for easy access of properties related to the neighborhood of the gene
				
				cosmicGenes.append(["chr" + chromosome, int(start), int(end), gene])
				
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
			geneObj = gene[3]
			causalGeneDict[geneObj.name] = 1			
		
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
				
				geneObj = Gene(geneID, chrom, int(start), int(end))
				
				if geneID not in causalGeneDict:
				
					nonCausalGeneList.append([chrom, int(start), int(end), geneObj])
				
		nonCausalGenes = np.array(nonCausalGeneList, dtype="object")
	
		
		return nonCausalGenes 
	
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
				
				
				TADObject = TAD(splitLine[0], int(splitLine[1]), int(splitLine[2]))
				
				#chr, start, end
				tadData.append([splitLine[0], int(splitLine[1]), int(splitLine[2]), TADObject])
		
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
	def getEQTLsFromFile(self, eQTLFile, genes, neighborhoodDefiner):
		"""
			Read the eQTLs from the file.
			
			eQTLFile: (string) Location of the eQTL file on disk.
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			neighborhoodDefiner: (object) neighborhoodDefiner class. Used to assign the elements to genes to determine losses later on. 
			
			return
			eQTLs: (numpy array) array with eQTL elements. chr, start, end, ElementObject
			
		"""
		#Filter the eQTLs that do not have a match
		geneDict = dict()
		
		for gene in genes:
			if gene not in geneDict:
				geneDict[gene.name] = gene
		
		
		eQTLs = []
		with open(eQTLFile, 'rb') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				
				if splitLine[3] not in geneDict:
					continue
				
				#Add the chr notation for uniformity. 		
				chrMatch = re.search("chr", splitLine[0], re.IGNORECASE)
				chrName = ""
				if chrMatch is None:
					chrName = "chr" + splitLine[0]
				else:
					chrName = splitLine[0]
				# eQTLObject = Element(chrName, int(splitLine[1]), int(splitLine[2])) #chr, start, end
				# eQTLObject.type = 'eQTL' #set the correct type
				# #The mapping information is in the file, so we can already do it here
				#This function belongs more to the neighborhood definer, so we use the function from there. 
				#neighborhoodDefiner.mapElementsToGenes(eQTLObject, geneDict, splitLine[3])
				eQTL = [chrName, int(splitLine[1]), int(splitLine[2]), "eQTL", splitLine[3]]
				neighborhoodDefiner.mapElementsToGenes(eQTL, geneDict, splitLine[3])

				#eQTLs.append([chrName, int(splitLine[1]), int(splitLine[2]), eQTLObject, "eQTL"]) #Keep the eQTL information raw as well for quick overlapping.
				eQTLs.append(eQTL) #Keep the eQTL information raw as well for quick overlapping. 
		
		
		return np.array(eQTLs, dtype='object')
	
	def getLncRNAsFromFile(self, lncRNAFile):
		"""
			TO DO:
			- Currently unused function
			- Needs to either be tested or removed
		
		"""
		
		
		
		
		lncRNAs = []
		with open(lncRNAFile, 'rb') as f:
			
			for line in f:
				line = line.strip()
				splitLine = line.split("\t")
				
				#Quick and dirty
				lncRNAObject = EQTL(splitLine[0], int(splitLine[1]), int(splitLine[2]))
				
				lncRNAs.append([splitLine[0], int(splitLine[1]), int(splitLine[2]), lncRNAObject])					
		
		return np.array(lncRNAs, dtype="object")
	
	
	def getEnhancersFromFile(self, enhancerFile, genes, neighborhoodDefiner):
		"""
			Read the enhancers from the file.
			
			enhancerFile: (string) Location of the enhancer file on disk.
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			neighborhoodDefiner: (object) neighborhoodDefiner class. Used to assign the elements to genes to determine losses later on. 
			
			return
			enhancers: (numpy array) array with enhancer elements. chr, start, end, ElementObject
		
		"""
		
		
		geneDict = dict()
		
		for gene in genes:
			if gene not in geneDict:
				geneDict[gene.name] = gene
		
		
		enhancers = []
		with open(enhancerFile, 'rb') as f:
			
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
				
				if geneName not in geneDict:
					continue
				
				
				# elementObject = Element(chrName, start, end)
				# elementObject.type = "enhancer"
				
				#The mapping information is in the file, so we can already do it here
				
						
				element = [chrName, start, end, "enhancer", geneName]
				neighborhoodDefiner.mapElementsToGenes(element, geneDict, geneName)
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
		with open(promoterFile, 'rb') as f:
			
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
				neighborhoodDefiner.mapElementsToGenes(promoter, geneDict, finalGeneName)
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
		with open(cpgFile, 'rb') as f:
			
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
		with open(tfFile, 'rb') as f:
			print tfFile
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
		with open(histoneFile, 'rb') as f:
			
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
		with open(dnaseIFile, 'rb') as f:
			
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
			
	