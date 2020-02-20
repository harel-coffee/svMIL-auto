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
from gene import Gene
from sv import SV
# from snv import SNV
from tad import TAD
from element import Element
import numpy as np
from random import randint
import glob
import re
import settings

class InputParser:
	
	def getSVsFromFile(self, svFile, typeFilter, excludedSVs): #excluded SVs no longer used here!!!
		"""
			Parse the SV data from the SV input file.
			
			svFile: (string) location of the SV file to read
			typeFilter: (string) meant to filter by which types of SVs to include in the model, but does not work. 
			
			return
			regions: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
		"""
		
			
		variantsList = []
		
		with open(svFile, 'r') as f:
			
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
				
				
			
				
				if settings.general['cancerType'] == "germline":
					
					if svType != "deletion" and svType != "inversion" and svType != "duplication":
						continue
			
				if settings.general['cancerType'] == "BRCA":
					
					if svType != "del" and svType != "invers" and svType != "tandem_dup" and svType != "DEL" and svType != "INV" and svType != "DUP":
						
						interChrTypeMatch = re.search("chr", svType, re.IGNORECASE)
						transTypeMatch = re.search("trans", svType, re.IGNORECASE)
						rangeTypeMatch = re.search("range", svType, re.IGNORECASE)
						itxTypeMatch = re.search("ITX", svType, re.IGNORECASE)
						ctxTypeMatch = re.search("CTX", svType, re.IGNORECASE)
						if interChrTypeMatch is None and transTypeMatch is None and rangeTypeMatch is None and itxTypeMatch is None and ctxTypeMatch is None:
							continue
					#if svType != 'del':
					# 	continue
					
				
			
		
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
				e1 = int(float(splitLine[e1Index]))
				s2 = int(splitLine[s2Index])
				e2 = int(float(splitLine[e2Index]))
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
				
				#remove the underscore from SV type
				svType = svType.replace("_", ".")
	
				if list(chr1)[0] == "c": #add the chr notation only when it is not already there
					svObject = SV(chr1, s1, e1, o1, chr2, s2, e2, o2, sampleName, cancerType, svType)
					
					#Check if this SV needs to be exluded
					svStr = chr1 + "_" + str(s1) + "_" + str(e1) + "_" + chr2 + "_" + str(s2) + "_" + str(e2) + "_" + sampleName + '_' + svType
					# if svStr not in excludedSVs:
					# 	continue

					variantsList.append([chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject])
					
				else:
					svObject = SV('chr' + chr1, s1, e1, o1, 'chr' + chr2, s2, e2, o2, sampleName, cancerType, svType)
					
					#Check if this SV needs to be exluded
					svStr = 'chr' + chr1 + "_" + str(s1) + "_" + str(e1) + "_" + 'chr' + chr2 + "_" + str(s2) + "_" + str(e2) + "_" + sampleName + '_' + svType

					# if svStr not in excludedSVs:
					# 	continue

					variantsList.append(['chr' + chr1, s1, e1, 'chr' + chr2, s2, e2, cancerType, sampleName, svObject])
				#chr 1, start, end, chr2, start2, end2
				
				
			
		regions = np.array(variantsList, dtype='object')
		
		return regions
	
	def getSVsFromFile_hmf(self, svDir): 
		"""
			Parse the SV data from the SV input file.
			
			svFile: (string) location of the SV file to read
			typeFilter: (string) meant to filter by which types of SVs to include in the model, but does not work. 
			
			return
			regions: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
		"""
		
		#Get all parsed and annotated SV type files from the main dir
		
		vcfs = glob.glob(svDir + '/**/*.svTypes.passed', recursive=True)
		
		variantsList = []
		addedVariants = [] #check if based on pairs no duplicates are added. 
		count = 0
		for vcf in vcfs:
			print(vcf)
			
			#get the samplename from the vcf
			sampleName = re.search('.*\/([A-Z\d]+)\.', vcf).group(1)
			
			#open the .gz file
			with open(vcf, 'r') as inF:
				
				for line in inF:

					if re.search('^#', line): #skip header
						continue
					
					#skip the SV if it did not pass.
					splitLine = line.split("\t")
					filterInfo = splitLine[6]
					if filterInfo != 'PASS':
						continue
					
					chr1 = splitLine[0]
					pos1 = int(splitLine[1])
					pos2Info = splitLine[4]
					
					#match the end position and orientation. if there is no orientation info, this is an insertion, which we can skip.
					if not re.search(':', pos2Info):
						continue

					if re.match('[A-Z]*\[.*\:\d+\[$', pos2Info):
						o1 = '+'
						o2 = '-'
					elif re.match('[A-Z]*\].*\:\d+\]$', pos2Info):
						o1 = '-'
						o2 = '+'
					elif re.match('^\].*\:\d+\][A-Z]*', pos2Info):
						o1 = '+'
						o2 = '+'
					elif re.match('^\[.*\:\d+\[[A-Z]*', pos2Info):
						o1 = '-'
						o2 = '-'
					else:
						print('unmatched: ', pos2Info)
						print(line)
						exit()
					
					#get the chr2 information
					chr2 = re.search('[\[\]]+(.*):(\d+).*', pos2Info).group(1)
					pos2 = int(re.search('.*\:(\d+).*', pos2Info).group(1))
					
					infoField = splitLine[7]
					splitInfoField = infoField.split(";")
					svType = ''
					for field in splitInfoField:
						
						splitField = field.split("=")
						if splitField[0] == 'SIMPLE_TYPE':
							svType = splitField[1]
					
					#skip SV types that we do not consider
					if svType not in ['DEL', 'DUP', 'ITX', 'INV']:
						continue
					
					#default positions
					s1 = pos1
					e1 = pos1
					s2 = pos2
					e2 = pos2
					orderedChr1 = chr1
					orderedChr2 = chr2
					
					#switch chromosomes if necessary
					if chr1 != chr2:
						if chr1 == 'Y' and chr2 == 'X':
							orderedChr1 = chr2
							orderedChr2 = chr1
						if (chr1 == 'X' or chr1 == 'Y' or chr1 == 'MT') and (chr2 != 'X' and chr2 != 'Y' and chr2 != 'MT'):
							orderedChr1 = chr2
							orderedChr2 = chr1
						if (chr1 != 'X' and chr1 != 'Y' and chr1 != 'MT') and (chr2 != 'X' and chr2 != 'Y' and chr2 != 'MT'):
							if int(chr1) > int(chr2):
								orderedChr1 = chr2
								orderedChr2 = chr1
						if (chr1 in ['X', 'Y', 'MT']) and (chr2 in ['X', 'Y', 'MT']): #order these as well
							if chr1 == 'Y' and chr2 == 'X':
								orderedChr1 = chr2
								orderedChr2 = chr1
							if chr1 == 'MT' and chr2 in ['X', 'Y']:
								orderedChr1 = chr2
								orderedChr2 = chr1
	
						
						#always switch the coordinates as well if chromosomes are switched.
						if orderedChr1 == chr2:
							s1 = pos2
							e1 = pos2
							s2  = pos1
							e2 = pos1	
					
					else: #if the chr are the same but the positions are reversed, change these as well. 
						if pos2 < pos1:
							s1 = pos2
							e1 = pos2
							s2  = pos1
							e2 = pos1	

					finalChr1 = 'chr' + orderedChr1
					finalChr2 = 'chr' + orderedChr2
					svObject = SV(finalChr1, s1, e1, o1, finalChr2, s2, e2, o2, sampleName, settings.general['cancerType'], svType)
					
					#check if this SV is already added. That may happen with the pair IDs. 
					svStr = finalChr1 + "_" + str(s1) + "_" + str(e1) + "_" + finalChr2 + "_" + str(s2) + "_" + str(e2) + "_" + sampleName
					
					
					
					if svStr in addedVariants:
						continue

					variantsList.append([finalChr1, s1, e1, finalChr2, s2, e2, settings.general['cancerType'], sampleName, svObject])
					addedVariants.append(svStr)

				#Get the required information for this SV:
				
				#chr1, s1, e1, o1, chr2, s2, e2, o2, samplename, cancer type, svType
	
		svs = np.array(variantsList, dtype='object')
		return svs
		#np.savetxt('hmf_svs.txt', svs, fmt='%s', delimiter='\t')
		
		# print(np.unique(svs[:,7]).shape)
		# 
		# svTypes = dict()
		# patientDistr = dict()
		# for sv in svs:
		# 	
		# 	if sv[8].svType not in svTypes:
		# 		svTypes[sv[8].svType] = 0
		# 	svTypes[sv[8].svType] += 1
		# 	
		# 	if sv[7] not in patientDistr:
		# 		patientDistr[sv[7]] = 0
		# 	patientDistr[sv[7]] += 1
		# 	
		# print(svTypes)
		# print(patientDistr)
		# import matplotlib.pyplot as plt
		# plt.hist(list(patientDistr.values()))
		# plt.show()
		# 	
		# exit()
		

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
		
		with open(snvFile, 'r') as f:
			
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
				gene.cosmic = 1
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
	
	def getCTCFSitesFromFile(self, ctcfFile):
		
		ctcfSites = []
		with open(ctcfFile, 'r') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				
				ctcf = [splitLine[0], int(splitLine[1]), int(splitLine[2]), "ctcf", None, float(splitLine[4])]

				ctcfSites.append(ctcf) #Keep the eQTL information raw as well for quick overlapping. 
		
		
		return np.array(ctcfSites, dtype='object')
		
		

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
		with open(eQTLFile, 'r') as f:
			
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
				eQTL = [chrName, int(splitLine[1]), int(splitLine[2]), "eQTL", splitLine[3], None]
				neighborhoodDefiner.mapElementsToGenes(eQTL, geneDict, splitLine[3])

				#eQTLs.append([chrName, int(splitLine[1]), int(splitLine[2]), eQTLObject, "eQTL"]) #Keep the eQTL information raw as well for quick overlapping.
				eQTLs.append(eQTL) #Keep the eQTL information raw as well for quick overlapping. 
		
		
		return np.array(eQTLs, dtype='object')
	
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
				#splitInteraction = interaction.split("_")
				
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
				
				score = float(splitInteraction[2])
				
				#Get the gene name
				splitGeneInfo = splitInteraction[1].split("$")
				geneName = splitGeneInfo[1]

				
				if geneName not in geneDict:
					continue
				
				
				# elementObject = Element(chrName, start, end)
				# elementObject.type = "enhancer"
				
				#The mapping information is in the file, so we can already do it here
				
						
				element = [chrName, start, end, "enhancer", geneName, score]
				neighborhoodDefiner.mapElementsToGenes(element, geneDict, geneName)
				enhancers.append(element)
		
		
		return np.array(enhancers, dtype='object')
		
		#unlinked enhancers for the time being
		enhancers = []
		with open(enhancerFile, 'r') as f:
			
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				
				element = [splitLine[0], int(splitLine[1]), int(splitLine[2]), 'enhancer', None]
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
				
				promoter = [chrName, start, end, "promoter", finalGeneName, None]
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
				
				cpgIsland = [chrName, start, end, "cpg", None, None] #None because it is not associated with a gene
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

				tfs.append([chrName, start, end, "tf", None, None])
		
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
				
			
				histoneMarks.append([chrName, start, end, histoneType, None, float(splitLine[4])])
		
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
				
			
				dnaseISites.append([chrName, start, end, "dnaseI", None, float(splitLine[4])])
		
		return np.array(dnaseISites, dtype='object')	
	
	def getChromHmmFromFile(self, chromHmmFile):
		
		
		chromHmmSites = []
		with open(chromHmmFile, 'r') as f:
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				chrName = splitLine[0]
				start = int(splitLine[1])
				end = int(splitLine[2])
				chromState = splitLine[3]
				
				chromHmmSites.append([chrName, start, end, chromState, None, None])
				
		
		return np.array(chromHmmSites, dtype='object')

	def getRnaPolFromFile(self, rnaPolFile):

		rnaPolSites = []
		with open(rnaPolFile, 'r') as f:
			lineCount = 0
			for line in f:
				if lineCount < 1:
					lineCount += 1
					continue
				
				line = line.strip()
				splitLine = line.split("\t")
				
				chrName = splitLine[0]
				start = int(splitLine[1])
				end = int(splitLine[2])
				
				rnaPolSites.append([chrName, start, end, 'rnaPol', None, float(splitLine[4])])
				
		
		return np.array(rnaPolSites, dtype='object')

	def getSuperEnhancersFromFile(self, superEnhancerFile):
	
		superEnhancers = []
		with open(superEnhancerFile, 'r') as f:
			
			for line in f:
				line = line.strip()
				splitLine = line.split("\t")
				
				superEnhancers.append([splitLine[0], int(splitLine[1]), int(splitLine[2]), 'superEnhancer', None, None])					
		
		return np.array(superEnhancers, dtype="object")

	def getMethylationFromFile(self, methylationFile, genes):
		"""
			Add the methylation information of the genes to the sv-gene pair altered elements for the given genes. This will be used in the MIL feature vector. 
		"""
		
		#get a list of the genes that have altered elements. These are the ones that we need to get from the methylation file
		affectedGenes = []
		for gene in genes:
			if len(gene[3].alteredElements) > 0:
				affectedGenes.append([gene[3].name, gene[3]])
		affectedGenes = np.array(affectedGenes, dtype='object')
		
		methylation = dict()
		patients = dict() #make a lookup for patients, so that we can easily calculate the line number at which the beta value for that patient/gene pair is located. 
		with open(methylationFile, 'r') as f:
			lineCount = 0
			for line in f:
				
				line = line.strip()
				splitLine = line.split("\t")
				if lineCount < 1:
					patientsList = splitLine[1:] #skip hybridization ref
					for patientInd in range(0, len(patientsList)):
						splitPatientName = patientsList[patientInd].split("-")
						patientID = 'brca' + splitPatientName[2] #we need a way around this later, because now it is specific for brca... due to the bad naming in the SV file.
						if patientID not in patients: #make sure to keep this unique, and always start with the first position of that patient in the line for later lookup.
							patients[patientID] = patientInd + 1 # +1 to avoid the hybrid ref position
					
					lineCount += 1
					continue
				if lineCount < 2:
					lineCount += 1
					continue
		
				#this information is repeated across the entire line.
				geneName = splitLine[2]
				if geneName not in methylation:
					methylation[geneName] = [] #only take the first patient as an example
				methylation[geneName].append(splitLine[1])
				
				
				if geneName not in affectedGenes[:,0]:
					continue #skip this line if the gene is not affected
				
				#get the affected gene.
				affectedGene = affectedGenes[affectedGenes[:,0] == geneName][0][1]
				
				#get the index of the patients that affect this gene.
				for sv in affectedGene.alteredElements:
					splitSV = sv.split("_")
					patientID = splitSV[6]
					#use the lookup to determine the first position in this line that we find data for this patient
					if patientID not in patients: #some patients appear to not have methylation data
						betaValue = 0
					else:
						patientLineInd = patients[patientID]
						
						#then the beta value should be the first entry. Gene is the second, then the chromosome, then the coordinates of the methylation. 
						betaValue = splitLine[patientLineInd]
					
						if betaValue == 'NA':
							betaValue = 0 #for now
					
					#Add the beta value to the altered elements
					for element in affectedGene.alteredElements[sv]:
						affectedGene.alteredElements[sv][element] += [float(betaValue)]

				# 
				# #there are 4 values per patient
				# #skipping the hybrid ref, we can start at pos 1
				# lineInd = 1
				# for patientInd in range(0, len(patients)):
				# 	patient = patients[patientInd]
				# 	#encode as: chromosome, coordinate, beta value, locus name, patient
				# 	#print([splitLine[lineInd+2], splitLine[lineInd+3], splitLine[lineInd], splitLine[lineInd+1], patient])
				# 	#methylation.append(['chr' + splitLine[lineInd+2], int(splitLine[lineInd+3]), splitLine[lineInd], splitLine[lineInd+1], patient])
				# 	if splitLine[lineInd+1] not in methylation:
				# 		methylation[splitLine[lineInd+1]] = []
				# 	methylation[splitLine[lineInd+1]].append(splitLine[lineInd])
				# 	
				# 	lineInd += 4
		for gene in methylation:
			if len(methylation[gene]) > 1:
				print(gene, methylation[gene])
		exit()		
		#return np.array(methylation, dtype='object')