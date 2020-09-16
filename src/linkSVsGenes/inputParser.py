"""
	Class intended to read all files that are provided to the program.
	There are functions to:
	- Read SV files
	- Read gene files
	- Read all the individual data types (e.g. eQTLs, enhancers, Hi-C data)
	
	Each of these input files will require a specific format with at least some required fields. 

"""
from __future__ import absolute_import
from __future__ import print_function
from gene import Gene
from sv import SV
from tad import TAD
import numpy as np
from random import randint
import glob
import re
import settings
import gzip

class InputParser:
	
	def readSVFile_pcawg(self, svFile, sampleName):
		
		addedVariants = [] # keep track of already read SVs to avoid duplicates
		variantsList = [] #all final variants in the file
		#open the .gz file
		with gzip.open(svFile, 'rb') as inF:
			
			for line in inF:
				line = line.strip().decode('utf-8')

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
					if splitField[0] == 'SVCLASS':
						
						if re.search('DEL', splitField[1]):
							svType = 'DEL'
						elif re.search('DUP', splitField[1]):
							svType = 'DUP'
						elif re.search('INV', splitField[1]):
							svType = 'INV'
						elif re.search('TRA', splitField[1]):
							svType = 'ITX'
						else:
							print('unknown: ', field)
							exit()

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
		#svs = np.array(variantsList, dtype='object')
		return variantsList
	
	def getSVsFromFile_pcawg(self, svDir):
		
		#get the cancer type from the settings
		cancerType = settings.general['cancerType']
		#read in the metadata file and get the right file identifiers
		metadataFile = settings.files['pcawgMetadata']
		
		#save the IDs of the patients with this cancer type
		cancerTypeIds = dict()
		with open(metadataFile, 'r') as inF:
			
			for line in inF:
				
				if re.search(cancerType, line):
					splitLine = line.split('\t')
					sampleId = splitLine[1]
					cancerTypeIds[sampleId] = 0
					
		
		#Then read the SV files for the right cancer type
		allSVs = []
		count = 0
		for sampleId in cancerTypeIds:
			
			#use glob to find the right file
			matchedFiles = glob.glob(svDir + '/' + sampleId + '*sv.vcf.gz')

			#if we don't have SVs for this sample, skip it. 
			if len(matchedFiles) < 1:
				print('Skipping ', sampleId, ' which has no SVs')
				continue

			#### here check if we have CNVs and SNVs for this sample. If these are not there,
			#we cannot say anything about the causality of SVs of this patient, so it is
			#better to skip it.
			cnvDir = settings.files['cnvDir']
			snvDir = settings.files['snvDir']
			matchedSNVFiles = glob.glob(snvDir + '/' + sampleId + '*.somatic.snv_mnv')
			matchedCNVFiles = glob.glob(cnvDir + '/' + sampleId + '*copyNumberEstimation*.somatic.cnv.vcf.gz')

			### also, not all CNVs are in the PCAWg directory, so use the ones from TCGA
			#where necessary. So we check for samples in there too.
			tcgaCNVFile = settings.files['tcgaCNVFile']
			tcgaCNVs = np.loadtxt(tcgaCNVFile, dtype='object')

			if len(matchedSNVFiles) < 1:
				print('Skipping ', sampleId, ' due to missing SNVs')
				continue

			if len(matchedCNVFiles) < 1 and sampleId not in tcgaCNVs[:,0]:
				print('Skipping ', sampleId, ' due to missing CNVs')
				continue
			
			#there should be just 1 file
			sampleSVFile = matchedFiles[0]
			print('Parsing SVs for file: ', sampleSVFile)
			
			#read in the SVs from this file
			sampleSVs = self.readSVFile_pcawg(sampleSVFile, sampleId)
			allSVs = allSVs + sampleSVs
			
		allSVs = np.array(allSVs, dtype='object')
		print(allSVs.shape)
		
		return allSVs
	
	
	def getSVsFromFile(self, svFile, typeFilter):
		"""
			Parse the SV data from the SV input file, specific for PCAWG/TCGA/germline data.
			
			svFile: (string) location of the SV file to read
			typeFilter: (string) meant to filter by which types of SVs to include in the model, is not used anymore.
			
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

				
				
				svTypeIndex = header.index("sv_type")
				svType = splitLine[svTypeIndex]
				
				#If the coordinates are missing on the second chromosome, we use the coordinates of the first chromosome unless chr 1 and chr 2 are different.
				if splitLine[chr1Index] == splitLine[chr2Index]:
					if splitLine[s2Index] == 'NaN':
						splitLine[s2Index] = int(splitLine[s1Index])
						
					if splitLine[e2Index] == 'NaN':
						splitLine[e2Index] = int(splitLine[e1Index])
				else:
					if splitLine[chr2Index] == 'NaN':
						continue # This line does not have correct chromosome 2 information so we don't know what this is
				
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

				#For TCGA/germline/PCAWG the SV naming is always different, so map everything back here to 1 format. 
				if svType == 'del':
					svType = 'DEL'
				elif svType == 'tandem_dup':
					svType = 'DUP'
				elif svType == 'invers':
					svType = 'INV'
				elif svType == 'transl_inter':
					svType = 'ITX'
				elif svType == 'INV':
					svType = 'INV'

				#germline
				elif svType == 'deletion':
					svType = 'DEL'
				elif svType == 'duplication':
					svType = 'DUP'
				elif svType == 'inversion':
					svType = 'INV'

				#also for pcawg data
				elif svType == 'DEL':
					svType = 'DEL'
				elif svType == 'DUP':
					svType = 'DUP'
				elif svType == 't2tINV':
					svType = 'INV'
				elif svType == 'h2hINV':
					svType = 'INV'
				elif svType == 'TRA':
					svType = 'ITX'

				else:
					print(svType)
					continue #unmatched SV type
	
				if list(chr1)[0] == "c": #add the chr notation only when it is not already there
					svObject = SV(chr1, s1, e1, o1, chr2, s2, e2, o2, sampleName, cancerType, svType)
					variantsList.append([chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject])
					
				else:
					svObject = SV('chr' + chr1, s1, e1, o1, 'chr' + chr2, s2, e2, o2, sampleName, cancerType, svType)
					variantsList.append(['chr' + chr1, s1, e1, 'chr' + chr2, s2, e2, cancerType, sampleName, svObject])

			
		regions = np.array(variantsList, dtype='object')
		return regions
	
	def getSVsFromFile_hmf(self, svDir): 
		"""
			Parse the SV data from the SV input file, specific for HMF
			
			svDir: (string) location of the SV vcf files to read
			
			return
			regions: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
		"""
		
		#Get all parsed and annotated SV type files from the main dir
		
		vcfs = glob.glob(svDir + '/**/*.svTypes.passed', recursive=True)
		
		variantsList = []
		addedVariants = [] #check if based on pairs no duplicates are added. 
		count = 0
		for vcf in vcfs:
			# count += 1
			# if count > 3:
			# 	continue

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

		svs = np.array(variantsList, dtype='object')

		return svs

	def getSNVsFromFile(self, snvFile):
		"""
			DEPRECATED

			Used to read SNV data from the TCGA. No longer used in the model at this stage.
		
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
		
		#If CRDs are provided, read those instead and treat them as TADs.
		if settings.general['crd'] == True:
			print('Reading CRDs instead of TADs and using them as if they were TADs')
			sortedTads = self.getCRDsFromFile(tadFile)
			return sortedTads

		
		#Read the TAD data into a list
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

	def getCRDsFromFile(self, crdFile):
		"""
			Get the CRDs from the provided CRD file. CRDs are similar to TADs, so we
			use them as if they were TADs.

			crdFile: (string) location of the crd file on disk

			return:
			sortedCRDs: (numpy array) array with the CRDs and their information. Sorted by chromosome & start position. chr, start, end, tadObject

		"""

		crdData = []
		with open(crdFile, "r") as f:
			for line in f:
				line = line.strip()
				splitLine = line.split("\t")

				#check for TRUE/FALSE to see if the entry is labeled as CRD
				crdLabel = splitLine[7]

				if crdLabel == 'FALSE':
					continue
				
				TADObject = TAD(splitLine[0], int(splitLine[1]), int(splitLine[2]))
				
				#chr, start, end
				crdData.append([splitLine[0], int(splitLine[1]), int(splitLine[2]), TADObject])
		
		#Also convert the other dataset to numpy
		crdData = np.array(crdData, dtype='object')
		
		#Make sure to sort the tads oer chromosome
		sortedCRDs = []
		chroms = np.unique(crdData[:,0])
		for chromosome in chroms:
			crdSubset = crdData[crdData[:,0] == chromosome]

			sortedSubset = crdSubset[crdSubset[:,1].argsort()]
			for crd in sortedSubset:

				sortedCRDs.append(crd)
			
		return np.array(sortedCRDs)

	def getCTCFSitesFromFile(self, ctcfFile):
		"""
			Read the CTCF sites from the provided CTCF file.

			ctcfFile: (string) location of the CTCF file on disk

			return:
			ctcfSites: (numpy array) array with the CTCF sites and their information. chr, start, end, elementType, object reference (unused), strength (peak intensity)

		"""
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

				ctcfSites.append(ctcf) 

		
		return np.array(ctcfSites, dtype='object')

	#Reading eQTL file
	def getEQTLsFromFile(self, eQTLFile, genes, neighborhoodDefiner):
		"""
			Read the eQTLs from the file.

			eQTLFile: (string) Location of the eQTL file on disk.
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			neighborhoodDefiner: (object) neighborhoodDefiner class. Used to assign the elements to genes to determine losses later on.
			
			return
			eQTLs: (numpy array) array with eQTL elements. chr, start, end, elementType, gene object, strength (None for eQTLs)
			
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
				
				
				eQTL = [chrName, int(splitLine[1]), int(splitLine[2]), "eQTL", splitLine[3], None]
				#This function belongs more to the neighborhood definer, so we use the function from there.
				#We use it because eQTLs, enhancers and promoters are linked to genes in the data, so a loss is only a loss if it was linked to the gene originally. 
				neighborhoodDefiner.mapElementsToGenes(eQTL, geneDict, splitLine[3])

				eQTLs.append(eQTL) #Keep the eQTL information in numpy array instead of objects for quick overlapping. 
		
		return np.array(eQTLs, dtype='object')
	
	def getEnhancersFromFile(self, enhancerFile, genes, neighborhoodDefiner):
		"""
			Read the enhancers from the file.
			
			enhancerFile: (string) Location of the enhancer file on disk.
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			neighborhoodDefiner: (object) neighborhoodDefiner class. Used to assign the elements to genes to determine losses later on. 
			
			return
			enhancers: (numpy array) array with enhancer elements. chr, start, end, gene object, strength
		
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
				
				score = float(splitInteraction[2]) #score is here the enhancer confidence. 
				
				#Get the gene name
				splitGeneInfo = splitInteraction[1].split("$")
				geneName = splitGeneInfo[1]

				
				if geneName not in geneDict:
					continue
				
				
				element = [chrName, start, end, "enhancer", geneName, score]
				#Map enhancers to the genes that these are linked to in the data. 
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
			promoters: (numpy array) array with promoter elements. chr, start, end, gene object, strength (None)
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
				
				if finalGeneName not in geneDict: #if this gene is not in our causal/non-causal gene list, we can skip it
					continue
				
				
				promoter = [chrName, start, end, "promoter", finalGeneName, None]
				#Map promoters to the genes that these belong to. 
				neighborhoodDefiner.mapElementsToGenes(promoter, geneDict, finalGeneName)
				promoters.append(promoter)  
		
		return np.array(promoters, dtype='object')	

	def getCpgIslandsFromFile(self, cpgFile):
		"""
			Read the CpG islands from the file.
			
			cpgFile: (string) Location of the CpG file on disk.
			
			
			return
			cpgIslands: (numpy array) array with CpG elements. chr, start, end, geneObject (None), strength (None)
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
				
				cpgIsland = [chrName, start, end, "cpg", None, None] #None because it is not associated with a gene
				cpgIslands.append(cpgIsland) 
		
		return np.array(cpgIslands, dtype='object')	

	def getTranscriptionFactorsFromFile(self, tfFile):
		"""
			Read the transcription factors from the file.
			
			tfFile: (string) Location of the transcription factor file on disk.
			
			return
			tfs: (numpy array) array with TF elements. chr, start, end, geneObject (None), strength (None)
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
				
				tfs.append([chrName, start, end, "tf", None, None])
		
		return np.array(tfs, dtype='object')	

	def getHiCInteractionsFromFile(self, interactionsFile):
		"""
			Read the Hi-C interactions from the file. To speed this up, the interactions file should already be linked to TADs to skip an additional step in the tool. (see preprocessing.sh)
			
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
			Read the histone marks from the file. The histone marks are across multiple files in the same format, so we provide the type of histone as well. 
			
			histoneFile: (string) Location of the histone marks file on disk.
			
			return
			histoneMarks: (numpy array) array with histone mark elements. chr, start, end, geneObject (None), strength (peak intensity)
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
			dnaseISites: (numpy array) array with DNAse I sites elements. chr, start, end, geneObject (None), strength (peak intensity)
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
		"""
			Read the ChromHMM states from the file.

			chromHmmFile: (string) Location of the chromHMM file on disk.

			return
			chromHmmSites: (numpy array) array with chromHmm states elements. chr, start, end, geneObject (None), strength (None)
		"""
		
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
		"""
			Read the RNA pol II sites from the file.

			rnaPolFile: (string) Location of the RNA pol II sites file on disk.

			return
			rnaPolSites: (numpy array) array with RNA pol elements. chr, start, end, geneObject (None), strength (peak intensity)
		"""

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
		"""
			Read the super enhancers from the file.

			superEnhancerFile: (string) Location of the super enhancers file on disk.

			return
			superEnhancers: (numpy array) array with RNA pol elements. chr, start, end, geneObject (None), strength (None)
		"""
	
		superEnhancers = []
		with open(superEnhancerFile, 'r') as f:
			
			for line in f:
				line = line.strip()
				splitLine = line.split("\t")
				
				superEnhancers.append([splitLine[0], int(splitLine[1]), int(splitLine[2]), 'superEnhancer', None, None])					
		
		return np.array(superEnhancers, dtype="object")
