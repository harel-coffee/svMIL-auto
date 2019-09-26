import settings
from inputParser import InputParser
import numpy as np



class FeatureMatrixDefiner:
	enhancerData = []
	tadData = []
	geneData = []
	
	"""
		Make a feature matrix for deep learning
	
	"""
	
	def defineFeatureMatrix(self, pairs, pairLabels, svGenePairsRules):
			
		bags = []
		labels = []
		currentChr = None
		pairInd = -1
		posBags = []
		negBags = []
		for pair in pairs:
			
			pairStr = "_".join([str(i) for i in pair])
			if pairStr not in svGenePairsRules[:,0]:
				continue
			#look in a 1 MB window around the SV of the pair.
			#find enhancers, genes and TADs that are inside this window
			pairInd += 1
			
			if pairInd > 5:
				continue
			
			sv = pair[1:len(pair)]
			
			instances = []
			
			#only get a new chr subset if the chr changes
			if sv[0] != currentChr:
				#Get the subset of eQTls on this chromosome
				chr1SubsetTAD = self.tadData[np.where(self.tadData[:,0] == sv[0])]
				print(chr1SubsetTAD.shape)
				currentChr = sv[0]
				
			startMatches = sv[1] <= chr1SubsetTAD[:,2]
			endMatches = sv[5] >= chr1SubsetTAD[:,1]
					
			tadMatches = chr1SubsetTAD[startMatches * endMatches]
			
			leftTAD = tadMatches[0]
			rightTAD = tadMatches[len(tadMatches)-1]
			
			if len(tadMatches) < 2: #skip deletions in tads for now
				continue
			
			
			#Get all the elements that are in these TADs
			leftElements = leftTAD[3].getElementsByRange(leftTAD[1], leftTAD[2])
			rightElements = rightTAD[3].getElementsByRange(rightTAD[1], rightTAD[2])
			
			for element in leftElements:
				instances.append([element[1], element[2]])
			for element in rightElements:
				instances.append([element[1], element[2]])
			
			
			
	
			#TADs
			
			instances.append([rightTAD[2], leftTAD[1]])
		
			#SV
			instances.append([sv[1], sv[5]])
			if pairLabels[pairInd] == 0:
				negBags.append(np.array(instances))
			else:
				posBags.append(np.array(instances))
			
		posBags = np.array(posBags)
		negBags = np.array(negBags)
		negBagsSubsampled = np.random.choice(negBags, posBags.shape[0])	
		
		bags = np.concatenate((posBags, negBagsSubsampled))
		instances = np.vstack(bags)
		labels = [1]*posBags.shape[0] + [0]*negBagsSubsampled.shape[0]
		
		from random import shuffle
		#shuffle(labels)
		
		labels = np.array(labels)
		
		return bags, instances, labels
		
	
	def setFeatureData(self):
		"""
			Get the feature data from all files. Re-use the original settings file from rule-based for now
		"""
		
		self.tadData = InputParser().getTADsFromFile(settings.files['tadFile'])
		#2. Get eQTLs from the eQTL file, and map eQTLs to genes. 
		eQTLData = [] #Keep empty by default in case we do not use eQTLs
		if settings.general['eQTLs'] == True: #Always check if eQTLs are enabled in the settings

			eQTLFile = settings.files['eQTLFile']
			print("getting eQTLs")
			self.eQTLData = InputParser().getEQTLsFromFile(eQTLFile, genes[:,3], self)

			# with open('eQTLData.pkl', 'wb') as h:
			# 	pkl.dump(eQTLData, h, protocol=pkl.HIGHEST_PROTOCOL)
			self.tadData = self.mapElementsToTads(eQTLData, tadData)
		
		
		#Get genes
		causalGenes = InputParser().readCausalGeneFile(settings.files['causalGenesFile'])
		nonCausalGenes = InputParser().readNonCausalGeneFile(settings.files['nonCausalGenesFile'], causalGenes) #In the same format as the causal genes.
		genes = np.concatenate((causalGenes, nonCausalGenes), axis=0)
		
		self.tadData = self.mapGenesToTads(genes, self.tadData) 
		
		#3. Get enhancers
		
		if settings.general['enhancers'] == True:
			print("getting enhancers")
			self.enhancerData = InputParser().getEnhancersFromFile(settings.files['enhancerFile'], genes[:,3], self)
			#Add the enhancers to TADs & genes as well	
			self.tadData = self.mapElementsToTads(self.enhancerData, self.tadData)	
		
		#4. Get promoters
		
		if settings.general['promoters'] == True:
			print("getting promoters")
			self.promoterData = InputParser().getPromotersFromFile(settings.files['promoterFile'], genes[:,3], self)
			
			#Add the promoters to the TADs
			self.tadData = self.mapElementsToTads(self.promoterData, self.tadData)
		
		#5. Get CpG islands
		if settings.general['cpgIslands'] == True:
			print("Getting cpg islands")
			self.cpgData = InputParser().getCpgIslandsFromFile(settings.files['cpgFile'])
		
			#Add the CpG sites to the TADs
			self.tadData = self.mapElementsToTads(self.cpgData, self.tadData)
		
		#6. Get Transcription factors
		if settings.general['transcriptionFactors'] == True:
			print("Getting transcription factors")

			self.tfData = InputParser().getTranscriptionFactorsFromFile(settings.files['tfFile'])
	
			#Add the CpG sites to the TADs
			self.tadData = self.mapElementsToTads(self.tfData, self.tadData)
		
		
		#7. Get Hi-C data
		if settings.general['hiC'] == True:
			print("Getting Hi-C data")
			self.hicData = InputParser().getHiCInteractionsFromFile(settings.files['hicFile'])
			
			#Map the interactions to TADs as elements
			self.tadData = self.mapInteractionsToTads(self.hicData, self.tadData)
		
		#8. Get histone marks
		if settings.general['histones'] == True:
			print("Getting histone marks")
			files = [settings.files['h3k9me3'], settings.files['h3k4me3'], settings.files['h3k27ac'], settings.files['h3k27me3'],
					 settings.files['h3k4me1'], settings.files['h3k36me3']]
			types = ['h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
			for histoneFileInd in range(0, len(files)):
				self.histoneData = InputParser().getHistoneMarksFromFile(files[histoneFileInd], types[histoneFileInd])
			
				#map the histone marks to the TADs
				self.tadData = self.mapElementsToTads(self.histoneData, self.tadData)
		
		#9. Get DNAse I hypersensitivty sites
		if settings.general['dnaseI'] == True:
			print("Getting DNAse I hypersensitivity sites")
			
			self.dnaseIData = InputParser().getDNAseIFromFile(settings.files['dnaseIFile'])
			
			self.tadData = self.mapElementsToTads(self.dnaseIData, self.tadData)
	
	def mapGenesToTads(self, genes, tadData):
		"""
			For computing effects of disruptions on genes, it is convenient to know which genes are located in which TADs.
			Find out for every TAD which genes are located inside of these, and map them to the TADs. 
		
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
		"""
		
		previousChromosome = 0
		for tad in tadData:
		
			if tad[0] != previousChromosome:
				previousChromosome = tad[0]
				geneChrSubset = genes[np.where(genes[:,0] == tad[0])] #TADs are intrachromosomal, so looking at 1 chromosome is sufficient
			
			
			#Find all genes that are within the TAD.
			#Because some genes are in two TADs, and the TAD definition is likely not entirely correct, we will count overlap of TADs with any part of the gene even if it is just 1 bp for now.
			#So the start of the gene must be before the end of the TAD, and the end of the gene after the start of the TAD. 
			startMatches = geneChrSubset[:,1] <= tad[2]
			endMatches = geneChrSubset[:,2] >= tad[1]
			
			allMatches = startMatches * endMatches
			matchingGenes = geneChrSubset[allMatches,:]
		
			#Add these genes to the TADs if any.
			if matchingGenes.shape[0] < 1:
				
				#If there are no genes, we can for now add the genes immediately to the left and right of the TAD.
				#With better data, this step can be removed, but since the TADs are so sparse, it is a temporary solution.
				
				#1. Compute the distance from either TAD boundary to all genes
				#We look at all genes to the left, so we can use the end of a gene to the start of a TAD, and the start of a gene to the end of a TAD. 
				#2. Get the gene with the smallest distance to each boundary
				
				startDistances = tad[1] - geneChrSubset[:,2]
				#Exclude all negative elements, these are on the right of the TAD while we want genes on the left of the TAD
				negativeElements = startDistances < 0

				startDistances[negativeElements] = float("inf") #skip the ones that are on the wrong side of the TAD boundary, so set them to inf and they will never be selected. 
				if startDistances[negativeElements].shape[0] == startDistances.shape[0]: #but if everything is inf, we skip this TAD.
					continue
								
				#Repeat but then for the other TAD end. 
				endDistances = geneChrSubset[:,1] - tad[2]
				#Exclude all negative elements
				negativeElements = endDistances < 0

				endDistances[negativeElements] = float("inf") #skip the ones that are on the wrong side of the TAD boundary, so set them to inf and they will never be selected. 
				if endDistances[negativeElements].shape[0] == endDistances.shape[0]: #but if everything is inf, we skip this TAD.
					continue
				
			else:
				tad[3].setGenes(matchingGenes[:,3]) #Add the eQTL objects to the TAD objects. 
	
		return tadData		
	
	def mapElementsToTads(self, elementData, tadData):
		"""
			Determine which eQTL interactions take place within which TAD.
			eQTL interactions will be mapped to TAD objects.
			
			elementData: (numpy array) array with elements. chr, start, end, ElementObject
			
		"""
		
		#Go through each TAD and find which elements are within the TAD.
		#List these elements as interaction objects in the tAD object.
		previousChromosome = 0
		for tad in tadData:
			if tad[0] != previousChromosome:
				previousChromosome = tad[0]
				elementChrSubset = elementData[np.where(elementData[:,0] == tad[0])] #Get all elements that are on this chromosome. All eQTLs are CIS, so intrachromosomal is fine for now

			#Find the eQTLs where the start and end of the eQTL are within the start and end of the TAD. 
			startMatches = elementChrSubset[:,1] >= tad[1]
			endMatches = elementChrSubset[:,2] <= tad[2]
			
			allMatches = startMatches * endMatches
			matchingElements = elementChrSubset[allMatches,:]

			#Add these elements to the TADs if any.
			if matchingElements.shape[0] < 1:
				
				continue
			else:
				
				#tad[3].addElements(matchingElements[:,3]) #Add the elements objects to the TAD objects.
				tad[3].addElements(matchingElements) #Add the elements objects to the TAD objects.
			
		return tadData
	
	def mapInteractionsToTads(self, interactionsByTad, tadData):
		"""
			
		"""
		
		hicOut = "../../../data/hic/hic.bed"
		with open(hicOut, 'w') as outF:
			
			print("mapping to tads")
			uniqueElements = dict()
			for tadStr in interactionsByTad:
				
				#1. Get the right TAD object
				splitTadStr = tadStr.split("_")
				
				
				tadChrMatch = tadData[:,0] == splitTadStr[0]
				tadStartMatch = tadData[:,1] == int(splitTadStr[1])
				tadEndMatch = tadData[:,2] == int(splitTadStr[2])
				
				matchingTad = tadData[tadChrMatch * tadStartMatch * tadEndMatch]
				
				if len(matchingTad) < 1: #This should not happen, but it does anyway? 
					continue
				
				matchingTad = matchingTad[0]
				
				#Assign the interactions to this TAD.
				#matchingTad[3].addElements(interactionsByTad[tadStr])
				interactionLines = []
				binSize = 5000
				
				for lineInd in range(0, len(interactionsByTad[tadStr])-1):
					line = interactionsByTad[tadStr][lineInd]
					# el = Element(splitTadStr[0], int(line), int(line)+binSize)
					# el.type = "hic"
					element = [splitTadStr[0], int(line), int(line)+binSize, "hic", None] #None because it is not associated with any gene
					interactionLines.append(element)
					outF.write(splitTadStr[0] + "\t" + str(line) + "\t" + str(int(line)+binSize) + "\n")
											
					#interactionLines.append([line, "hic"])
				# 	if line != "":
				# 		elementStr = splitTadStr[0] + "_" + str(line) + "_" + str(int(line)+binSize)
				# 		#if elementStr not in uniqueElements:
				# 		uniqueElements[elementStr] = [splitTadStr[0], int(line), int(line)+binSize]
				# 			#Add this element
				# 		#elementObject = Element(splitTadStr[0], int(line), int(line)+binSize)
				# 		#elementObject.type = "hic"
				# 		#matchingTad[3].addElements([elementObject])
				# 
				#matchingTad[3].addElements(np.array(interactionLines, dtype="object"))
				matchingTad[3].addElements(interactionLines)
		
		return tadData 
		