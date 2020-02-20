from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import json
import pickle as pkl
import re
import os.path
from os import listdir
from os.path import isfile, join
import glob
import gzip

from tad import TAD
from sv import SV
from gene import Gene
from element import Element
# from snv import SNV
from derivativeTADMaker import DerivativeTADMaker
# from genome import Genome
from genomicShuffler import GenomicShuffler
from inputParser import InputParser

import settings
from six.moves import range

class NeighborhoodDefiner:
	"""
		Class responsible for defining the neighborhood of genes.
		
		Currently, the neighborhood consists of:
		
		- All TADs in the genome
		- nearest TADs on the left and right of the gene
		- all eQTLs mapped to the gene and TADs (can be switched for lncRNAs)
		- all enhancers mapped to the gene and TADs
		- SVs (and SNVs) overlapping either the gene directly, or other elements in the neighborhood
		- Genes, with gains and losses of specific elements annotated as a result of SVs
		
	"""
	
	def checkIfSettingsAreSame(self):
		"""
			The purpose of this function is to check if the settings were updated since last run (currently only support for the data). If these did not change, we can load a pre-existing pkl to
			speed up the runs. 
			
			The settings are stored in a separate file. If the settings in the actual settings file are the same as in this cached file, return false.
			Otherwise, return true. 
			
			*** Currently, this function is not used becuase with pkl we exceed the recursion depth if we store the genes. I cannot change this depth. 
		"""
		
		#If the settings file is not on disk, create it
		settingsCacheFile = 'settings.cache'

		if os.path.isfile(settingsCacheFile) == False:
					
			with open(settingsCacheFile, 'w') as cacheFile:
				 json.dump(settings.general, cacheFile)
		
			#In this case, return true because the settings were just initialized and so we have no neighborhood yet
			return True
		else: #In this case, the file exists, so we can start checking if the settings are the same
			with open(settingsCacheFile, 'r') as cacheFile:
				cachedSettings = json.load(cacheFile)
				matchedSettings = True #if even 1 is different, the data should be re-loaded
				for setting in cachedSettings:
					if cachedSettings[setting] != settings.general[setting]:
						matchedSettings = False
						break
				
				return matchedSettings

	def __init__(self, genes, svData, snvData, mode, excludedSVs):
		"""
			Initialize the neighborhood defining. This involves gathering all required data types, mapping these to TADs/genes, and associating the effects of SVs to genes. 

			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject. Can be empty in SNV mode. 
			snvData: (numpy array) array with the SNVs and their information. chr, start, end, None, None, None, sampleName, cancerType. Can be empty in SV mode (NOT WORKING ATM)
			mode: (string) SV, SNV or SV+SNV.
		"""
		
		#1. Get TADs from the TAD file, and then map TADs to genes (left/right TAD).
		tadData = []
		if settings.general['tads'] == True:
			tadFile = settings.files['tadFile']
			
			print("Getting TADs")
			tadData = InputParser().getTADsFromFile(tadFile)
			
			#Filter the SVs out that overlap more than 10 TADs. (temporarily)
			print("original number of svs:", svData.shape)
			
			
			# filteredSvData = []
			# types = []
			# for sv in svData:
			# 	if sv[8].svType == "del" or sv[8].svType == "invers" or sv[8].svType == "tandem_dup": #only do the threshold thing for dels, invers and dups
			# 		tadSubset = tadData[tadData[:,0] == sv[0]]
			# 	
			# 		startMatches = sv[1] < tadSubset[:,1]
			# 		endMatches = sv[5] > tadSubset[:,2]
			# 		
			# 		allMatches = startMatches * endMatches
			# 		overlappedTads = tadSubset[allMatches]
			# 		if len(overlappedTads) > 2:
			# 			continue
			# 		filteredSvData.append(sv)
			# 	else:
			# 		filteredSvData.append(sv)
			# 	types.append(sv[8].svType)
			# 	
			# print np.unique(types)	
			# svData = np.array(filteredSvData, dtype="object")
			# print "number of svs: ", svData.shape
			# 
			# #temporarily write to a file
			# testOut = "unshuffledSVs.txt"
			# header = "chr1\ts1\te1\to1\tchr2\ts2\te2\to2\tsource\tsample_name\tsv_type\tcancer_type\n"
			# 
			# with open(testOut, 'w') as outF:
			# 	outF.write(header)
			# 	for sv in svData:
			# 		line = sv[0] + "\t" + str(sv[1]) + "\t" + str(sv[2]) + "\t" + sv[8].o1 + "\t" + sv[3] + "\t" + str(sv[4]) + "\t" + str(sv[5]) + "\t" + sv[8].o2 + "\t-\t" + sv[8].sampleName + "\t" + sv[8].svType + "\t" + sv[8].cancerType + "\n"
			# 		outF.write(line)
			# exit()
			
			#Plot the number of SVs that start or end in a TAD
			# svCount = []
			# for sv in svData:
			# 	
			# 	tadChrSubset = tadData[tadData[:,0] == sv[0]]
			# 	
			# 	#Number of TADs that are overlapped by an SV
			# 	startMatches = sv[1] < tadChrSubset[:,2]
			# 	endMatches = sv[5] > tadChrSubset[:,1]
			# 	
			# 	tadMatches = tadChrSubset[startMatches * endMatches]
			# 	# 
			# 	# startMatches = (sv[1] > tadChrSubset[:,1]) * (sv[1] < tadChrSubset[:,2])
			# 	# endMatches = (sv[5] > tadChrSubset[:,1]) * (sv[5] < tadChrSubset[:,2])
			# 	# 
			# 	# tadMatches = tadChrSubset[startMatches + endMatches]
			# 	
			# 	# if tadMatches.shape[0] not in svCount:
			# 	# 	svCount[tadMatches.shape[0]] = 0
			# 	# svCount[tadMatches.shape[0]] += 1	
			# 	#
			# 	if tadMatches.shape[0] < 50:
			# 		svCount.append(tadMatches.shape[0])
			# 
			# print svCount
			# 
			# import matplotlib.pyplot as plt
			# plt.hist(svCount)
			# #plt.bar(svCount.keys(), svCount.values())
			# plt.show()
			# exit()
				
			
			
			if settings.general['shuffleTads'] == True:
				#Shuffle the TADs. Assign random genomic positions to the TADs, but keep the same length. 
				genomicShuffler = GenomicShuffler()
				tadData = genomicShuffler.shuffleTADs(tadData)

			print("mapping TADs to genes")
			self.mapTADsToGenes(genes[:,3], tadData)
		
		#For every SV, find the TADs on the left and right
		
		#Compute how large this window is in total and report that (sum left TAD + right TAD for translocations)
		# windowSizes = []
		# for sv in svData:
		# 	
		# 	tadChr1Subset = tadData[tadData[:,0] == sv[0]]
		# 	tadChr2Subset = tadData[tadData[:,0] == sv[3]]
		# 
		# 	#Get the TAD that is overlapped by the left part of the SV
		# 	leftTad = tadChr1Subset[(tadChr1Subset[:,1] <= sv[1]) * (tadChr1Subset[:,2] >= sv[1])]
		# 	rightTad = tadChr2Subset[(tadChr2Subset[:,1] <= sv[4]) * (tadChr2Subset[:,2] >= sv[4])]
		# 	
		# 	if len(leftTad) < 1 or len(rightTad) < 1:
		# 		continue
		# 	
		# 	leftTad = leftTad[0]
		# 	rightTad = rightTad[0]
		# 	
		# 	if leftTad[3] == rightTad[3]:
		# 		continue
		# 	
		# 	leftTadSize = leftTad[2] - leftTad[1]
		# 	rightTadSize = rightTad[2] - rightTad[1]
		# 	
		# 	windowSizes.append(leftTadSize + rightTadSize)
		# 
		# 
		# print windowSizes
		# import matplotlib.pyplot as plt
		# 
		# plt.hist(windowSizes)
		# plt.show()
		# exit()
		
		
		#2. Get eQTLs from the eQTL file, and map eQTLs to genes. 
		eQTLData = [] #Keep empty by default in case we do not use eQTLs
		if settings.general['eQTLs'] == True: #Always check if eQTLs are enabled in the settings
			#Save the processed data, only needs to be done once
			
			# import os.path
			# 
			# if os.path.exists('eQTLData.pkl'):
			# 	print "loading eQTLs from pkl"
			# 	#Load the eqtls
			# 	with open('eQTLData.pkl', 'rb') as h:
			# 		eQTLData = pkl.load(h)
		
			eQTLFile = settings.files['eQTLFile']
			print("getting eQTLs")
			eQTLData = InputParser().getEQTLsFromFile(eQTLFile, genes[:,3], self)

			# with open('eQTLData.pkl', 'wb') as h:
			# 	pkl.dump(eQTLData, h, protocol=pkl.HIGHEST_PROTOCOL)
			tadData = self.mapElementsToTads(eQTLData, tadData)
		
		#Read the lncRNA data. If we enable lncRNAs, for now we switch that for eQTLs. Will be fixed in a later version. 
		if settings.general['lncRNA'] == True:
			lncRNAData = InputParser().getLncRNAsFromFile(settings.files['lncRNAFile'])
			eQTLData = lncRNAData 

		#Map the elements to the TADs, and map the genes to the TADs. 
		
		tadData = self.mapGenesToTads(genes, tadData) 
		
		#3. Get enhancers
		
		if settings.general['enhancers'] == True:
			print("getting enhancers")
			enhancerData = InputParser().getEnhancersFromFile(settings.files['enhancerFile'], genes[:,3], self)
			#Add the enhancers to TADs & genes as well	
			tadData = self.mapElementsToTads(enhancerData, tadData)	
		
		#4. Get promoters
		
		if settings.general['promoters'] == True:
			print("getting promoters")
			promoterData = InputParser().getPromotersFromFile(settings.files['promoterFile'], genes[:,3], self)
			
			#Add the promoters to the TADs
			tadData = self.mapElementsToTads(promoterData, tadData)
		
		#5. Get CpG islands
		if settings.general['cpgIslands'] == True:
			print("Getting cpg islands")
			cpgData = InputParser().getCpgIslandsFromFile(settings.files['cpgFile'])
		
			#Add the CpG sites to the TADs
			tadData = self.mapElementsToTads(cpgData, tadData)
		
		#6. Get Transcription factors
		if settings.general['transcriptionFactors'] == True:
			print("Getting transcription factors")

			tfData = InputParser().getTranscriptionFactorsFromFile(settings.files['tfFile'])
	
			#Add the CpG sites to the TADs
			tadData = self.mapElementsToTads(tfData, tadData)
		
		
		#7. Get Hi-C data
		if settings.general['hiC'] == True:
			print("Getting Hi-C data")
			hicData = InputParser().getHiCInteractionsFromFile(settings.files['hicFile'])
			
			#Map the interactions to TADs as elements
			tadData = self.mapInteractionsToTads(hicData, tadData)
		
		#8. Get histone marks
		if settings.general['histones'] == True:
			print("Getting histone marks")
			files = [settings.files['h3k9me3'], settings.files['h3k4me3'], settings.files['h3k27ac'], settings.files['h3k27me3'],
					 settings.files['h3k4me1'], settings.files['h3k36me3']]
			types = ['h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
			for histoneFileInd in range(0, len(files)):
				histoneData = InputParser().getHistoneMarksFromFile(files[histoneFileInd], types[histoneFileInd])
			
				#map the histone marks to the TADs
				tadData = self.mapElementsToTads(histoneData, tadData)
		
		#9. Get DNAse I hypersensitivty sites
		if settings.general['dnaseI'] == True:
			print("Getting DNAse I hypersensitivity sites")
			
			dnaseIData = InputParser().getDNAseIFromFile(settings.files['dnaseIFile'])
			
			tadData = self.mapElementsToTads(dnaseIData, tadData)
		
		#10. get chromHMM states
		if settings.general['chromHMM'] == True:
			print("Getting chromHMM states")
			chromHmmData = InputParser().getChromHmmFromFile(settings.files['chromHmmFile'])
			
			tadData = self.mapElementsToTads(chromHmmData, tadData)
		
		#11. get RNAPolII peaks
		if settings.general['rnaPol'] == True:
			print("Getting rnaPol binding sites")
			rnaPolData = InputParser().getRnaPolFromFile(settings.files['rnaPolFile'])
			
			tadData = self.mapElementsToTads(rnaPolData, tadData)
		
		#12. get super enhancers
		if settings.general['superEnhancers'] == True:
			print("Getting super enhancers")
			superEnhancerData = InputParser().getSuperEnhancersFromFile(settings.files['superEnhancerFile'])
			
			tadData = self.mapElementsToTads(superEnhancerData, tadData)
		
		if settings.general['ctcfSites'] == True:
			print("Getting ctcf sites")
			ctcfData = InputParser().getCTCFSitesFromFile(settings.files['ctcfFile'])
			
			tadData = self.mapElementsToTads(ctcfData, tadData)
			tadData = self.mapCTCFStrengthToTads(ctcfData, tadData)
		
		
		#3. Map SVs to all neighborhood elements
		if mode == "SV":
			print("Mapping SVs to the neighborhood")
			self.mapSVsToNeighborhood(genes, svData, tadData, excludedSVs)

		# if settings.general['snvs'] == False: #in this case, we want to filter out pairs that have SNV effects.
		# 	print('Mapping SNVs to genes')
		# 	self.mapSNVsToNeighborhood(genes, settings.files['snvDir'])
		# 	
		# if settings.general['cnvs'] == False: #in this case, we want to filter out pairs that have SNV effects.
		# 	print('Mapping CNVs to genes')
		# 	self.mapCNVsToNeighborhood(genes, settings.files['cnvDir'])
		# 		
			

		#Add the gene methylation to genes for MIL. Do this step here, because the file is huge and takes a long time to process.
		if settings.general['methylation'] == True:
			InputParser().getMethylationFromFile(settings.files['methylationFile'], genes)
		
			
	def mapCTCFStrengthToTads(self, ctcfData, tadData):
		
		#for each tad, check which ctcf site is overlapping the boundaries
		for tad in tadData:
			
			ctcfChrSubset = ctcfData[ctcfData[:,0] == tad[0]]
			
			startMatches = (tad[1] >= ctcfChrSubset[:,1]) * (tad[1] <= ctcfChrSubset[:,2])
			endMatches = (tad[2] >= ctcfChrSubset[:,1]) * (tad[2] <= ctcfChrSubset[:,2])
			
			if len(ctcfChrSubset[startMatches]) > 0:
				
				startSite = ctcfChrSubset[startMatches][0]
				tad[3].startStrength = startSite[5]
			
			if len(ctcfChrSubset[endMatches]) > 0:	
				endSite = ctcfChrSubset[endMatches][0]
				tad[3].endStrength = endSite[5]
			
		
		return tadData
		
	
	def mapTADsToGenes(self, genes, tadData):
		"""
			Adds the left and right TAD to each gene object.
			
			TO DO:
			- This function may be deprecated, I don't really use it anymore in the code but leave it here in case it is useful at a later point. 
			
		"""
		print(tadData)
		#For each gene, check which TADs are on the correct chromosome.
		
		#Then see which ones are directly on the left and directly on the right of the gene.
		previousChr = None
		for gene in genes:
			
			
			
			#1. Make a subset of TADs on the right chromosome. There should be only 1 chromosome
			
			if str(gene.chromosome) != previousChr:
				
				#Find the two subsets that match on both chromosomes. 
				matchingChrInd = tadData[:,0] == str(gene.chromosome)
				
				#It is not even necessary to make 2 lists if both chromosomes are the same, we could use a reference without re-allocating
				chrSubset = tadData[np.where(matchingChrInd)]
				
				#Make sure to update the previous chromosome when it changes
				previousChr = str(gene.chromosome)
				
			if np.size(chrSubset) < 1:
				continue #no need to compute the distance, there are no genes on these chromosomes
			
			#Within this subset, check which TADs are on the right and left of the current gene
			
			#TADs on the left have their start position before the gene. But their start must be closest to the gene to be the left TAD. 
			#TADs on the right have their start after the gene start. 
			
			leftTADs = chrSubset[np.where(chrSubset[:,1] <= gene.start)]
			rightTADs = chrSubset[np.where((chrSubset[:,1] >= gene.start) & (chrSubset[:,2] >= gene.end))] #make sure that the TAD is not entirely within the gene. 
			
			#Compute the distance to each of these TADs and take the TADs with the minimum distance
			if leftTADs.shape[0] > 0:
				leftTADsDistance = np.abs(leftTADs[:,1] - gene.start)
				nearestLeftTAD = leftTADs[np.argmin(leftTADsDistance),3]
				gene.setLeftTAD(nearestLeftTAD)
				
			else:
				gene.setLeftTAD(None)
			if rightTADs.shape[0] > 0:
				rightTADsDistance = np.abs(rightTADs[:,1] - gene.start)
				nearestRightTAD = rightTADs[np.argmin(rightTADsDistance),3]
				gene.setRightTAD(nearestRightTAD)
				
			else:
				gene.setRightTAD(None)

		
	def mapElementsToGenes(self, element, geneDict, geneSymbol):
		"""
			Map the right gene object to the eQTL object.
			
			element (object): Element object to which we add the gene that the element regulates.
			geneDict (dictionary): Dictionary with gene names and gene objects in the values.
			geneSymbol (string): Name of the gene that we want to map the element to
		
		"""
		
		geneDict[geneSymbol].addElement(element)
		#element.addGene(geneDict[geneSymbol])
		
	def mapInteractionsToGenes(self, genes, interactionData):
		"""
			Take Hi-C interactions as input (see getInteractionsFromFile for the format) and link these to the causal genes.
			
			TO DO:
			- Deprecated function, may be removed or re-added later. 
		
		"""
		
		#1. For each gene, find the non-coding regions that have an interaction with this gene. 
		
		previousChr = None
		for gene in genes:
			
			#Get the interactions bby the name of the gene
			interactionsSubset = interactionData[interactionData[:,3] == gene.name,5] #get the objects
			
			#add interaction objects to the genes
			gene.setInteractions(interactionsSubset)
	
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
		
		hicOut = "../../data/hic/hic.bed"
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
		
	def determineGainedInteractions(self, svData, tadData):
		"""
			Determine for each SV what the TADs are that the SV ends in on the left or the right. Determine which elements are gained as a result of a deletion and assign these to the genes in the respective TADs. 
			
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
			
			TO DO:
			- Move function to derivative TAD function
			
		"""
		
		#Loop through all SVs and see which ones have overlap with any TAD.
		#The SVs are sorted per cancer type, so taking chromosome subsets may not gain as much speed

		for sv in svData:
			
			#Only focus on deletions
			typeMatch = re.search("del", sv[8].svType, re.IGNORECASE)
			if typeMatch is None:
				continue
			
			##Here we need to take translocations into account as well.
			#If the SV is intrachromosomal, we should find the left TAD by chr1, s1 and the right TAD by e2.
			#if the SV is interchromosomal, we should find the left TAD by chr1, s1 and the right TAD by chr2 and e2. 
			#We may slightly overshoot TADs if there is a difference in the start and end coordinates, but this may do for now.

			if sv[0] == sv[3]: #intrachromosomal
				leftTadChr = sv[0]
				svStart = sv[1] #s1
				svEnd = sv[5] #e2
				
				#Now search for the matches on the left and right TAD by these coordinates
				#First get the right chromosome subset
				tadChromosomeSubset = tadData[np.where(tadData[:,0] == leftTadChr)]
				
				#The problem here is that SVs that are entirely within one TAD are also counted, which should not happen.
				#The SV should end in the left and right TAD, but not be the same TAD. 
				
				#First find all matches overlapping the right side of the TAD
				startMatches = svStart > tadChromosomeSubset[:,1]
				endMatches = svStart <= tadChromosomeSubset[:,2]
				allStartMatches = startMatches * endMatches
				
				
				#Then find the matches overlapping the left side of the TAD
				startMatches = svEnd < tadChromosomeSubset[:,2]
				endMatches = svEnd >= tadChromosomeSubset[:,1]
				allEndMatches = startMatches * endMatches
				
				#Every TAD that is matched for both the start and end should be excluded
				allMatches = np.logical_xor(allStartMatches, allEndMatches)
				
				overlappingTADsRight = tadChromosomeSubset[allStartMatches]
				overlappingTADsLeft = tadChromosomeSubset[allEndMatches]
				
				#Both ends should be in a TAD for this approach to work
				if overlappingTADsRight.shape[0] < 1:
					continue
				if overlappingTADsLeft.shape[0] < 1:
					continue

				#Make the total subset
				#overlappingTads = np.concatenate((overlappingTADsLeft, overlappingTADsRight), axis=0)
				overlappingTads = tadChromosomeSubset[allMatches]

				
				if overlappingTads.shape[0] < 1:
					continue
			else:
				continue
			
			#For the overlapping TADs, we can get the TAD on the right and the TAD on the left.
			#The SV could end on the boundary and then the desired result is not in our match set, so we skip these for now.
			
			#Get the TADs on the far left and far right of the SV.
			#These are the TADs that are still in tact.
			#We get the eQTL interactions that are in these two TADs (have been pre-mapped to the TAD that these take place in)
			#Then we find the genes that are present within the TAD (has also been pre-mapped to the TAD)
			#To each gene in one TAD, we add the number of gained interactions equal to the number of eQTL interactions in the other TAD.
				
			farLeftTad = overlappingTads[0] #This list is sorted
			farRightTad = overlappingTads[overlappingTads.shape[0]-1,:]

			
			#The genes in the far left TAD, only in the part that is not overlapped by the deletion, gain the elements that are not overlapped by the deletion
			#in the far right tad.
			remainingLeftGenes = farLeftTad[3].getGenesByRange(farLeftTad[1], sv[1])
			remainingLeftElements = farLeftTad[3].getElementsByRange(farLeftTad[1], sv[1])
			
			remainingRightGenes = farLeftTad[3].getGenesByRange(sv[5], farRightTad[2])
			remainingRightElements = farLeftTad[3].getElementsByRange(sv[5], farRightTad[2])
			
			for gene in remainingLeftGenes:
				
				if len(remainingRightElements) > 0:
					gene.addGainedElements(remainingRightElements, sv[7])
					gene.addGainedElementsSVs(remainingRightElements, sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[8].sampleName)
					
			
			for gene in remainingRightGenes:

				if len(remainingLeftElements) > 0:	
					gene.addGainedElements(remainingLeftElements, sv[7])
					gene.addGainedElementsSVs(remainingLeftElements, sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[8].sampleName)
		
		return 0
		
		
	def mapSVsToNeighborhood(self, genes, svData, tadData, excludedSVs):
		"""
			Take as input gene objects with the neighborhood pre-set, and search through the SVs to find which SVs overlap the genes, TADs and eQTLs in the gene neighborhood
		
			The TADs are also parsed as input, because when we compute gained interactions we search through other TADs that are on the other end of the SV breakpoint. This is the TAD from which we can
			expect interactions in our current TAD. 
		
			!!!!!! This function will need to be properly split into multiple functions for the different data types.
			TODO:
			- Split into multiple functions
			Move to derivative TAD maker for deletions as well. Here we now only look at losses, and gains are in a separate function. This should be moved to the derivative TAD class and do everything for all SVs at once. 
			- Move the gained interactions outside of the left/right TAD check. See comments below. 
		
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
		
		"""
		
		print("mapping SVs to genes")
		#Map SVs to genes that these overlap
		for sv in svData:
			
			#1. Get the genes on the right chromosome
			
			if sv[0] == sv[3]: #intrachromosomal SV
				
				geneChrSubset = genes[sv[0] == genes[:,0]]
				
				#Find all genes overlapped by this SV. 
				startMatches = sv[1] <= geneChrSubset[:,2]
				endMatches = sv[5] >= geneChrSubset[:,1]
				
				matchingGenes = geneChrSubset[startMatches * endMatches]
				
				
				
				for gene in matchingGenes[:,3]:
					svEntry = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[8].sampleName + '_' + sv[8].svType
					gene.SVs[svEntry] = 1
				
			else: #interchromosomal SV
				
				
				
				geneChr1Subset = genes[sv[0] == genes[:,0]]
				#Find all genes overlapped by this SV. 
				startMatches = sv[1] <= geneChr1Subset[:,2]
				endMatches = sv[2] >= geneChr1Subset[:,1]
				
				matchingGenes = geneChr1Subset[startMatches * endMatches]
				for gene in matchingGenes[:,3]:
					svEntry = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[8].sampleName + '_' + sv[8].svType
					gene.SVs[svEntry] = 1
				
				geneChr2Subset = genes[sv[3] == genes[:,0]]
				#Find all genes overlapped by this SV.  
				startMatches = sv[4] <= geneChr2Subset[:,2]
				endMatches = sv[5] >= geneChr2Subset[:,1]
				
				matchingGenes = geneChr2Subset[startMatches * endMatches]
				
				for gene in matchingGenes[:,3]:
					svEntry = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[8].sampleName + '_' + sv[8].svType
					gene.SVs[svEntry] = 1
			
		
		print("Done mapping SVs")
		
		
		#dirty solution here now, exclude the SVs.
		
		filteredSVs = []
		for sv in svData:
			svStr = sv[0] + "_" + str(sv[1]) + "_" + str(sv[2]) + "_" + sv[3] + "_" + str(sv[4]) + "_" + str(sv[5]) + "_" + sv[8].sampleName
			#if svEntry not in excludedSVs:
			#	filteredSVs.append(sv)
			#if svEntry in excludedSVs:
			#	filteredSVs.append(sv)
			
			filteredSVs.append(sv)

		filteredSVs = np.array(filteredSVs, dtype='object')
		
		#filteredSVs = svData
		DerivativeTADMaker(filteredSVs, genes, tadData)
		
		
		#Specific for deletions, this will need to be part of the derivative TAD maker later on
		#if settings.general['gainOfInteractions'] == True:
		#	self.determineGainedInteractions(filteredSVs, tadData)
		
	def mapSNVsToNeighborhood(self, genes, snvDir):
		
		geneMap = dict() #know which gene object to add to by name
		for geneInd in range(0, genes.shape[0]):
			gene = genes[geneInd]
			geneName = gene[3].name
			geneMap[geneName] = geneInd
		
		
		if settings.general['source'] == 'TCGA':
			
			#go through the snv dir, and make a note per gene if it has an SNV in that patient.

			allFiles = [f for f in listdir(snvDir) if isfile(join(snvDir, f))]
			
			for currentFile in allFiles:
				
				if currentFile == "MANIFEST.txt":
					continue
				splitFileName = currentFile.split(".")
				patientID = splitFileName[0]
				splitPatientID = patientID.split("-")
				shortPatientID = 'brca' + splitPatientID[2]
			
				#Load the contents of the file
				with open(snvDir + "/" + currentFile, 'r') as inF:
					lineCount = 0
					for line in inF:
						line = line.strip() #remove newlines
						if lineCount < 1: #only read the line if it is not a header line
							lineCount += 1
							continue
			
						splitLine = line.split("\t")
						geneName = splitLine[0]
						
						if splitLine[8] == 'Silent':
							continue
						
						if geneName not in geneMap:
							continue
						geneInd = geneMap[geneName]
						geneObj = genes[geneInd,3]
						geneObj.addSNV(shortPatientID)
		if settings.general['source'] == 'HMF':
			
			#Make a map from ENSG identifiers to our gene names
			geneNameConversionFile = settings.files['geneNameConversionFile']
			
			geneNameConversionMap = dict()
			with open(geneNameConversionFile, 'r') as inF:
				
				for line in inF:
					
					line = line.strip()
					splitLine = line.split("\t")
					ensgId = splitLine[3]
					splitEnsgId = ensgId.split('.') #we only keep everything before the dot
					geneName = splitLine[4]
					geneNameConversionMap[splitEnsgId[0]] = geneName

			#search through the SNVs and link these to genes.
			vcfs = glob.glob(snvDir + '/**/*.somatic.vcf.gz', recursive=True)
		
			variantsList = []
			addedVariants = [] #check if based on pairs no duplicates are added. 
			
			for vcf in vcfs:
				
				
				#get the samplename from the vcf
				sampleName = re.search('.*\/([A-Z\d]+)\.', vcf).group(1)
				
				#open the .gz file
				with gzip.open(vcf, 'rb') as inF:
					
					for line in inF:
						line = line.strip().decode('utf-8')
	
						if re.search('^#', line): #skip header
							continue
						
						#skip the SV if it did not pass.
						splitLine = line.split("\t")
						filterInfo = splitLine[6]
						if filterInfo != 'PASS':
							continue
				
						#print(line)
						
						#Check if this SNV has any affiliation with a gene. This means that in the info field, a gene is mentioned somewhere. That is, there is an ENSG identifier.
						infoField = splitLine[7]
						
						geneSearch = re.search('(ENSG\d+)', infoField)
						if geneSearch:
							geneMatch = re.search('(ENSG\d+)', infoField).group(1)
							#skip genes for which we do not know the name
							if geneMatch not in geneNameConversionMap:
								continue
							geneName = geneNameConversionMap[geneMatch]
							
							if geneName not in geneMap:
								continue
							
							geneInd = geneMap[geneName]
							geneObj = genes[geneInd,3]
							geneObj.addSNV(sampleName)
	
	
	def mapCNVsToNeighborhood(self, genes, cnvDir):
		
		geneMap = dict() #know which gene object to add to by name
		for geneInd in range(0, genes.shape[0]):
			gene = genes[geneInd]
			geneName = gene[3].name
			geneMap[geneName] = geneInd
		
		if settings.general['source'] == 'HMF':
			
			tsvs = glob.glob(cnvDir + '/**/*.gene.tsv', recursive=True)
			
			for tsv in tsvs:
				
				#get the samplename from the vcf
				sampleName = re.search('.*\/([A-Z\d]+)\.', tsv).group(1)
				
				#open the .gz file
				with open(tsv, 'r') as inF:
					
					lineCount = 0
					for line in inF:
		
						if lineCount < 1: #skip header
							lineCount += 1
							continue
			
						splitLine = line.split("\t")
						
						gene = splitLine[3]
						
						if float(splitLine[5]) > 1.7 and float(splitLine[5]) < 2.3: #these are not CNVs
							continue
						
						if gene not in geneMap:
							continue
						
						geneInd = geneMap[gene]
						geneObj = genes[geneInd,3]
						geneObj.addCNV(sampleName)