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
import sys

from tad import TAD
from sv import SV
from gene import Gene
from derivativeTADMaker import DerivativeTADMaker
from genomicShuffler import GenomicShuffler
from inputParser import InputParser

import settings
from six.moves import range

class NeighborhoodDefiner:
	"""
		Class responsible for defining the neighborhood, 'regulator set', of genes.
		
	"""
	
	def __init__(self, genes, svData):
		"""
			Initialize the neighborhood defining. This involves gathering all required data types, mapping these to TADs/genes, and associating the effects of SVs to genes. 

			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
		"""
		
		#1. Get TADs from the TAD file, and then map TADs to genes (left/right TAD).
		tadData = []

		tadFile = settings.files['tadFile']


		print("Getting TADs")
		tadData = InputParser().getTADsFromFile(tadFile)

		print("original number of svs:", svData.shape)

		if settings.general['shuffleTads'] == True:
			#Shuffle the TADs. Assign random genomic positions to the TADs, but keep the same length.
			genomicShuffler = GenomicShuffler()
			tadData = genomicShuffler.shuffleTADs(tadData)

		print("mapping TADs to genes")
		self.mapTADsToGenes(genes[:,3], tadData)

		
		#2. Get eQTLs from the eQTL file, and map eQTLs to TADs. 			
		eQTLFile = settings.files['eQTLFile']
		print("getting eQTLs")
		eQTLData = InputParser().getEQTLsFromFile(eQTLFile, genes[:,3], self)
		#map the regulatory elements to the TADs so that we can later on when looking at disrupted TADs easily find which elements are affected.
		tadData = self.mapElementsToTads(eQTLData, tadData)
		
		#map the genes to TADs. These are all the gene objects that we can then access when looking at disrupted TADs. 
		tadData = self.mapGenesToTads(genes, tadData) 
		
		#3. Get enhancers

		print("getting enhancers")
		enhancerData = InputParser().getEnhancersFromFile(settings.files['enhancerFile'], genes[:,3], self)
		#Add the enhancers to TADs & genes as well	
		tadData = self.mapElementsToTads(enhancerData, tadData)	
		
		#4. Get promoters
		
		print("getting promoters")
		promoterData = InputParser().getPromotersFromFile(settings.files['promoterFile'], genes[:,3], self)
			
		#Add the promoters to the TADs
		tadData = self.mapElementsToTads(promoterData, tadData)
		
		#5. Get CpG islands
		print("Getting cpg islands")
		cpgData = InputParser().getCpgIslandsFromFile(settings.files['cpgFile'])
		
		#Add the CpG sites to the TADs
		tadData = self.mapElementsToTads(cpgData, tadData)
		
		#6. Get Transcription factors
		print("Getting transcription factors")

		tfData = InputParser().getTranscriptionFactorsFromFile(settings.files['tfFile'])
	
		#Add the CpG sites to the TADs
		tadData = self.mapElementsToTads(tfData, tadData)
		
		
		#7. Get Hi-C data
		print("Getting Hi-C data")
		hicData = InputParser().getHiCInteractionsFromFile(settings.files['hicFile'])
			
		#Map the interactions to TADs as elements
		tadData = self.mapInteractionsToTads(hicData, tadData)
		
		#8. Get histone marks
		
		print("Getting histone marks")
		files = [settings.files['h3k9me3'], settings.files['h3k4me3'], settings.files['h3k27ac'], settings.files['h3k27me3'],
					settings.files['h3k4me1'], settings.files['h3k36me3']]
		types = ['h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
		for histoneFileInd in range(0, len(files)):
			histoneData = InputParser().getHistoneMarksFromFile(files[histoneFileInd], types[histoneFileInd])
			
			#map the histone marks to the TADs
			tadData = self.mapElementsToTads(histoneData, tadData)
		
		#9. Get DNAse I hypersensitivty sites
		print("Getting DNAse I hypersensitivity sites")
			
		dnaseIData = InputParser().getDNAseIFromFile(settings.files['dnaseIFile'])
			
		tadData = self.mapElementsToTads(dnaseIData, tadData)
		
		#10. get chromHMM states
		print("Getting chromHMM states")
		chromHmmData = InputParser().getChromHmmFromFile(settings.files['chromHmmFile'])
			
		tadData = self.mapElementsToTads(chromHmmData, tadData)
		
		#11. get RNAPolII peaks
		print("Getting rnaPol binding sites")
		rnaPolData = InputParser().getRnaPolFromFile(settings.files['rnaPolFile'])
			
		tadData = self.mapElementsToTads(rnaPolData, tadData)
		
		#12. get super enhancers
		print("Getting super enhancers")
		superEnhancerData = InputParser().getSuperEnhancersFromFile(settings.files['superEnhancerFile'])
			
		tadData = self.mapElementsToTads(superEnhancerData, tadData)
		
		#13. get CTCF sites
		print("Getting ctcf sites")
		ctcfData = InputParser().getCTCFSitesFromFile(settings.files['ctcfFile'])
			
		tadData = self.mapElementsToTads(ctcfData, tadData)
		tadData = self.mapCTCFStrengthToTads(ctcfData, tadData)
		
		
		#3. Determine the effect of the SVs on the neighborhood/regulator set
		print("Mapping SVs to the neighborhood")
		self.mapSVsToNeighborhood(genes, svData, tadData)
			
	def mapCTCFStrengthToTads(self, ctcfData, tadData):
		"""
			Function to determine the strength of a TAD boundary. Is not used anymore in later predictions.

			There are 2 ways to look at TAD boundary strength:
			1. The maximum strength of the CTCF sites within a window around the TAD boundary
			2. The number of CTCF sites in a window around the TAD boundary

			ctcfData (numpy array): CTCF data as parsed by InputParser
			tadData (numpt array); TAD data as parsed by InputParser

			return: tadData (np array), input TAD data, but then with the CTCF strengths mapped to them
		"""
		
		#for each tad, check which ctcf site is overlapping the boundaries
		for tad in tadData:
			
			ctcfChrSubset = ctcfData[ctcfData[:,0] == tad[0]]

			# annotate by number of CTCF sites around boundary
			#100kb according to PCAWG 2020 TAD nature paper should be enough
			window = 100000

			startMatches = ctcfChrSubset[(ctcfChrSubset[:,2] >= tad[1]-window) * (ctcfChrSubset[:,1] <= tad[1]+window)]
			endMatches = ctcfChrSubset[(ctcfChrSubset[:,1] >= tad[2]-window) * (ctcfChrSubset[:,2] <= tad[2]+window)]

			if len(startMatches) > 0:
				tad[3].startStrength = startMatches.shape[0]
			if len(endMatches) > 0:
				tad[3].endStrength = endMatches.shape[0]

			#also annotate the CTCF site with the strongest signal in this region
			maxStart = 0
			if len(startMatches) > 0:

				for match in startMatches:
					if match[5] > maxStart:
						maxStart = match[5]

			tad[3].startStrengthSignal = maxStart

			maxEnd = 0
			if len(endMatches) > 0:
				for match in endMatches:

					if match[5] > maxEnd:
						maxEnd = match[5]

			tad[3].endStrengthSignal = maxEnd
		
		return tadData
		
	
	def mapTADsToGenes(self, genes, tadData):
		"""
			Adds the left and right TAD to each gene object.
			
			DEPRECATED:
			- This function may be deprecated, I don't really use it anymore in the code but leave it here in case it is useful at a later point. 
			
		"""
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

			element (np array): element to be added to the gene
			geneDict (dictionary): Dictionary with gene names and gene objects in the values.
			geneSymbol (string): Name of the gene that we want to map the element to

		"""

		geneDict[geneSymbol].addElement(element)

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
			Determine which element interactions take place within which TAD.
			element interactions will be mapped to TAD objects.

			elementData: (numpy array) array with elements. chr, start, end, geneObject, strength

			return:
			tadData (numpy array): array with TAD information (from InputParser), updated with elements mapped to the TADs.
			
		"""
		
		#Go through each TAD and find which elements are within the TAD.
		#List these elements as interaction objects in the tAD object.
		previousChromosome = 0
		for tad in tadData:
			if tad[0] != previousChromosome:
				previousChromosome = tad[0]
				elementChrSubset = elementData[np.where(elementData[:,0] == tad[0])] #Get all elements that are on this chromosome. All elements are CIS, so intrachromosomal is fine for now

			#Find the element where the start and end of the element are within the start and end of the TAD. 
			startMatches = elementChrSubset[:,1] >= tad[1]
			endMatches = elementChrSubset[:,2] <= tad[2]
			
			allMatches = startMatches * endMatches
			matchingElements = elementChrSubset[allMatches,:]

			#Add these elements to the TADs if any.
			if matchingElements.shape[0] < 1:
				
				continue
			else:

				tad[3].addElements(matchingElements) #Add the elements objects to the TAD objects.

		return tadData		
	
	def mapInteractionsToTads(self, interactionsByTad, tadData):
		"""
			Function to specifically map Hi-C interactions to the TADs. Hi-C interactions have 2 sides, each of which is considered a separate regulatory element that can be
			gained or lost.

			interactionsByTad: for each interaction side, we have pre-mapped it to the TAD that it takes place in (see preprocessing.sh). Use that here for a quick mapping.
			tadData (numpy array): array with TAD information (from InputParser)

			return:
			array with TAD information (from InputParser), updated with elements mapped to the TADs.
			
		"""
		
		print("mapping to tads")
		uniqueElements = dict()
		for tadStr in interactionsByTad:

			#1. Get the right TAD object that this interaction bin is in. 
			splitTadStr = tadStr.split("_")

			tadChrMatch = tadData[:,0] == splitTadStr[0]
			tadStartMatch = tadData[:,1] == int(splitTadStr[1])
			tadEndMatch = tadData[:,2] == int(splitTadStr[2])

			matchingTad = tadData[tadChrMatch * tadStartMatch * tadEndMatch]

			if len(matchingTad) < 1: #Sometimes the interaction are outside TADs
				continue

			matchingTad = matchingTad[0] #assume that there is 1 TAD the bin falls into. We have 1 start coordinate for the interaction, so we assume that it is that one that is within the right TAD.
			#Assign the interactions to this TAD.
			interactionLines = []
			binSize = 5000 #should be a setting to work with other Hi-C data than we currently use.

			for lineInd in range(0, len(interactionsByTad[tadStr])-1):
				line = interactionsByTad[tadStr][lineInd]
				element = [splitTadStr[0], int(line), int(line)+binSize, "hic", None, None] #add the bin size to indicate the start/end of the interaction to use for the element. 
				interactionLines.append(element)

			#Assign the interactions to this TAD.
			matchingTad[3].addElements(interactionLines)
		
		return tadData 
		
	def mapSVsToNeighborhood(self, genes, svData, tadData):
		"""
			Take as input gene objects with all regulatory elements mapped to them that could potentially be disrupted, and search through the SVs to find which SVs overlap the genes,
			and what the changes to the regulatory set are as a result of the SVs.

			First looks at 'coding' effects; what are the SVs that directly overlap genes? This is output to a file, but we don't use it anymore anywhere.
			Then, it calls the DerivativeTADMaker, which looks at each SV and the TADs that these disrupt, and then calculate gains and losses for each gene.
		
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
		
		#filteredSVs = svData
		DerivativeTADMaker(svData, genes, tadData)
