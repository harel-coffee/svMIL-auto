import numpy as np
import json
import pickle as pkl
import re

from tad import TAD
from sv import SV
from gene import Gene
from eQTL import EQTL
from interaction import Interaction
from snv import SNV
from derivativeTADMaker import DerivativeTADMaker
from genome import Genome
from genomicShuffler import GenomicShuffler


import settings

class NeighborhoodDefiner:
	"""
		Class responsible for defining the neighborhood of causal genes.
		
		Currently, the neighborhood consists of:
		
		- nearest TADs on the left and right of the gene
		- all eQTLs mapped to the gene (and if these are enhancers or not)
		- Hi-C interactions. These are mapped to the TAD that these are present in. We use this to determine which interactions are potentially gained. 
		- SVs (and SNVs) overlapping either the gene directly, or other elements in the neighborhood
	
	"""

	
	def __init__(self, genes, svData, snvData, mode, genome):
		
		#1. Map genes to TADs
		
		tadData = []
		if settings.general['tads'] == True or settings.general['gainOfInteractions'] == True: #Gain of interactions is dependent on TADs
			
			#Make these pats a setting!
			tadFile = "../../data/tads/HMEC_Lieberman-raw_TADs.bed" #These files will need to go to the settings!
			#tadFile = "../../data/tads/prostate_tads.txt"
			#tadFile = "../../data/tads/ovary_hg38_liftover.bed"
			#tadFile = "../../data/tads/tadsNoCellType.bed"
			print "Getting TADs"
			tadData = self.getTADsFromFile(tadFile)
			
			
			# maxTadLength = 0
			# largestTad = None
			# minTadLength = float("inf")
			# smallestTad = None
			# for tad in tadData:
			# 	tadLength = tad[2] - tad[1]
			# 	if tadLength > maxTadLength:
			# 		maxTadLength = tadLength
			# 		largestTad = tad
			# 	if tadLength < minTadLength:
			# 		minTadLength = tadLength
			# 		smallestTad = tad
			# 
			# print largestTad, largestTad[2] - largestTad[1]
			# print smallestTad, smallestTad[2] - smallestTad[1]
			# exit()
			
			if settings.general['shuffleTads'] == True:
				genomicShuffler = GenomicShuffler()
				#Shuffle the TADs. Assign random genomic positions to the TADs, but keep the same length. 
				tadData = genomicShuffler.shuffleTADs(tadData)

				
			
			print "mapping TADs to genes"
			self.mapTADsToGenes(genes[:,3], tadData) #only pass the gene objects will suffice
			
		#2. Map genes to eQTLs
		
				
		eQTLData = [] #Keep empty by default
		if settings.general['eQTLs'] == True or settings.general['gainOfInteractions'] == True: #Gain of eQTL interactions depends on eQTLs. 
			#Save the processed data, only needs to be done once
			
			# import os.path
			# 
			# if os.path.exists('eQTLData.pkl'):
			# 	print "loading eQTLs from pkl"
			# 	#Load the eqtls
			# 	with open('eQTLData.pkl', 'rb') as h:
			# 		eQTLData = pkl.load(h)
			# else:
			print "re-creating eQTLs"
			#eQTLFile = "../../data/eQTLs/eQTLsFilteredForCausalGenes.txt" #These files will need to go to the settings!
			eQTLFile = "../../data/eQTLs/breast_eQTLs.txt" #These files will need to go to the settings!
			#eQTLFile = "../../data/eQTLs/prostate_eQTLs.txt"
			#eQTLFile = "../../data/eQTLs/ovarian_eQTLs.txt"
			
			print "getting eQTLs"
			eQTLData = self.getEQTLsFromFile(eQTLFile, genes[:,3])
		
			
			with open('eQTLData.pkl', 'wb') as h:
				pkl.dump(eQTLData, h, protocol=pkl.HIGHEST_PROTOCOL)
		
	
		
		#Add reading/parsing of 3D genome information
		#Disable the 3D information for now
		
		# if settings.general['interactionChains'] == True: 
		# 	interactionsFile = "../../data/HiC/intrachromosomal_geneNonGeneInteractions.csv" #should be setting!!!!
		# 	print "Getting interactions"
		# 	interactionData = self.getInteractionsFromFile(interactionsFile)
		# 	print "mapping interactions to genes"
		# 	self.mapInteractionsToGenes(genes[:,3], interactionData)
			
		#Read/parse Hi-C interactions (intrachromosomal) to later compute gains of interactions
		#We ignore the Hi-C information in this for now, we now only use eQTLs for gains of interactions
		# if settings.general['gainOfInteractions'] == True:
		# 	interactionsFile = settings.files['hiCInteractionsFile']
		# 	print "Reading all Hi-C interactions"
		# 	#First get all Hi-C interactions
		# 	interactions, regions = self.getHiCInteractionsFromFile(interactionsFile)
		# 	#Map the Hi-C interactions to the respective TADs
		# 	tadData = self.mapInteractionsToTads(interactions, regions, tadData)

		#lncRNA data
		if settings.general['lncRNA'] == True:
			lncRNAData = self.getLncRNAsFromFile(settings.files['lncRNAFile'])
			eQTLData = lncRNAData 

		#For now, we use eQTLs for gains of interactions rather than Hi-C.
		if settings.general['gainOfInteractions'] == True:
			tadData = self.mapEQTLInteractionsToTads(eQTLData, tadData)
			tadData = self.mapGenesToTads(genes, tadData) 
		
		
		
		
		#First define the Genome object with all data that we collected to far
		print "Defining genomic bins"
		genome.defineBins(genes, tadData, eQTLData)
		
		
		#3. Map SVs to all neighborhood elements
		if mode == "SV":
			print "Mapping SVs to the neighborhood"
			self.mapSVsToNeighborhood(genes, svData, tadData, genome)
		if mode == "SNV":
			print "Mapping SNVs to the neighborhood"
			self.mapSNVsToNeighborhood(genes, snvData, eQTLData)
		if mode == "SV+SNV": #in this case map both
			print "Mapping SVs to the neighborhood"
			self.mapSVsToNeighborhood(genes, svData, tadData, genome)
			
			print "Mapping SNVs to the neighborhood"
			self.mapSNVsToNeighborhood(genes, snvData, eQTLData)
		
		
	def getTADsFromFile(self, tadFile):
		"""
			Read the TADs into NumPy format. I don't read the TADs into objects immediately, because the NumPy matrices work very well for
			quick overlapping. I add an object reference to the matrix so that we can later add the right TAD object to the genes. 
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
				
				
				#Convert the numbers to integers for quicker comparison. 0 is chromosome, 1 is start, 2 is end. Celltypes are not used for now. 
				# splitLine[1] = int(splitLine[1])
				# splitLine[2] = int(splitLine[2])
				
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
		
		
	def getEQTLsFromFile(self, eQTLFile, genes):
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
				eQTLObject = EQTL(chrName, int(splitLine[1]), int(splitLine[2])) #chr, start, end
				
				#The mapping information is in the file, so we can already do it here
				self.mapEQTLsToGenes(eQTLObject, geneDict, splitLine[3])
						
				
				eQTLs.append([chrName, int(splitLine[1]), int(splitLine[2]), eQTLObject]) #Keep the eQTL information raw as well for quick overlapping. 
		
		
		return np.array(eQTLs, dtype='object')
	
	def getLncRNAsFromFile(self, lncRNAFile):
		
		lncRNAs = []
		with open(lncRNAFile, 'rb') as f:
			
			for line in f:
				line = line.strip()
				splitLine = line.split("\t")
				
				#Quick and dirty
				lncRNAObject = EQTL(splitLine[0], int(splitLine[1]), int(splitLine[2]))
				
				lncRNAs.append([splitLine[0], int(splitLine[1]), int(splitLine[2]), lncRNAObject])					
		
		return np.array(lncRNAs, dtype="object")
		

	def getInteractionsFromFile(self, interactionsFile):
		"""
			!!! DEPRECATED, Hi-C BASED HEAT DIFFUSION IS NOT WORKING
		
			Read the HiC interactions between genes and interaction pairs into numPy array format. Each numpy array will be in the format of:
			0: chromosome of non-coding interaction
			1: start of non-coding interaction
			2: end of non-coding interaction
			3: name of gene that is interacted with. If there are multiple genes, this is split across multiple lines
			4: bin in which the gene resides in the format of chr_binStartPosition
			
			Notes:
			There is an issue with the csv export in Neo4J and now there are double quotation marks which I cannot fix. Thus the string is not real json. 
		"""
		
		#Read the gene list data into a list
		interactionData = []
		with open(interactionsFile, "r") as f:
			lineCount = 0
			for line in f:
				if lineCount < 1: #skip header
					lineCount += 1
					continue
				line = line.strip()
				splitLine = line.split('}","{') #split on every column which is encoded strangely now
				
				encodedNonCodingRegion = splitLine[0].split('"')[7] #this is where the region ID will always be
				nonCodingRegionChr = encodedNonCodingRegion.split("_")[0]
				nonCodingRegionStart = int(encodedNonCodingRegion.split("_")[1])
				nonCodingRegionEnd = nonCodingRegionStart + 5000 #bin of 5kb, should be a setting
				
				encodedGeneData = splitLine[1].split('"')
				encodedGeneRegion = encodedGeneData[14]
				encodedGeneRegionChr = encodedGeneRegion.split("_")[0]
				encodedGeneRegionStart = int(encodedGeneRegion.split("_")[1])
				encodedGeneRegionEnd = encodedGeneRegionStart + 5000
				encodedGeneNames = encodedGeneData[6]
				
				splitGeneNames = encodedGeneNames.split(";")
				for geneName in splitGeneNames: #Add the region separately for each gene in the coding region. 
					
					interactionObj = Interaction(nonCodingRegionChr, nonCodingRegionStart, nonCodingRegionEnd, encodedGeneRegionChr, encodedGeneRegionStart, encodedGeneRegionEnd, geneName)
					#For now I will not make objects because working with the numpy arrays is much much faster
					interactionData.append([nonCodingRegionChr, nonCodingRegionStart, nonCodingRegionEnd, geneName, encodedGeneRegion, interactionObj])
					
	
		#Also convert the dataset to numpy
		interactionData = np.array(interactionData, dtype='object')

		return interactionData
	
	def getHiCInteractionsFromFile(self, interactionsFile):
		"""
			Read all Hi-C interactions from the interactions file
			
			- Column 1 is the starting region of the interaction
			- Column 2 is the ending region of the interaction
			
			
			
		"""
		seenRegions = dict() #use a dictionary to quickly determine if we have added this region before to keep the regions unique
		regions = []
		interactions = dict() #for now I won't make objects for interactions, do we really need them? 
		with open(interactionsFile, 'r') as inF:
			
			lineCount = 0
			for line in inF:
				line = line.strip()
				
				if lineCount < 1: #skip header
					lineCount += 1
					continue
				
				splitLine = line.split(",") #csv format

				interactionStart = splitLine[0]
				interactionEnd = splitLine[1]
				
				#Split the regions into the chromosome and region/bin
				splitInteractionStart = interactionStart.split("_")
				splitInteractionEnd = interactionEnd.split("_")
				
				chr1 = splitInteractionStart[0]
				start1 = int(splitInteractionStart[1])
				end1 = start1 + int(settings.interactions['binSize'])
				
				chr2 = splitInteractionEnd[0]
				start2 = int(splitInteractionEnd[1])
				end2 = start2 + int(settings.interactions['binSize'])
				
				if interactionStart not in seenRegions:
					regions.append([chr1, start1, end1, interactionStart])
					seenRegions[interactionStart] = len(regions) #keep the index at which the region is
				if interactionEnd not in seenRegions:
					regions.append([chr2, start2, end2, interactionEnd])
					seenRegions[interactionEnd] = len(regions)
				
				if interactionStart not in interactions:
					interactions[interactionStart] = []
				if interactionEnd not in interactions:
					interactions[interactionEnd] = [] #Some interactions are only in the end region
				
				interactions[interactionStart].append(interactionEnd)
				interactions[interactionEnd].append(interactionStart)
				
		
		regions = np.array(regions, dtype="object")
		#interactions = np.array(interactions, dtype="object")

		return interactions, regions
	
	def mapTADsToGenes(self, genes, tadData):
		"""
			Adds the left and right TAD to each gene object.
		"""
		print tadData
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
			# 	
			# if gene.name == "ARID1A":
			# 	print nearestRightTAD.start, nearestRightTAD.end
			# 	print nearestLeftTAD.start, nearestLeftTAD.end
			# 	exit()

		
	def mapEQTLsToGenes(self, eQTL, geneDict, geneSymbol):
		"""
			Map the right gene object to the eQTL object. 
		
		"""
		
		geneDict[geneSymbol].addEQTL(eQTL)
		eQTL.addGene(geneDict[geneSymbol])
		# 
		# for gene in genes:
		# 	
		# 	if gene.name == geneSymbol:
		# 		
		# 		gene.addEQTL(eQTL)
		# 		eQTL.addGene(gene) #Also update the gene object in the eQTL object so that we can find the gene of the eQTL back later by eQTL
	
	def mapInteractionsToGenes(self, genes, interactionData):
		"""
			Take Hi-C interactions as input (see getInteractionsFromFile for the format) and link these to the causal genes. 
		
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
				else:

					nearestGeneInd = np.argmin(startDistances)
					nearestGene = geneChrSubset[nearestGeneInd,:]
					#tad[3].addGene(nearestGene[3]) #Check the performance with genes inside TADs only
					
				
				#Repeat but then for the other TAD end. 
				endDistances = geneChrSubset[:,1] - tad[2]
				#Exclude all negative elements
				negativeElements = endDistances < 0

				endDistances[negativeElements] = float("inf") #skip the ones that are on the wrong side of the TAD boundary, so set them to inf and they will never be selected. 
				if endDistances[negativeElements].shape[0] == endDistances.shape[0]: #but if everything is inf, we skip this TAD.
					continue
				else:
					nearestGeneInd = np.argmin(endDistances)
					nearestGene = geneChrSubset[nearestGeneInd,:]
					#tad[3].addGene(nearestGene[3])
				
			else:
				tad[3].setGenes(matchingGenes[:,3]) #Add the eQTL objects to the TAD objects. 
	
		return tadData		
	
	def mapEQTLInteractionsToTads(self, eQTLData, tadData):
		"""
			Determine which eQTL interactions take place within which TAD.
			eQTL interactions will be mapped to TAD objects.
			
		"""
		
		#Go through each TAD and find which eQTLs are within the TAD.
		#List these eQTLs as interaction objects in the tAD object.
		
		
		previousChromosome = 0
		for tad in tadData:
			
		
			if tad[0] != previousChromosome:
				previousChromosome = tad[0]
				eQTLChrSubset = eQTLData[np.where(eQTLData[:,0] == tad[0])] #Get all eQTLs that are on this chromosome. All eQTLs are CIS, so intrachromosomal is fine for now



			#Find the eQTLs where the start and end of the eQTL are within the start and end of the TAD. 
			startMatches = eQTLChrSubset[:,1] >= tad[1]
			endMatches = eQTLChrSubset[:,2] <= tad[2]
			
			allMatches = startMatches * endMatches
			matchingEQTLs = eQTLChrSubset[allMatches,:]
			
			#Add these eQTLs to the TADs if any.
			if matchingEQTLs.shape[0] < 1:
				
				continue
			else:
				
				tad[3].setEQTLInteractions(matchingEQTLs[:,3]) #Add the eQTL objects to the TAD objects.
			
		return tadData		
	
	def mapInteractionsToTads(self, interactions, regions, tadData):
		"""
			Determine which interactions take place within which TAD.
			Interactions will be mapped to TAD objects.
			
			- For now, I remove all interactions that take place with regions outside of the TAD. For the purpose of this script, we only look at interactions that are gained as a
			result of disrupted TAD boundaries, so regions that are already outside of TAD boundaries are questionable. Are these errors? Are the TAD boundary definitions not correct?
			
			
			Returns the tadData, where the objects now have interactions mapped to them (but these are objects, so by reference this should not be necessary)
			
		"""
		print "Mapping interactions to TADs"
		previousChromosome = 0
		for tad in tadData:
			
			#Find all interactions taking place within this TAD.
			
			#The interactions are currently intrachromosomal, but may later be interchromosomal
			#Thus for now we can just match on 1 chromosome
			
			#Get all regions that are within the TAD
			#Find all corresponding interactions between the regions within the TAD (make sure that these do not go outside of the boundary, this may indicate data errors)
			
			if tad[0] != previousChromosome:
				previousChromosome = tad[0]
				regionChrSubset = regions[np.where(regions[:,0] == tad[0])] #First get the subset of all regions on the current chromosome, TADs are sorted

			
			#Find which regions are within this TAD
			#All regions must have their start position larger than the TAD start and the end smaller than the TAD end
			startMatches = regionChrSubset[:,1] >= tad[1]
			endMatches = regionChrSubset[:,2] <= tad[2]
			
			allMatches = startMatches * endMatches
			matchingRegions = regionChrSubset[allMatches,:]
			
			#Get the interactions of these matching regions and ensure that these are within the TAD boundaries
			
			#Find all the interactions mapping to this region.
			 
			matchingRegionsInteractions = []
			for region in matchingRegions:
				regionInteractions = interactions[region[3]] #the interacting regions
				for interaction in regionInteractions: #make sure that the interacting regions are actually within the TAD
					
					splitInteraction = interaction.split("_")
					interactionChrom = splitInteraction[0]
					splitInteractionStartRegion = int(splitInteraction[1])
					splitInteractionEndRegion = splitInteractionStartRegion + settings.interactions['binSize']
					
					if splitInteractionStartRegion < tad[2]: #remove the interacting regions where the start of the interaction is outside of the TAD
						#Make the final format for the interactions and then add them to the TAD
						markedUpInteraction = [region[0], region[1], region[2], interactionChrom, splitInteractionStartRegion, splitInteractionEndRegion]
						matchingRegionsInteractions.append(markedUpInteraction)
						
						
			#Add the interactions to this TAD			
			tad[3].setInteractions(matchingRegionsInteractions)
			
		
		
		return tadData 
		
	def determineGainedInteractions(self, svData, tadData):
		"""
			- Function name??
			- Split better into respective parts
		"""
		
		#Loop through all SVs and see which ones have overlap with any TAD.
		#The SVs are sorted per cancer type, so taking chromosome subsets may not gain as much speed

		for sv in svData:
			
			#Only focus on deletions
			typeMatch = re.search("del", sv[8].svType, re.IGNORECASE)
			if typeMatch is None:
				continue
			
			#print "sv: ", sv
			
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

				#Find all SVs where the start of the SV is before the TAD end, and the end of SV after the TAD start. These are boundary overlapping SVs.
				#Then make sure to remove all SVs that are entirely within a TAD, we cannot say anything about the effect of these
				
				# startMatches = svStart <= tadChromosomeSubset[:,2]
				# endMatches = svEnd >= tadChromosomeSubset[:,1]
				# 
				# allMatches = startMatches * endMatches
				# 
				# 
				# #Find the SVs that are entirely within TADs
				# withinStartMatches = svStart > tadChromosomeSubset[:,1]
				# withinEndMatches = svEnd < tadChromosomeSubset[:,2]
				# 
				# withinTadMatches = withinStartMatches * withinEndMatches
				# #print withinTadMatches.shape
				# #Now everything that is true in the within vector should be FALSE in the allMatches vector.
				# #We can make the within vector negative, and then use *. Everything that is false in the previous vector will remain false, but what is true may become false.
				# 
				# filteredMatches  = -withinTadMatches * allMatches
				# #print filteredMatches.shape
				# 
				# overlappingTads = tadChromosomeSubset[allMatches]
				
				if overlappingTads.shape[0] < 1:
					continue
				
			if sv[0] != sv[3]: #interchromosomal
				continue #skip these translocations for now
				leftTadChr = sv[0]
				rightTadChr = sv[3]
				svStartLeft = sv[1] #s1
				svEndLeft = sv[2] #e1
				svStartRight = sv[4] #s2
				svEndRight = sv[5] #e2

				#Now search for the matches on the left and right TAD by these coordinates
				#First get the right chromosome subset
				tadChromosomeSubset = tadData[np.where(tadData[:,0] == leftTadChr)] #Limit to intrachromosomal for now for testing purposes

				#For the interchromosomal case, we match chromosome 1 by all SVs disrupting the TAD boundary with the s1 and e1 coordinates
				
				startMatches = svStartLeft <= tadChromosomeSubset[:,2]
				endMatches = svEndLeft >= tadChromosomeSubset[:,1]
				
				allMatchesLeft = startMatches * endMatches
				
				overlappingTadsLeft = tadChromosomeSubset[allMatchesLeft] #There should be no need to remove intra-TAD SVs in the interchromosomal case. 
				if overlappingTadsLeft.shape[0] < 1:
					continue
				
				#Repeat for chromosome 2
				tadChromosomeSubset = tadData[np.where(tadData[:,0] == rightTadChr)] #Limit to intrachromosomal for now for testing purposes

				#For the interchromosomal case, we match chromosome 1 by all SVs disrupting the TAD boundary with the s1 and e1 coordinates
				
				startMatches = svStartRight <= tadChromosomeSubset[:,2]
				endMatches = svEndRight >= tadChromosomeSubset[:,1]
				
				allMatchesRight = startMatches * endMatches
				
				overlappingTadsRight = tadChromosomeSubset[allMatchesRight] #There should be no need to remove intra-TAD SVs in the interchromosomal case. 
				if overlappingTadsRight.shape[0] < 1:
					continue
				
				#Combine the two
				overlappingTads = np.concatenate((overlappingTadsLeft, overlappingTadsRight), axis=0)
			
			
				
			#For the overlapping TADs, we can get the TAD on the right and the TAD on the left.
			#The SV could end on the boundary and then the desired result is not in our match set, so we skip these for now.
			
			#Get the TADs on the far left and far right of the SV.
			#These are the TADs that are still in tact.
			#We get the eQTL interactions that are in these two TADs (have been pre-mapped to the TAD that these take place in)
			#Then we find the genes that are present within the TAD (has also been pre-mapped to the TAD)
			#To each gene in one TAD, we add the number of gained interactions equal to the number of eQTL interactions in the other TAD.
				
			farLeftTad = overlappingTads[0] #This list is sorted
			farRightTad = overlappingTads[overlappingTads.shape[0]-1,:]
			
			
			
	
			#For every gene in the TAD, add the eQTLs of the other TAD as potentially gained interactions.
			for gene in farLeftTad[3].genes:
				
				if gene.name == "PTK6":
					print "TADs:"
					print farLeftTad
					print farRightTad
					print sv
				
				# if gene.name == "SYT14":
				#  	print "left tad: ", farLeftTad
				#  	print "sv: ", sv
				if len(farRightTad[3].eQTLInteractions) > 0:
					gene.setGainedEQTLs(farRightTad[3].eQTLInteractions, sv[7])
					
			
			for gene in farRightTad[3].genes:
				
				if gene.name == "PTK6":
					print "TADs:"
					print farLeftTad
					print farRightTad
					print sv
				# if gene.name == "SYT14":
				# 	print "right tad: ", farRightTad
				# 	print "sv: ", sv
				# 	print "gained eQTLs: ", len(farLeftTad[3].eQTLInteractions)
				# 	print "gaining from tad: ", farLeftTad
				if len(farLeftTad[3].eQTLInteractions) > 0:	
					gene.setGainedEQTLs(farLeftTad[3].eQTLInteractions, sv[7])
				
		
		return 0
		
		
	def mapSVsToNeighborhood(self, genes, svData, tadData, genome):
		"""
			Take as input gene objects with the neighborhood pre-set, and search through the SVs to find which SVs overlap the genes, TADs and eQTLs in the gene neighborhood
		
			The TADs are also parsed as input, because when we compute gained interactions we search through other TADs that are on the other end of the SV breakpoint. This is the TAD from which we can
			expect interactions in our current TAD. 
		
			!!!!!! This function will need to be properly split into multiple functions for the different data types.
			TODO:
			- Split into multiple functions
			- Move the gained interactions outside of the left/right TAD check. See comments below. 
		
		"""
		
		DerivativeTADMaker(svData, genes, tadData, genome)
		
		
		#First map the SVs to TADs to see if we can infer gained interactions
		if settings.general['gainOfInteractions'] == True:
			self.determineGainedInteractions(svData, tadData)
		
		#The code below is now very specifically for all elements, but I will initially remove the TAD part and focus just on eQTLs.
		#eQTLs will only be labelled as lost when these are targeted by a deletion, and also when these are within the nearest TAD boundaries of their gene. 
		
		previousChr = None
		for gene in genes[:,3]:
			
			#We first split the data per chromosome for faster overlapping. Thus, if the chromosome of the gene is 1, then we find all SVs for which either chromosome 1 or chromosome 2 (translcoations)
			#are on chromosome 1. The genes are sorted, so this saves time when doing overlap with large sets of SVs or SNVs. 
			
			#To overlap, there are SVs that either have overlap with elements in the gene neighborhood on chr1 or chr2 (if translocation). Thus, we make 3 different sets
			#for SVs that are not translocations (checking if overlapping with chromosome 1, start 1 and end 2 (intraChrSubset)).
			#If the SV is a translocation, we can match either on chromosome 1 or chromosome 2. Thus we have 2 sets of SVs on chr1, where we overlap with s1 and e1, and on chr2, where we overlap with
			#s2 and e2. These are then all combined into one large set, so that we can immediately do the overlap at once for all elements in the neighborhood on the chromosomes of the SVs. 
			
			if gene.chromosome != previousChr:
				
				matchingSvChr1Ind = svData[:,0] == str(gene.chromosome)
				matchingSvChr2Ind = svData[:,3] == str(gene.chromosome)
				
				#Intra and chr1 and chr2 will overlap if we don't exclude the positions where both chr1 and chr2 are the same. 
		
				notChr1Matches = svData[:,3] != str(gene.chromosome)
				chr1OnlyInd = matchingSvChr1Ind * notChr1Matches
				
				notChr2Matches = svData[:,0] != str(gene.chromosome)
				chr2OnlyInd = matchingSvChr2Ind * notChr2Matches
				
				#intraSubset: chr1 and chr2 both match
				matchingChr1AndChr2Ind = matchingSvChr1Ind * matchingSvChr2Ind
				intraSubset = svData[matchingChr1AndChr2Ind]
				
				#interChr1Subset: only chr1 matches
				interChr1Subset = svData[chr1OnlyInd]
				
				#interChr2Subset: only chr2 matches
				interChr2Subset = svData[chr2OnlyInd]
				
		
				#Now concatenate them into one set, but keep the formatting the same as: chr, start, end
				
				svChr1Subset = np.empty([interChr1Subset.shape[0],11],dtype='object')
				svChr1Subset[:,0] = interChr1Subset[:,0] #For chromosome 1, we use just the first chromosome, s1 and e1.
				svChr1Subset[:,1] = interChr1Subset[:,1] #For chromosome 1, we use just the first chromosome, s1 and e1.
				svChr1Subset[:,2] = interChr1Subset[:,2] #For chromosome 1, we use just the first chromosome, s1 and e1.
				svChr1Subset[:,3] = None #Here fill with None because the SVs need to have the cancer type and sample name in the same place in the array as the SNVs, but the SNVs don't have this info. Also just use None because we won't use the other position anymore.
				svChr1Subset[:,4] = None
				svChr1Subset[:,5] = None
				svChr1Subset[:,6] = interChr1Subset[:,7]
				svChr1Subset[:,7] = interChr1Subset[:,6]
				
				#Make the subset for the chr2 matches
				svChr2Subset = np.empty([interChr2Subset.shape[0],11], dtype='object')
				svChr2Subset[:,0] = interChr2Subset[:,0] #For chromosome 2, we use just the second chromosome, s2 and e2.
				svChr2Subset[:,1] = interChr2Subset[:,4] 
				svChr2Subset[:,2] = interChr2Subset[:,5] 
				svChr2Subset[:,3] = None
				svChr2Subset[:,4] = None
				svChr2Subset[:,5] = None
				svChr2Subset[:,6] = interChr2Subset[:,7] 
				svChr2Subset[:,7] = interChr2Subset[:,6] 
				
				
				#For the intra subset, we need to use s1 and e2.
				svIntraSubset = np.empty([intraSubset.shape[0],11], dtype='object')
				svIntraSubset[:,0] = intraSubset[:,0] #For chromosome 2, we use chromosome 1, s1 and e2.
				svIntraSubset[:,1] = intraSubset[:,1] 
				svIntraSubset[:,2] = intraSubset[:,5] 
				svIntraSubset[:,3] = None
				svIntraSubset[:,4] = None
				svIntraSubset[:,5] = None
				svIntraSubset[:,6] = intraSubset[:,7] 
				svIntraSubset[:,7] = intraSubset[:,6] 
				
				#Now concatenate the arrays
				svSubset = np.concatenate((svChr1Subset, svChr2Subset, svIntraSubset))
			
				previousChr = gene.chromosome
			
			if np.size(svSubset) < 1:
				continue #no need to compute the distance, there are no TADs on this chromosome
			
			#Check which of these SVs overlap with the gene itself
			
			geneStartMatches = gene.start <= svSubset[:,2]
			geneEndMatches = gene.end >= svSubset[:,1]
			
			geneMatches = geneStartMatches * geneEndMatches #both should be true, use AND for concatenating
		
			svsOverlappingGenes = svSubset[geneMatches]
			
			
			#Get the SV objects and link them to the gene
			gene.setSVs(svsOverlappingGenes)
			
			#Check which SVs overlap with the eQTLs
			
			#Repeat for eQTLs. Is the gene on the same chromosome as the eQTL? Then use the above chromosome subset.
			
			#In the input parser we already removed all non-deletions. This should be a check around here later on. 
			
			geneEQTLs = gene.eQTLs
			
			for eQTL in geneEQTLs: #only if the gene has eQTLs
				
				startMatches = eQTL.start <= svSubset[:,2]
				endMatches = eQTL.end >= svSubset[:,1]
				
				allMatches = startMatches * endMatches

				svsOverlappingEQTL = svSubset[allMatches]
				
				#Filter the matches further on if these are within the gene's TAD boundary. Position of eQTL > left TAD, < right TAD
				
				if gene.leftTAD == None or gene.rightTAD == None or len(svsOverlappingEQTL) == 0: #if there are no TAD boundaries, we cannot check losses
					continue
				
				#Assume that start and end are so small for eQTLs that usig only the start is enough
				startMatches = eQTL.start >= gene.leftTAD.end
				endMatches = eQTL.start <= gene.rightTAD.start

				#if this eQTL is within the TAD boundaries of the gene, count the samples in which it is disrupted by an SV. 
				match = startMatches * endMatches
				
				
				if match == True:
					
					#Make sure that the SVs that overlap the gene as well are not counted
					for sv in svsOverlappingEQTL:
						for sv2 in svsOverlappingGenes:
							#if sv not in svsOverlappingGenes:
							if sv[0] != sv2[0] and sv[1] != sv2[1] and sv[5] != sv2[5]:
								gene.addLostEQTL(eQTL, sv[6])
					
					# 
					# samples = np.unique(svsOverlappingEQTL[:,6])
					# 
					# #Add the eQTL to the list of lost eQTLs for this gene
					# ##Temporarily turned this off for testing without deletions
					# for sample in samples:
					# 	# if gene.name == "EPS15":
					# 	# 	print "gene: ", gene.name, " loses eQTL in sample ", sample
					# 	# 	print gene.leftTAD.start, gene.leftTAD.end
					# 	# 	print gene.rightTAD.start, gene.rightTAD.end
					# 	# 	print "eQTL: ", eQTL.start
					# 	# exit()	
					# 	gene.addLostEQTL(eQTL, sample)
				
				#If there is anything in svsWithinTAD, then we can count ths eQTL as a lost eQTL
				
				eQTL.setSVs(svsOverlappingEQTL)
				
				
			
			#Check which SVs overlap with the interactions.
			# for interaction in gene.interactions:
			# 	
			# 	startMatches = interaction.start1 <= svSubset[:,2]
			# 	endMatches = interaction.end1 >= svSubset[:,1]
			# 	
			# 	allMatches = startMatches * endMatches
			# 	
			# 	svsOverlappingInteraction = svSubset[allMatches]
			# 	
			# 	interaction.setSVs(svsOverlappingInteraction)
		
		
		
	def mapSNVsToNeighborhood(self, genes, snvData, eQTLData):
		"""
			Same as the function for mapping SVs to the neighborhood, but then for SNVs.
			!!!!!  Will also need to be properly split into multiple functions, many pieces of code can be re-used. 
		
			TODO:
			- update the documentation in this code
			- Split into functions, re-use code for the SVs
		
		"""
		import time
		import math
		
		
		# startTime = time.time()
		# previousChr = None
		# for gene in genes[:,3]:
		# 	
		# 	
		# 	#Make a subset of SNVs on the right chromosome
		# 	
		# 	if gene.chromosome != previousChr:
		# 		
		# 		endTime = time.time()
		# 		
		# 		print "new chromosome: ", gene.chromosome
		# 		print "time for one chromosome: ", endTime - startTime
		# 		startTime = time.time()
		# 		
		# 		#check if snv start and end are within gene boundaries
		# 		
		# 		snvSubset = snvData[snvData[:,0] == gene.chromosome]
		# 		
		# 		previousChr = gene.chromosome
		# 		
		# 	startMatchesStart = snvSubset[:,1] < gene.end
		# 	startMatchesEnd = snvSubset[:,1] > gene.start
		# 	
		# 	startMatches = startMatchesStart * startMatchesEnd
		# 	
		# 	endMatchesStart = snvSubset[:,2] < gene.end
		# 	endMatchesEnd = snvSubset[:,2] > gene.start
		# 	
		# 	endMatches = endMatchesStart * endMatchesEnd
		# 	
		# 	allMatches = startMatches * endMatches
		# 	
		# 	snvsOverlappingGene = snvSubset[allMatches]
		# 	gene.setSNVs(snvsOverlappingGene)
		
		# 
		#### Some testing to link SNVs to eQTLs quickly
		#Do a pre-filtering to have a much shorter list of SNVs and make overlapping faster, there wil be many SNVs that never overlap any eQTL so we don't need to look at these. 
		
		#1. Define ranges for the eQTLs
		#2. Determine the boundaries from the range
		ranges = []
		boundaries = []
		# for eQTL in eQTLData:
		# 	#eQTLRange = (eQTL.start, eQTL.end)
		# 	#ranges.append(eQTLRange)
		# 	boundaries.append(eQTL.start)
		# 	boundaries.append(eQTL.end)
		
		#3. Obtain the SNV coordinates as one list
		snvEnds = []
		snvStarts = []
		for snv in snvData:
			snvStarts.append(snv[1])
			snvEnds.append(snv[2])
			
		boundaries = np.array(boundaries)
		snvStarts = np.array(snvStarts)
		snvEnds = np.array(snvEnds)
		
		#Do the overlap to get the eQTLs that have overlap with ANY SNV on ANY chromosome. From here we can further subset.
		startTime = time.time()
		
		
		#I didn't document this and then forgot what it does... 
		# startOverlaps = np.where(np.searchsorted(boundaries, snvStarts, side="right") %2)[0]
		# endOverlaps = np.where(np.searchsorted(boundaries, snvEnds, side="right") %2)[0]
		# 
		# allEQTLOverlappingSNVs = snvData[np.union1d(startOverlaps, endOverlaps)]
		# 
		#Remove any overlap between these sets
		
		
		#Then repeat filtering for genes and TADs.
		
		geneBoundaries = []
		leftTADBoundaries = []
		rightTADBoundaries = []
		for gene in genes:
			
			geneBoundaries.append(gene[1])
			geneBoundaries.append(gene[2])
			# 
			# if gene[3].leftTAD is not None:
			# 
			# 	leftTADBoundaries.append(gene[3].leftTAD.start)
			# 	leftTADBoundaries.append(gene[3].leftTAD.end)
			# 
			# if gene[3].rightTAD is not None:
			# 	
			# 	rightTADBoundaries.append(gene[3].rightTAD.start)
			# 	rightTADBoundaries.append(gene[3].rightTAD.end)
		
			
		startOverlaps = np.where(np.searchsorted(geneBoundaries, snvStarts, side="right") %2)[0]
		endOverlaps = np.where(np.searchsorted(geneBoundaries, snvEnds, side="right") %2)[0]
		
		allGeneOverlappingSNVs = snvData[np.union1d(startOverlaps, endOverlaps)]
		
		# startOverlaps = np.where(np.searchsorted(leftTADBoundaries, snvStarts, side="right") %2)[0]
		# endOverlaps = np.where(np.searchsorted(leftTADBoundaries, snvEnds, side="right") %2)[0]
		# 
		# allLeftTADOverlappingSNVs = snvData[np.union1d(startOverlaps, endOverlaps)]
		# 
		# startOverlaps = np.where(np.searchsorted(rightTADBoundaries, snvStarts, side="right") %2)[0]
		# endOverlaps = np.where(np.searchsorted(rightTADBoundaries, snvEnds, side="right") %2)[0]
		# 
		# allRightTADOverlappingSNVs = snvData[np.union1d(startOverlaps, endOverlaps)]
		# 
		
		
		startTime = time.time()
		previousChr = None
		for gene in genes[:,3]:
			
			#Make a subset of SNVs on the right chromosome
			
			if gene.chromosome != previousChr:
				
				endTime = time.time()
				
				print "new chromosome: ", gene.chromosome
				print "time for one chromosome: ", endTime - startTime
				startTime = time.time()
				# matchingChrInd = snvData[:,0] == str(gene.chromosome)
		
				# snvSubset = snvData[matchingChrInd]
				
				#Make the chr subsets for each element type
				# matchingChrIndEQTL = allEQTLOverlappingSNVs[:,0] == str(gene.chromosome)
				matchingChrIndGenes = allGeneOverlappingSNVs[:,0] == str(gene.chromosome)
				# matchingChrIndLeftTADs = allLeftTADOverlappingSNVs[:,0] == str(gene.chromosome)
				# matchingChrIndRightTADs = allRightTADOverlappingSNVs[:,0] == str(gene.chromosome)
				
				# eQTLSNVSubset = allEQTLOverlappingSNVs[matchingChrIndEQTL]
				geneSNVSubset = allGeneOverlappingSNVs[matchingChrIndGenes]
				# leftTADSNVSubset = allLeftTADOverlappingSNVs[matchingChrIndLeftTADs]
				# rightTADSNVSubset = allRightTADOverlappingSNVs[matchingChrIndRightTADs]
				
				previousChr = gene.chromosome
			# 
			# if np.size(snvSubset) < 1:
			# 	continue #no need to compute the distance, there are no TADs on this chromosome
			# 
			#Make a smaller subset for the interval. Is this speeding up the code?
			
			#Search through blocks instead of doing the overlap on the whole set at once.
		
			
		
		
			#Search through this smaller block for the gene, TADs and eQTLs at once
			geneStartMatches = gene.start <= geneSNVSubset[:,2]
			geneEndMatches = gene.end >= geneSNVSubset[:,1]
		
			geneMatches = geneStartMatches * geneEndMatches #both should be true, use AND for concatenating
		
			#Get the SNV objects of the overlapping SNVs
			
			snvsOverlappingGenes = geneSNVSubset[geneMatches]
			#snvsOverlappingGenes = snvSubset[geneMatches]
			
			#Get the SV objects and link them to the gene
			gene.setSNVs(snvsOverlappingGenes)
			# 
			# if gene.leftTAD != None:
			# 	
			# 	leftTADStartMatches = gene.leftTAD.start <= leftTADSNVSubset[:,2]
			# 	leftTADEndMatches = gene.leftTAD.end >= leftTADSNVSubset[:,1]
			# 	
			# 	
			# 	leftTADMatches = leftTADStartMatches * leftTADEndMatches
			# 	
			# 	snvsOverlappingLeftTAD = leftTADSNVSubset[leftTADMatches]
			# 	
			# 	gene.leftTAD.setSNVs(snvsOverlappingLeftTAD)
			# 
			# if gene.rightTAD != None:
			# 	
			# 	
			# 	rightTADStartMatches = gene.rightTAD.start <= rightTADSNVSubset[:,2]
			# 	rightTADEndMatches = gene.rightTAD.end >= rightTADSNVSubset[:,1]
			# 	 
			# 	rightTADMatches = rightTADStartMatches * rightTADEndMatches
			# 
			# 	snvsOverlappingRightTAD = rightTADSNVSubset[rightTADMatches]
			# 	gene.rightTAD.setSNVs(snvsOverlappingRightTAD)
			# 
			# #Check which SVs overlap with the eQTLs
			# 
			# #Repeat for eQTLs. Is the gene on the same chromosome as the eQTL? Then use the above chromosome subset.
			# 
			# geneEQTLs = gene.eQTLs
			# 
			# for eQTL in geneEQTLs: #only if the gene has eQTLs
			# 	
			# 	startMatches = eQTL.start <= eQTLSNVSubset[:,2]
			# 	endMatches = eQTL.end >= eQTLSNVSubset[:,1]
			# 	
			# 	allMatches = startMatches * endMatches
			# 	
			# 	
			# 	snvsOverlappingEQTL = eQTLSNVSubset[allMatches]
			# 
			# 	
			# 	eQTL.setSNVs(snvsOverlappingEQTL)
			# 	
		#Quick test to get all genes with SNVs
		# dummyOut = '../../data/Genes/genesWithSNVs.txt'
		# with open(dummyOut, 'w') as outF:
		# 	
		# 	for gene in genes[:,3]:
		# 		if gene.SNVs is not None and len(gene.SNVs) > 0:
		# 			outF.write(gene.name)
		# 			outF.write("\n")
		# exit()
		dummyOut = '../../data/Genes/genesWithSNVs_counts.txt'
		with open(dummyOut, 'w') as outF:
			
			for gene in genes[:,3]:
				if gene.SNVs is not None and len(gene.SNVs) > 0:
					outF.write(gene.name + "\t" + str(len(gene.SNVs)))
					outF.write("\n")
		exit()			
			#Interactions have not been implemented for SNVs