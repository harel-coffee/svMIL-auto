"""
	Make bins of the genome so that we can easily get all genomic elements within a specific range.
	
	Go through all TADs, these are already made objects that have elements assigned.
	Within regions that are not TADs, make a bin. Given genomic elements, add these to the correct bins. 
	
"""
from tad import TAD
import numpy as np

class Genome:
	
	bins = dict()

		
		
	def defineBins(self, genes, tadData, eQTLData):
		
		## ORIGINAL code to define a genomic bin as a region between TAD boundaries
		
		lastBoundary = 1
		
		for tad in tadData:
			
			if tad[0] not in self.bins:
				self.bins[tad[0]] = []
			
			if tad[1] > lastBoundary: #this menans that we need to introduce a bin.
				newBinObj = TAD(tad[0], lastBoundary, tad[1]) #just copy TAD. Can fix that later.
				
				#Assign genomic elements to bin
				self.assignGenomicElementsToBin(newBinObj, genes, eQTLData)
				
				newBin = [tad[0], lastBoundary, tad[1], newBinObj]
				self.bins[tad[0]].append(newBin)
			else:
				self.bins[tad[0]].append(tad)
			lastBoundary = tad[2]
				
		#Make every value a np array to easily access and overlap
		
		for binKey in self.bins:
			self.bins[binKey] = np.array(self.bins[binKey], dtype="object")
		
		
		## NEW code to divide the genomic bins into regions of max 40 KB
		
		
		#There should only be 78429 in the end. The genome has 3 billion base pairs, / 40000 is 78429 And since there are already 16747 TADs, we only need to make 61682 genomic bins. 
		# print tadData.shape
		# exit()
		# 
		# 
		# windowSize = 40000
		# lastBoundary = 1
		# totalBins = 0
		# for tad in tadData:
		# 	
		# 	if tad[0] not in self.bins:
		# 		self.bins[tad[0]] = []
		# 	
		# 	if tad[1] > lastBoundary: #this menans that we need to introduce a bin.
		# 		
		# 		# #If the region is larger than our max, divide it into multiple regions.
		# 		if tad[1] - lastBoundary > windowSize:
		# 			
		# 			
		# 			noOfBins = (tad[1] - lastBoundary) / windowSize
		# 
		# 			#print "Number of bins to make: ", noOfBins
		# 			currentPos = lastBoundary
		# 			for binInd in range(0, noOfBins):
		# 				newBinObj = TAD(tad[0], currentPos, currentPos + windowSize)
		# 				
		# 				#Assign genomic elements to bin
		# 				self.assignGenomicElementsToBin(newBinObj, genes, eQTLData)
		# 				
		# 				newBin = [tad[0], currentPos, currentPos + windowSize, newBinObj]
		# 				#print "made new bin: ", newBin
		# 				self.bins[tad[0]].append(newBin)
		# 				
		# 				currentPos += windowSize
		# 				totalBins += 1
		# 		# 	print "Added bins: ", totalBins
		# 		# 	exit()
		# 		else:
		# 			
		# 			newBinObj = TAD(tad[0], lastBoundary, tad[1]) #just copy TAD. Can fix that later.
		# 			
		# 			#Assign genomic elements to bin
		# 			self.assignGenomicElementsToBin(newBinObj, genes, eQTLData)
		# 			
		# 			newBin = [tad[0], lastBoundary, tad[1], newBinObj]
		# 			self.bins[tad[0]].append(newBin)
		# 			totalBins += 1
		# 	else:
		# 		
		# 		self.bins[tad[0]].append(tad)
		# 		totalBins += 1
		# 		
		# 	lastBoundary = tad[2]
		# 		
		# #Make every value a np array to easily access and overlap
		# 
		# for binKey in self.bins:
		# 	self.bins[binKey] = np.array(self.bins[binKey], dtype="object")
		# print "total bins: ", totalBins
		
	
	def assignGenomicElementsToBin(self, newBin, genes, eQTLData):
		
		
			
		
		#Subset the right chromosome
		chrGenes = genes[genes[:,0] == newBin.chromosome]
		
		startGeneMatches = (newBin.start < chrGenes[:,1]) * (newBin.end > chrGenes[:,1])
		endGeneMatches = (newBin.start < chrGenes[:,2]) * (newBin.end > chrGenes[:,2])
		
		matchingGenes = chrGenes[startGeneMatches + endGeneMatches,3]
		
		newBin.setGenes(matchingGenes)
		
		#Repeat for eQTLs
		chrEQTLs = eQTLData[eQTLData[:,0] == newBin.chromosome]
		
		startEQTLMatches = (newBin.start < chrEQTLs[:,1]) * (newBin.end > chrEQTLs[:,1])
		endEQTLMatches = (newBin.start < chrEQTLs[:,2]) * (newBin.end > chrEQTLs[:,2])
		
		matchingEQTLs = chrEQTLs[startEQTLMatches + endEQTLMatches, 3]
		
		
		newBin.setEQTLInteractions(matchingEQTLs)
	
	def collectGenomicBin(self, chromosome, start, end):
		"""
			Provide a chromosome and start and end, and return the bin that this range is in. Getting the specific elements should be possible through the bin itself.
		"""
		
		bins = self.bins[chromosome]
		
		matchingBinInd = (start >= bins[:,1]) * (end <= bins[:,2])
		
		if len(bins[matchingBinInd]) < 1:
			return None #Sometimes there is no bin for the eQTL. Mapping issues? 
		
		return bins[matchingBinInd][0] #In principle, there should be only one. But if the wrong coordinates are supplied, this can return multiple. 
	
	def getAllBinsInTADFormat(self):
		"""
			Get all the bins in the same format in which the TAD data is formatted (numpy array).
			Chromosome	Start	End	Object
		"""
		allBins = []
		for chromosome in self.bins:
			for region in self.bins[chromosome]:
				allBins.append(region)
		
		return np.array(allBins)
		
		
		
		