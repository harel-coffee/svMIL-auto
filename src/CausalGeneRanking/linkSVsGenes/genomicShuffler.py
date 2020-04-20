from __future__ import absolute_import
from __future__ import print_function
from sv import SV
from tad import TAD
import numpy as np

import random
import settings

class GenomicShuffler:
	"""
		The goal of this class is to provide functions to enable shuffling of genomic elements.
		
		Currently implemented:
		- Shuffling SVs
		- Shuffling TADs
	
	"""

	#By default set the hg19 coordinates for shuffling SVs
	def __init__(self):
		"""
			By default, we need the hg19 coordinates to make sure that we do not shuffle outside of the genome. 
		
		"""
		
		hg19CoordinatesFile = settings.files['hg19CoordinatesFile']
		
		self.hg19Coordinates = dict()
		with open(hg19CoordinatesFile, 'r') as f:
			lineCount = 0
			for line in f:
				line = line.strip()
				
				if lineCount < 1: #skip header
					lineCount += 1
					continue
				
				splitLine = line.split("\t")
				
				chromosome = 'chr' + splitLine[0]
				end = splitLine[3]
				
				self.hg19Coordinates[chromosome] = int(end) #use the outermost coordinate as a limit
	
	def shuffleTADs(self, tadData):
		"""
			Assign new random positions to the TADs. Similar to the shuffling of the SVs.
			The contents of the TADs should depend on where these are located on the genome, not be copied from the original TADs.
			
			We use a circular shuffling, where TADs are shifted by a specified offset to the right. This is 1 mb for now.
			If a TAD goes outside of the genome, we continue with the TAD at the beginning. This will be split into 2 TADs. 
			
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
			
			return
			shuffledTads: (numpy array) array with the TADs and their information, but then with shuffled coordinates. chr, start, end, tadObject
			
		"""
	
		shuffledTads = []
		
		#Instead of assigning random positions to the TADs, use an offset. At the end of the genome, start again at the beginning (cicular shuffling)
		offset = 100000 #Use a 1 mb offset.This is a random choice for now. 
	
		for tad in tadData:
			
			chrom = tad[0]
			
			if chrom not in self.hg19Coordinates: #skip strange chromosomes 
				continue
			
			
			chrLength = self.hg19Coordinates[chrom]
			
			tadLength = tad[2] - tad[1]

			offsetStart = tad[1] + offset
			offsetEnd = tad[2] + offset
			
			if offsetStart > chrLength: #In this case, start from the beginning of the chromosome. 
				overhang = chrLength - offsetStart
				newStart = overhang #simply count from 1 here
				newEnd = newStart + tadLength
			
				newTadObj = TAD(chrom, newStart, newEnd)
				shuffledTads.append([chrom, newStart, newEnd, newTadObj])
				continue
			
			if offsetEnd > chrLength: #In this case, split the TAD in 2 TADs.
				newTadObj = TAD(chrom, offsetStart, chrLength)
				shuffledTads.append([chrom, offsetStart, chrLength, newTadObj])
				
				overhang = chrLength - offsetStart
				newStart = 1
				newEnd = overhang #the remaining part of the TAD
			
				newTadObj = TAD(chrom, newStart, newEnd)
				shuffledTads.append([chrom, newStart, newEnd, newTadObj])
				continue
			
			newTadObj = TAD(chrom, offsetStart, offsetEnd)
			shuffledTads.append([chrom, offsetStart, offsetEnd, newTadObj])

		#Sort the TADs
		shuffledTads = np.array(shuffledTads, dtype="object")
		ind = np.lexsort((shuffledTads[:,1], shuffledTads[:,0]))

		shuffledTads = shuffledTads[ind,:]
		
		#Write the shuffled TADs to a file to check in IGV
		outFile = "shuffledTads.bed"
		with open(outFile, 'wb') as outF:
			for tad in shuffledTads:
				outF.write(tad[0] + "\t" + str(tad[1]) + "\t" + str(tad[2]) + "\n")

		return shuffledTads
		
	
	def shuffleSVs(self, svData):
		"""
			Shuffle the coordinates of the SVs across the genome. SVs are given random coordinates, but we keep their length and original chromosome.
			We make sure that the new coordinates are not outside of the genome. 
			
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
		
			return
			shuffledSvs: (numpy array) array with the SVs and their information, but with random coordinates. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
		"""
		
		shuffledSvs = []
		for sv in svData:
			
			#1. Get the min max bounds for the chromosomes (2 in case of translocation)
			#The minimum start coordinate is just 1, the maximum start is the length-start
			#The end coordinate is just start + length
			
			chromosome1 = sv[0]
			chromosome2 = sv[3]
			
			if chromosome1 not in self.hg19Coordinates:
				print(chromosome1)
				continue
			if chromosome2 not in self.hg19Coordinates:
				print(chromosome2)
				continue
			
			chr1Length = self.hg19Coordinates[chromosome1]
			chr2Length = self.hg19Coordinates[chromosome2]
			
			start = int(sv[1])
			end = int(sv[5])
			
			minimumStart = 1	
			
			###Based on random positions
			
			#If the chromosomes are the same, use the maximum start for chr1. Otherwise use different chromosomes.
			if chromosome1 == chromosome2 and start < end:
				
				svLength = end - start
				
				maximumStart = chr1Length - svLength #do not use the end here, because that is the end of the original SV. It can be anywhere, depending on the length. 
				
				newStart = random.randint(minimumStart, maximumStart)
				newEnd = newStart + svLength
	
			elif chromosome1 != chromosome2: #If the chromosomes are not the same, the start and end coordinates do not necessarily need to be equidistant. Here we can randomly sample the start on both chromosomes 		
			
				maximumStart1 = chr1Length #use the end on the first chromosome
				newStart = random.randint(minimumStart, maximumStart1)
				
				maximumStart2 = chr2Length
				
				#The end coordinate here is actually the start coordinate for chromosome 2.
				newEnd = random.randint(minimumStart, maximumStart2) #chr2 has the same start coordinate as chr 1
			else: #intrachromosomal translocation, also randomly sample new positions here, equidistant can mean that the SV will go outside of the chromosome
				maximumStart1 = chr1Length #use the end on the first chromosome
				newStart = random.randint(minimumStart, maximumStart1)
				
				maximumStart2 = chr2Length
				
				#The end coordinate here is actually the start coordinate for chromosome 2.
				newEnd = random.randint(minimumStart, maximumStart2) #chr2 has the same start coordinate as chr 1
			
			#Sample name and cancer type can be copied from the other SV. Chromosome information also remains the same.
			#Keep in mind that the sample name and cancer type are reversed in the SV object beause that order makes more sense.
			newStart1 = newStart
			newEnd1 = newStart
			newStart2 = newEnd
			newEnd2 = newEnd
			newSvObj = SV(chromosome1, newStart1, newEnd1, sv[8].o1, chromosome2, newStart2, newEnd2, sv[8].o2, sv[7], sv[6], sv[8].svType)
			newSv = [chromosome1, newStart1, newEnd1, chromosome2, newStart2, newEnd2, sv[6], sv[7], newSvObj]	
			shuffledSvs.append(newSv)	
		
		shuffledSvs = np.array(shuffledSvs, dtype="object")	
		
		return shuffledSvs
