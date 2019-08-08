from __future__ import absolute_import
from sv import SV
import numpy as np

class SVShuffler:
	
	
	#By default set the hg19 coordinates for shuffling SVs
	def __init__(self):
			
		###Make this a setting later
		hg19CoordinatesFile = "../../data/chromosomes/hg19Coordinates.txt"
		
		self.hg19Coordinates = dict()
		with open(hg19CoordinatesFile, 'r') as f:
			lineCount = 0
			for line in f:
				line = line.strip()
				
				if lineCount < 1: #skip header
					lineCount += 1
					continue
				
				splitLine = line.split("\t")
				
				chromosome = splitLine[0]
				end = splitLine[3]
				
				self.hg19Coordinates[chromosome] = int(end)
	
	def shuffleSVs(self, svData):
		
		shuffledSvs = []
		
		for sv in svData:
	
			#1. Get the min max bounds for the chromosomes (2 in case of translocation)
			#The minimum start coordinate is just 1, the maximum start is the length-start
			#The end coordinate is just start + length
			
			chromosome1 = sv[0]
			chromosome2 = sv[3]
			
			chr1Length = self.hg19Coordinates[chromosome1]
			chr2Length = self.hg19Coordinates[chromosome2]
			
			start1 = int(sv[1])
			start2 = int(sv[2])
			
			startDifference = start2 - start1
			
			end1 = int(sv[4])
			end2 = int(sv[5])
			
			endDifference = end2 - end1
			
			minimumStart1 = 1	
			maximumStart1 = chr1Length - start1
			
			#If the start position is after the bounds of the chromosome, the maximum start should be the length of the chromosome - the length of the SV.
		
			if maximumStart1 < 0:
				maximumStart1 = chr1Length - svLength
		
			newStart1 = random.randint(minimumStart1, maximumStart1)
			newStart2 = newStart1 + startDifference
			
			#If the chromosomes are the same, use the maximum start for chr1. Otherwise use different chromosomes.
			if chromosome1 == chromosome2:
			
				#difference between end1 and start 1
				svLength = end1 - start1
					
				newEnd1 = newStart1 + svLength
				newEnd2 = newEnd1 + endDifference
				
			else: #If the chromosomes are not the same, the start and end coordinates do not necessarily need to be equidistant. Here we can randomly sample the start on both chromosomes 		
				maximumStart2 = chr2Length - start2
				
				#If the chr2 start is outside of the bounds, we only need to subtract the difference between the positions on chromosome 2. 
				if maximumStart2 < 0:
					maximumStart2 = chr2Length - endDifference
				
				#The end coordinate here is actually the start coordinate for chromosome 2.
				newEnd1 = random.randint(minimumStart1, maximumStart2) #chr2 has the same start coordinate as chr 1
				newEnd2 = newEnd1 + endDifference
				
			#Sample name and cancer type can be copied from the other SV. Chromosome information also remains the same.
			#Keep in mind that the sample name and cancer type are reversed in the SV object beause that order makes more sense. 
			newSvObj = SV(chromosome1, newStart1, newStart2, chromosome2, newEnd1, newEnd2, sv[7], sv[6], sv[8].svType)
			newSv = [chromosome1, newStart1, newStart2, chromosome2, newEnd1, newEnd2, sv[6], sv[7], newSvObj]	
			shuffledSvs.append(newSv)	
		
		shuffledSvs = np.array(shuffledSvs, dtype="object")	
		
		return shuffledSvs
