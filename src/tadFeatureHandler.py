#!/usr/bin/env python

######### Documentation #########

__author__ = "Marleen Nieboer"
__credits__ = []
__maintainer__ = "Marleen Nieboer"
__email__ = "m.m.nieboer@umcutrecht.nl"
__status__ = "Development"


### Imports ###

import sys
sys.path.append('settings/') #Add settings path
import settings
import numpy as np
import time

### Code ###

class TADFeatureHandler:
	"""
		Handler for annotating regions with TAD-related information. Currently, it annotates with the following features:
			- number of overlapping TADs
	"""
	
	def annotate(self, regions, tadData, enabledFeatures = None):
		"""
			The annotate function is mandatory. This function should be called from a database handler, such as the FlatFileDb class. It should have one return value, which contains annotations.
			
			Parameters:
				regions (numpy matrix): a Numpy matrix of dtype 'object' with on the columns chr 1, s1, e1, chr2, s2, e2, and the regions on the rows.
				tadData (numpy matrix): matrix with columns chr, start, end, with each TAD on the rows. 
				enabledFeatures (list): optional. Not yet implemented. The idea is that a list can be provided that contains the names of features that are enabled in the settings. If necessary, the listed features
				can be excluded from computation to save computational time. It can be easily implemented by having an if statement checking for that in this function.
				
			Returns:
				annotations (dict): a dictionary with annotations. {'featureName' : [values]}, where the values are in the same order as the provided regions. The value must be a list to guarantee that data is
				written to the output file in the same manner, even if it has only 1 value. 
		"""
		annotations = self.computeTADFeatures(regions, tadData)
		
		return annotations
	
	def computeTADFeatures(self, regions, tadData):
		"""
			Function that will call individual functions for different features related to TADs.
			
			Parameters:
				regions (numpy matrix): a Numpy matrix of dtype 'object' with on the columns chr 1, s1, e1, chr2, s2, e2, and the regions on the rows.
				tadData (numpy matrix): matrix with columns chr, start, end, with each TAD on the rows.
				
			Returns:
				annotations (dict): a dictionary with annotations. Contains the features:
					- overlappingTadCount: see computeNumberOfOverlappingTADs
		"""
		annotations = dict()
		
		#1. How many tads overlap with each region?
		#annotations['overlappingTadCount'] = self.computeNumberOfOverlappingTADs(regions, tadData)
		if settings.features['numberOfDisruptedTadBoundaries'] == True:
			annotations['overlappingTadBoundaries'] = self.computeNumberOfOverlappingTadsByBoundaries(regions, tadData)
		
		#2. How many boundaries are disrupted by an SV?
		#- for this query, we need to check if the SV directly overlaps with either a start or end coordinate, so cases where the SV is within the TAD should not count. 
		
		return annotations
	
	def computeNumberOfOverlappingTADs(self, regions, tadData):
		"""
			Overlaps each region with all known TADs. Overlap is defined as overlap of at least 1 bp.
			For each region, chr1 and chr2 are obtained. We then subset the TAD data to look for TADs that are on chr1 and chr2.
			Then, we check each TAD in this list and count how many TADs are overlapping with the region.
			
			Parameters:
				regions (numpy matrix): a Numpy matrix of dtype 'object' with on the columns chr 1, s1, e1, chr2, s2, e2, and the regions on the rows.
				tadData (numpy matrix): matrix with columns chr, start, end, with each TAD on the rows.
				
			Returns:
				overlappingTads (list): a list in which each value is the number of TADs that region overlaps with. In the same order as the regions.
				
			TODO:
				- Fix the notations. We are no longer working with lines here, but with tadData. This could be from any source, now it is too tailored to flat files. 
		"""
		
		overlappingTads = []
		
		#Do the actual filtering. Make sure that all files have been pre-processed to save time here
		previousChr1 = None
		previousChr2 = None
		
		#len(regions)
		for lineCount in range(0, len(regions)):
			#print lineCount
			#Make sure to filter the dataset per chromosome! Good when sorted, we can keep it until we hit the next chromosome
			#The genes are not cross-chromosomal, but translocations can be, so we need to check both chromosome entries!
			
			#Obtain start and end positions and chromosome (chr 1 and 2)
			lineList = regions[lineCount,:]
			
			#We should check the chromosome of the previous line.
			#This would work nicely if the input file has been properly sorted! We probably need a pre-sorting to make this work. 
			if str(lineList[0]) != previousChr1 and str(lineList[3]) != previousChr2:
				
				#Find the two subsets that match on both chromosomes. 
				matchingChr1Ind = tadData[:,0] == 'chr' + lineList[0]
				matchingChr2Ind = tadData[:,0] == 'chr' + lineList[3] #The tad entry is here the same since the genes are only on 1 chromosome. 
				
				#It is not even necessary to make 2 lists if both chromosomes are the same, we could use a reference without re-allocating
				chr1Subset = tadData[np.where(matchingChr1Ind)]
				if lineList[0] == lineList[3]:
					chr2Subset = chr1Subset
				else:
					chr2Subset = tadData[np.where(matchingChr2Ind)]
				
				#Make sure to update the previous chromosome when it changes
				previousChr1 = 'chr' + str(lineList[0])
				previousChr2 = 'chr' + str(lineList[3])
			
			if np.size(chr1Subset) < 1 and np.size(chr2Subset) < 1:
				overlappingTads.append(0) #otherwise we will not have one annotation per SV
				continue #no need to compute the distance, there are no genes on these chromosomes
			
			startTime = time.time()
			#Now compute how many TADs overlap
			#The start should be before the end of the tad, and the end after the start of the tad. This detects 4 scenarios of overlap (currently without using the overlap parameter).
			
			#If the chromosomes are the same, the start of the TAD is in s1 for chr1, and the end is in e2 for chr2.
			if lineList[0] == lineList[3]:
				start = int(lineList[1])
				end = int(lineList[5])
			else: #otherwise, we use the positions for chromosome 1 only
				start = int(lineList[1])
				end = int(lineList[2])
				
			startOverlapChr1 = start < chr1Subset[:,2] #I'm not sure if these indices are correct? 0 is chromosome right? 
			endOverlapChr1 = end > chr1Subset[:,1]
			
			#Now find where both of these are true (multiply because both conditions need to be true)
			overlappingSvsChr1 = startOverlapChr1 * endOverlapChr1
			
			#Overlap chr2 as well. #If the chromosomes are the same, the start of the TAD is in s1 for chr1, and the end is in e2 for chr2.
			if lineList[0] == lineList[3]:
				start = int(lineList[1])
				end = int(lineList[5])
			else: #otherwise, we use the positions for chromosome 2 only
				start = int(lineList[4])
				end = int(lineList[5])
			
			startOverlapChr2 = start < chr2Subset[:,2] #I'm not sure if these indices are correct? 0 is chromosome right? 
			endOverlapChr2 = end > chr2Subset[:,1]
			
			#Now find where both of these are true (multiply because both conditions need to be true)
			overlappingSvsChr2 = startOverlapChr2 * endOverlapChr2
			
			#Get the indices where both the start and end match the overlap criteria
			overlappingSvIndsChr1 = np.argwhere(overlappingSvsChr1 == True)
			overlappingCountChr1 = np.size(overlappingSvIndsChr1)
			overlappingSvIndsChr2 = np.argwhere(overlappingSvsChr2 == True)
			overlappingCountChr2 = np.size(overlappingSvIndsChr2)
			
			overlappingCount = overlappingCountChr1 + overlappingCountChr2
			
			overlappingTads.append(overlappingCount)
			
		
		return overlappingTads
	

	def computeNumberOfOverlappingTadsByBoundaries(self, regions, tadData):
		"""
			Function to compute how many TADs are overlapped on their boundaries (if both, counts as 1) by a variant. This feature should be applicable to both SVs and SNVs.
			
			Parameters:
				regions (numpy matrix): a Numpy matrix of dtype 'object' with on the columns chr 1, s1, e1, chr2, s2, e2, and the regions on the rows.
				tadData (numpy matrix): matrix with columns chr, start, end, with each TAD on the rows.
				
			Returns:
				overlappingTadBoundaries (list): a list in which each value is the number of TADs in which the region overlaps with the TAD boundary. In the same order as the regions.
			
		
		"""
		
		overlappingTadBoundaries = []
		
		#Do the actual filtering. Make sure that all files have been pre-processed to save time here
		previousChr1 = None
		previousChr2 = None
		
		#len(regions)
		for lineCount in range(0, len(regions)):
			#print lineCount
			#Make sure to filter the dataset per chromosome! Good when sorted, we can keep it until we hit the next chromosome
			#The genes are not cross-chromosomal, but translocations can be, so we need to check both chromosome entries!
			
			#Obtain start and end positions and chromosome (chr 1 and 2)
			lineList = regions[lineCount,:]
			
			#We should check the chromosome of the previous line.
			#This would work nicely if the input file has been properly sorted! We probably need a pre-sorting to make this work. 
			if str(lineList[0]) != previousChr1 and str(lineList[3]) != previousChr2:
				
				#Find the two subsets that match on both chromosomes. 
				matchingChr1Ind = tadData[:,0] == 'chr' + lineList[0]
				matchingChr2Ind = tadData[:,0] == 'chr' + lineList[3] #The tad entry is here the same since the genes are only on 1 chromosome. 
				
				#It is not even necessary to make 2 lists if both chromosomes are the same, we could use a reference without re-allocating
				chr1Subset = tadData[np.where(matchingChr1Ind)]
				if lineList[0] == lineList[3]:
					chr2Subset = chr1Subset
				else:
					chr2Subset = tadData[np.where(matchingChr2Ind)]
				
				#Make sure to update the previous chromosome when it changes
				previousChr1 = 'chr' + str(lineList[0])
				previousChr2 = 'chr' + str(lineList[3])
			
			if np.size(chr1Subset) < 1 and np.size(chr2Subset) < 1:
				overlappingTadBoundaries.append(0) #otherwise we will not have one annotation per SV
				continue #no need to compute the distance, there are no genes on these chromosomes
			
			#Now compute how many TADs overlap
			#The start should be before the end of the tad, and the end after the start of the tad. This detects 4 scenarios of overlap (currently without using the overlap parameter).
			
			#If the chromosomes are the same, the start of the TAD is in s1 for chr1, and the end is in e2 for chr2.
			if lineList[0] == lineList[3]:
				start = int(lineList[1])
				end = int(lineList[5])
			else: #otherwise, we use the positions for chromosome 1 only
				start = int(lineList[1])
				end = int(lineList[2])
			
			
			#We check if the start position of the TAD (in chrSubset) is overlapped by the start and end of the SV on this chromosome.
			startOverlapChr1StartTad = start <= chr1Subset[:,1]
			endOverlapChr1StartTad = end >= chr1Subset[:,1]
			#Repeat the boundary overlap for the end position. 
			startOverlapChr1EndTad = start <= chr1Subset[:,2]
			endOverlapChr1EndTad = end >= chr1Subset[:,2]
			
			#Both start and the end must overlap the start/end of the TAD for a disruption to be true
			startOverlapChr1 = startOverlapChr1StartTad * endOverlapChr1StartTad
			endOverlapChr1 = startOverlapChr1EndTad * endOverlapChr1EndTad
			
			#Now find where either of these values is true
			overlappingSvsChr1 = startOverlapChr1 + endOverlapChr1
			
			#Overlap chr2 as well. #If the chromosomes are the same, the start of the TAD is in s1 for chr1, and the end is in e2 for chr2.
			if lineList[0] == lineList[3]:
				start = int(lineList[1])
				end = int(lineList[5])
			else: #otherwise, we use the positions for chromosome 2 only
				start = int(lineList[4])
				end = int(lineList[5])
			
			
			#We check if the start position of the TAD (in chrSubset) is overlapped by the start and end of the SV on this chromosome.
			startOverlapChr2StartTad = start <= chr2Subset[:,1]
			endOverlapChr2StartTad = end >= chr2Subset[:,1]
			#Repeat the boundary overlap for the end position. 
			startOverlapChr2EndTad = start <= chr2Subset[:,2]
			endOverlapChr2EndTad = end >= chr2Subset[:,2]
			
			startOverlapChr2 = startOverlapChr2StartTad * endOverlapChr2StartTad
			endOverlapChr2 = startOverlapChr2EndTad * endOverlapChr2EndTad
			
			#Now find where either of these is true
			overlappingSvsChr2 = startOverlapChr2 + endOverlapChr2
			
			#Get the indices where both the start and end match the overlap criteria
			overlappingSvIndsChr1 = np.argwhere(overlappingSvsChr1 == True)
			overlappingCountChr1 = np.size(overlappingSvIndsChr1)
			overlappingSvIndsChr2 = np.argwhere(overlappingSvsChr2 == True)
			overlappingCountChr2 = np.size(overlappingSvIndsChr2)
			
			overlappingCount = overlappingCountChr1 + overlappingCountChr2
			
			overlappingTadBoundaries.append(overlappingCount)
			
		
		return overlappingTadBoundaries

		
		