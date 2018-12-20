"""
	The goal of this script is to take a set of SVs as input, and then for each of these determine what the derivate TADs are on the part of the genome that they disrupt.
	
	For duplications, inversions, and translocations. 


"""

import re
import numpy as np

class DerivativeTADMaker:
	
	
	def __init__(self, svData, genes, tadData):
		
		
		self.linkSVEffectsToGenes(svData, genes, tadData)
		
	
	def linkSVEffectsToGenes(self, svData, genes, tadData):
		
		"""
			For every SV, determine the type. If this is an inversion or duplication, we only need this particular SV.
			If the SV is a translocation, we collect a set of SVs that 'work together'. Based on a window around their start/end positions.
			For the collected SV or SVs, we first make derivative TADs. 
			Then after these derivative TADs have been made, we go through the genes that are present within these new TADs, and add the newly gained or remove newly lost genomic elements for these genes.
			The updated genes are returned and will later be used to determine channels of gains/loss of genomic elements. 
		"""
		
		for sv in svData:
			# 
			# typeMatch = re.search("intra", sv[8].svType, re.IGNORECASE) #using 'chr' will match intra, inter
			# if typeMatch is not None: #skip the most complicated kind for now
			# 	continue
			
			typeMatch = re.search("dup", sv[8].svType, re.IGNORECASE)
			if typeMatch is None:
				continue
			
			self.determineDerivativeTADs(sv, tadData)
			
			
		exit()
			
		
		return 0
		
		
		
	def determineDerivativeTADs(self, svData, tadData):	
	
		"""
			Given an SV or a set of SVs, depending on the type of SVs, we compute how the affected region of the genome will look after the SV.
			We then make a set of new TAD objects that are located next to each other, and update all elements that are now inside these new/affected TADs. 
		"""
	
		#1. Get all TADs that are affected by the SV.
		#2. String all TADs together in their new formation
			#- Inversion: only 2 TADs need to be re-modeled. Make new objects for both. Assign to the first all elements until the SV start, and everything from the TAD boundary to the SV end. Vice versa for the other TAD. 
			#- Duplication: Get all TADs that are covered by the duplication. Keep the elements left and right of the boundary within duplication separate. Start from the left TAD. This remains the same. Get the TAD on the right side
			#  of the duplication.
			#  Get all the TAD boundaries covered by the duplication. Place them in order. The genomic positions should be correct. So within the right TAD, we take the insert position. We attach everything within the duplication to that.
			#  Also make sure that all elements are added to the derivative TADs correctly.
			
			
		
		### DUPLICATION ###
		
		#1. Determine which TADs are involved in the duplication (only the outmost 2 are affected, the rest can be kept in tact)
		tadChrSubsetInd = svData[0] == tadData[:,0]
		tadChrSubset = tadData[tadChrSubsetInd]
		
		startMatches = svData[1] < tadChrSubset[:,2]
		endMatches = svData[2] > tadChrSubset[:,1]
		
		matches = startMatches * endMatches

		matchingTads = tadChrSubset[matches]
		
		#Remove all matches where the SV is exclusively within a TAD
		filteredTads = []
		for tad in matchingTads:
			
			if svData[1] > tad[1] and svData[2] < tad[2]:
				continue
			filteredTads.append(tad)
		
		print np.array(filteredTads)
		
		#2. Make the derivative positions for the TAD boundaries
		
		#
		
		
		
		
		
		
	
	
		return 0

