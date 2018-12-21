"""
	The goal of this script is to take a set of SVs as input, and then for each of these determine what the derivate TADs are on the part of the genome that they disrupt.
	
	For duplications, inversions, and translocations. 


"""

import re
import numpy as np
from tad import TAD
from sv import SV

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
		
		
		#Make some dummy TADs and SVs for testing
		
		# tad1 = ["chr1", 1, 100, TAD("chr1", 1, 100)]
		# tad2 = ["chr1", 101, 200, TAD("chr1", 101, 200)]
		# tad3 = ["chr1", 201, 300, TAD("chr1", 201, 300)]
		# 
		# tadData = np.array([tad1, tad2, tad3], dtype="object")
		# 
		# tad1 = ["chr1", 1, 200, TAD("chr1", 1, 200)]
		# tadData = np.array([tad1], dtype="object")
		#
		
		tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		tadData = np.array([tad1], dtype="object")
		
		svData = ["chr1", 50, 50, "chr1", 250, 250, "dummy", "protective against cancer yay", SV("chr1", 50, 50, "chr1", 250, 250, "dummy", "protective against cancer yay", "duplication")]
			
			
		
		### DUPLICATION ###
		
		#1. Determine which TADs are involved in the duplication (only the outmost 2 are affected, the rest can be kept in tact)
		tadChrSubsetInd = svData[0] == tadData[:,0]
		tadChrSubset = tadData[tadChrSubsetInd]
		
		startMatches = svData[1] < tadChrSubset[:,2]
		endMatches = svData[2] > tadChrSubset[:,1]
		
		matches = startMatches + endMatches #either the start or end needs to match

		matchingTads = tadChrSubset[matches]

		#Remove all matches where the SV is exclusively within a TAD
		filteredTads = []
		for tad in matchingTads:
			if svData[1] > tad[1] and svData[5] < tad[2]:
				continue
			filteredTads.append(tad)
		
		print "filtered tads: ", np.array(filteredTads)
		
		if len(filteredTads) < 1:
			return
		
		#2. Make the derivative positions for the TAD boundaries
		
		#The original TAD boundaries of the left TAD remains where it is.
		#Then get the insert point in the right TAD (dup end position)
		#Derive the position of the duplicated boundary (piece of left TAD until TAD boundary)
		
		#The start position of the new TAD is the end of the leftmost TAD in which the duplication starts.
		#If the duplication overlaps with the start of the first TAD, so there is no TAD on the left, this should be the start of rightmost tad in which the duplication ends. 
		
		
		
		
		#In case only 1 boundary is overlapped and there is no TAD to the right, the first new TAD is also the last TAD.
		
		if len(filteredTads) > 1:
			
			#The next TAD boundary is the first boundary disrupted by the duplication, so the start of the last TAD disrupted
			#The position should be the insert position (duplication end) + (TAD boundary - duplication start)
			#The TAD boundary position is the end of the TAD if the duplication overlaps with the end of the TAD, otherwise it is the start
			newTad1Start = filteredTads[len(filteredTads)-1][1]
			
			
			if svData[1] < filteredTads[0][1] and svData[5] > filteredTads[0][1]:
				newTad1End = svData[5] + (filteredTads[0][1] - svData[1])	
			else:
				newTad1End = svData[5] + (filteredTads[0][2] - svData[1])
			
			newTad1 = TAD(svData[0], newTad1Start, newTad1End)
				
			#For every other TAD overlapped by the SV (except for the first and last), simply add the original TADs.
			followingTads = []
			print "TAD in middle: ", len(filteredTads)-1
			for tadInd in range(1, len(filteredTads)-1):
				followingTads.append(filteredTads[tadInd])
				
			
			print "new TAD 1 pos: ", newTad1.start, newTad1.end
			print "appended TADs: ", followingTads
		
		
			#Make the final affected TAD
			#The boundary starts at the end of the last appended TAD.
			#Then the duplication is the first part of the TAD.
			#So, the end of this TAD is the last TAD boundary + (end of duplication - last TAD boundary) + (end of original TAD - duplication end) + duplication size
			newLastTadStart = followingTads[len(followingTads)-1][2] + (svData[5] - svData[1])
			newLastTadEnd = followingTads[len(followingTads)-1][2] + (svData[5] - followingTads[len(followingTads)-1][2]) + (filteredTads[len(filteredTads)-1][2] - svData[5]) + (svData[5] - svData[1])
			
			newLastTad = TAD(svData[0], newLastTadStart, newLastTadEnd)
			print "new last tad pos: ", newLastTadStart, newLastTadEnd
		else: #There is only 1 overlapped boundary.
			
			#If the dup overlaps with the start of the TAD, we need different coordinates than at the end.
			#For the start, the start of the new TAD is the start of the original TAD
			#The end of the first TAD is the start of the new TAD +  + (start of TAD - duplication start)
			#Then the final TAD is the end of the first TAD until the original TAD end + duplication size - original TAD start
			
			if svData[1] < filteredTads[0][1] and svData[5] > filteredTads[0][1]:
				newTad1Start = filteredTads[0][1]
				newTad1End = filteredTads[0][1] + (filteredTads[0][1] - svData[1])
				
				newLastTadStart = newTad1End
				#The TAD end is the SV end - original TAD start + original TAD end - SV end + the new TAD end
				newLastTadEnd = (svData[5] - filteredTads[0][1]) + (filteredTads[0][2] - svData[5]) + newTad1End
			else:	
				newTad1Start = filteredTads[0][2]
				#The end of the new TAD1 is the leftmost TAD end + leftmostTAD end - duplication start
				newTad1End = filteredTads[0][2] + (filteredTads[0][2] - svData[1])
				
			print "new tad 1 pos: ", newTad1Start, newTad1End
			print "last tad pos: ", newLastTadStart, newLastTadEnd
			
		
		
		
		
		#Assign the elements to the right TADs. 
		#This only needs to be done for the derivative TADs.
		#The first TAD is a combination of the first part of the last TAD (last TAD boundary -> duplication end) and the last part of the first TAD (duplication start -> end of TAD)
		#The last TAD is a combination of the first part of the last TAD and the second part of the last TAD. So just itself, I guess.
		
		#Does this still work for duplications that overlap only 1 TAD boundary? 
		
		
		
		
		
		
		
	
	
		return 0

