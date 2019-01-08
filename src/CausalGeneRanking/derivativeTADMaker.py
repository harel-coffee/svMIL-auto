"""
	The goal of this script is to take a set of SVs as input, and then for each of these determine what the derivate TADs are on the part of the genome that they disrupt.
	
	For duplications, inversions, and translocations. 


"""

import re
import numpy as np
from tad import TAD
from sv import SV

class DerivativeTADMaker:
	
	
	def __init__(self, svData, genes, tadData, genome):
		
		
		self.linkSVEffectsToGenes(svData, genes, tadData, genome)
		
	
	def linkSVEffectsToGenes(self, svData, genes, tadData, genome):
		
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
			else;
				self.determineDerivativeTADs(sv, tadData, genome, "dup")
			
			
		
			
		
		return 0
		
		
		
	def determineDerivativeTADs(self, svData, tadData, genome, svType):	
	
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
		# 
		# # tad1 = ["chr1", 1, 100, TAD("chr1", 1, 100)]
		# # tad2 = ["chr1", 101, 200, TAD("chr1", 101, 200)]
		# # tad3 = ["chr1", 201, 300, TAD("chr1", 201, 300)]
		# # 
		# # tadData = np.array([tad1, tad2, tad3], dtype="object")
		# 
		# #tad1 = ["chr1", 1, 200, TAD("chr1", 1, 200)]
		# #tadData = np.array([tad1], dtype="object")
		# 
		# 
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tadData = np.array([tad1], dtype="object")
		# 
		# svData = ["chr1", 50, 50, "chr1", 250, 250, "dummy", "protective against cancer yay", SV("chr1", 50, 50, "chr1", 250, 250, "dummy", "protective against cancer yay", "duplication")]
		# 	
			
		
		### DUPLICATION ###
		if svType == "dup":
			
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
				
				for tadInd in range(1, len(filteredTads)-1):
					followingTads.append(filteredTads[tadInd])
				
			
			
				#Make the final affected TAD
				#The boundary starts at the end of the last appended TAD.
				#Then the duplication is the first part of the TAD.
				#So, the end of this TAD is the last TAD boundary + (end of duplication - last TAD boundary) + (end of original TAD - duplication end) + duplication size
				newLastTadStart = followingTads[len(followingTads)-1][2] + (svData[5] - svData[1])
				newLastTadEnd = followingTads[len(followingTads)-1][2] + (svData[5] - followingTads[len(followingTads)-1][2]) + (filteredTads[len(filteredTads)-1][2] - svData[5]) + (svData[5] - svData[1])
				
				newLastTad = TAD(svData[0], newLastTadStart, newLastTadEnd)
				
				
				#Assign the gained elements to the TAD.
				#1. Get all eQTLs that are until the SV (unaffected) in both TADs. 
				firstTadInteractions = filteredTads[0][3].getElementsByRange(filteredTads[0][1], svData[1])
				lastTadInteractions = filteredTads[len(filteredTads)-1][3].getElementsByRange(svData[5], filteredTads[len(filteredTads)-1][2])
				
				#Assign the elements to the new TADs in the right order.
				#The first TAD gets the eQTLs within the SV of the last TAD.
				#The last TAD gets the eQTLs within the SV of the last TAD.
				
				svInteractionsFirstTad = filteredTads[0][3].getElementsByRange(svData[1], filteredTads[0][2])
				svInteractionsLastTad = filteredTads[len(filteredTads)-1][3].getElementsByRange(filteredTads[len(filteredTads)-1][2], svData[5])
				
				newTad1.eQTLInteractions = svInteractionsFirstTad + firstTadInteractions
				newLastTad.eQTLInteractions = svInteractionsLastTad + lastTadInteractions #Actually, this is the original TAD!
				
				#Determine the gains for every gene. Also for the copied TADs, there are now multiple of these genes. 
				
				#For the new TADs, this is the same principle as for the eQTLs.
				#For the duplicated TADs, we can do * 2 of the elements
				
				#For TAD 1, the first part of C can interact with the second half of A.
				svGenesFirstTad = filteredTads[0][3].getGenesByRange(svData[1], filteredTads[0][2])
				svGenesLastTad = filteredTads[len(filteredTads)-1][3].getGenesByRange(filteredTads[len(filteredTads)-1][2], svData[5])
				
				for gene in svGenesFirstTad:
					gene.addGainedEQTLs(svInteractionsLastTad, svData[7])
				for gene in svGenesLastTad:
					gene.addGainedEQTLs(svInteractionsFirstTad, svData[7])
				
				#The last TAD remains the same overall.
				#Only the TADs in the middle are duplicated.
				for tad in followingTads:
					for gene in tad[3].genes:
						for sample in gene.gainedEQTLs:
							gene.addGainedEQTLs(gene.gainedEQTLs[sample], sample)
						#gene.addGainedEQTLs(gene.gainedEQTLs, svData[7]) #just duplicate the interactions
						
				
				
				
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
					
					#The first part of the TAD will gain interactions from the SV in the bin on the left.
					#The genes in the SV from the bin will gain interactions from the first part of the TAD.
					#The genes in the first part of the TAD will gainn interactions from the SV from the bin. 
					
					#The rest of the TAD remains in tact and does not lose anything. 
					
					genomicBin = genome.collectGenomicBin(svData[0], svData[5], filteredTads[0][1])
					
					svInteractionsBin = genomicBin[3].getElementsByRange(svData[1], filteredTads[0][2])
					svGenesBin = genomicBin[3].getGenesByRange(svData[1], filteredTads[0][2])
					
					svInteractionsFirstTad = filteredTads[0][3].getElementsByRange(filteredTads[0][1], svData[5])
					svGenesFirstTad = filteredTads[0][3].getGenesByRange(filteredTads[0][1], svData[5])
					
					for gene in svGenesBin:
						#These genes gain interactions from the TAD.
						gene.addGainedEQTLs(svInteractionsFirstTad, svData[7])
						
					
					for gene in svGenesFirstTad:
						#These genes gain interactions from the SV part in the bin.
						gene.addGainedEQTLs(svInteractionsBin, svData[7])
						
					
					
				else:	
					#In this first case, we only have a new part after the original TAD, where genes can gain. there are no losses.
					
					#First part of DUP until TAD boundary is in an extra TAD. All genes in there get eQTLs in there.
					#The second part of the DUP is duplicated into a bin, so all genes in the SV and in the bin gain eQTLs in the bin and in the SV. 
					
					newTad1Start = filteredTads[0][2]
					#The end of the new TAD1 is the leftmost TAD end + leftmostTAD end - duplication start
					newTad1End = filteredTads[0][2] + (filteredTads[0][2] - svData[1])
					
					#Get all interactions from the start of the SV until the end of the first TAD.
					svInteractionsFirstTad = filteredTads[0][3].getElementsByRange(svData[1], filteredTads[0][2])
					svGenesFirstTad = filteredTads[0][3].getGenesByRange(svData[1], filteredTads[0][2])
					
					#Get the bin and the interactions inside the SV part specifically
					genomicBin = genome.collectGenomicBin(svData[0], filteredTads[0][2], svData[5])
					
					svInteractionsBin = genomicBin[3].getElementsByRange(svData[1], filteredTads[0][2])
					svGenesBin = genomicBin[3].getGenesByRange(svData[1], filteredTads[0][2])
					
					#Get the genes in the bin outside of the SV
					#Get the interactions in the bin outside of the SV
					interactionsBin = genomicBin[3].getElementsByRange(svData[5], genomicBin[2])
					genesBin = genomicBin[3].getGenesByRange(svData[5], genomicBin[2])
					
					for gene in svGenesFirstTad:
						#Each gene in this bin gets all eQTLs that are within the SV.
						gene.addGainedEQTLs(svInteractionsFirstTad, svData[7])
						
					
					for gene in svGenesBin:
						#Each gene here gains eQTLs from outside of the SV in the bin.
						gene.addGainedEQTLs(interactionsBin, svData[7])
	
						
					for gene in genesBin:
						#Each gene here gains eQTLs from inside the SV.
						gene.addGainedEQTLs(svInteractionsBin, svData[7])
						
					
					
				
				
			

