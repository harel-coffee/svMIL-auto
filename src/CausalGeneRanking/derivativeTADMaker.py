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
		
		allBins = genome.getAllBinsInTADFormat()
		#Combine the TADs with the bins
		combinedData = np.concatenate((tadData, allBins))
		#Should the combined data be sorted?
		print combinedData.shape
		print combinedData
		
		#Not sure if sorting here is useful, since the data would need to be sorted by chromosome as well, but since we take subsets later, it is easier to sort at that point. 
		#combinedData = combinedData[combinedData[:,1].argsort()]
		
		
		
		
		self.linkSVEffectsToGenes(svData, genes, tadData, genome)
		
	
	def linkSVEffectsToGenes(self, svData, genes, tadData, genome):
		
		"""
			For every SV, determine the type. If this is an inversion or duplication, we only need this particular SV.
			If the SV is a translocation, we collect a set of SVs that 'work together'. Based on a window around their start/end positions.
			For the collected SV or SVs, we first make derivative TADs. 
			Then after these derivative TADs have been made, we go through the genes that are present within these new TADs, and add the newly gained or remove newly lost genomic elements for these genes.
			The updated genes are returned and will later be used to determine channels of gains/loss of genomic elements. 
		"""
		print "Linking SV effects to genes"
		invCount = 0
		dupCount = 0
		
		
# 		
		
		for sv in svData:
			# 
			# typeMatch = re.search("intra", sv[8].svType, re.IGNORECASE) #using 'chr' will match intra, inter
			# if typeMatch is not None: #skip the most complicated kind for now
			# 	continue
			# 
			
			# 
			
			
			typeMatch = re.search("inv", sv[8].svType, re.IGNORECASE)
			if typeMatch is not None:
				self.determineDerivativeTADs(sv, tadData, genome, "inv")
				invCount += 1
				print "inversion count: ", invCount
		
		for sv in svData:
			
			typeMatch = re.search("dup", sv[8].svType, re.IGNORECASE)
			if typeMatch is not None:
				
				self.determineDerivativeTADs(sv, tadData, genome, "dup")
				dupCount += 1
				print "duplication count: ", dupCount
				
		
		
		print "done making derivative TADs"
		
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
		
		# svData = ["chr1", 50, 50, "chr1", 250, 250, "dummy", "protective against cancer yay", SV("chr1", 50, 50, "chr1", 250, 250, "dummy", "protective against cancer yay", "duplication")]
		# 
		# tad1 = ["chr1", 1, 100, TAD("chr1", 1, 100)]
		# tad3 = ["chr1", 200, 300, TAD("chr1", 20, 300)]
		# 
		# tadData = np.array([tad1, tad3], dtype="object")
		# 
		# svData = ["chr1", 60, 60, "chr1", 210, 210, "dummy", "protective against cancer yay", SV("chr1", 60, 60, "chr1", 210, 210, "dummy", "protective against cancer yay", "inversion")]
		# 
		# 
		
		# 	
			
		### INVERSION ###
		if svType == "inv":
			print "Subsetting TADs"
			
			#1. Get the two TADs at the start and end of the inversions (these may be genomic bins)
			
			tadChrSubsetInd = svData[0] == tadData[:,0]
			tadChrSubset = tadData[tadChrSubsetInd]
			
			
			
			#Get all TADs overlapped by the inversion
			#First get all TADs overlapped by the start of the inversion
			startMatches = svData[1] <= tadChrSubset[:,2]
			endMatches = svData[1] >= tadChrSubset[:,1]
			
			invStartMatches = startMatches * endMatches #either the start or end needs to match
			
			#Then get the TADs overlapped by the end of the inversion
			startMatches = svData[5] >= tadChrSubset[:,1]
			endMatches = svData[5] <= tadChrSubset[:,2]
			
			invEndMatches = startMatches * endMatches
			
			leftMostTad = tadChrSubset[invStartMatches]
			rightMostTad = tadChrSubset[invEndMatches]
			
				
			if len(leftMostTad) < 1 or len(rightMostTad) < 1:
				return #for now skip all inversions that do not end in a TAD on either side. 

			#These are only the cases where the inversion ends in a TAD on both sides. 
			if len(leftMostTad) > 0 and len(rightMostTad) > 0:
				print "SV ends in two TADs"
				
				
				if leftMostTad[0][1] == rightMostTad[0][1] and leftMostTad[0][2] == rightMostTad[0][2]: #skip if the SV is within a TAD entirely
					return
			
				
				leftMostTad = leftMostTad[0]
				rightMostTad = rightMostTad[0]
				
				#2. Collect all elements until the right TAD boundary inside the inversion.
				
				leftSideElements = leftMostTad[3].getElementsByRange(svData[1], leftMostTad[2]) #From the start of the inversion until the end of the left most TAD
				unaffectedElementsLeft = leftMostTad[3].getElementsByRange(leftMostTad[1], svData[1])
	
				#Also get the genes
				leftSideGenes = leftMostTad[3].getGenesByRange(svData[1], leftMostTad[2]) #From the start of the inversion until the end of the left most TAD
				unaffectedGenesLeft = leftMostTad[3].getGenesByRange(leftMostTad[1], svData[1])
				
				#3. Collect all elements from the left TAD boundary until the end of the inversion.
				
				rightSideElements = rightMostTad[3].getElementsByRange(rightMostTad[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
				
				unaffectedElementsRight = rightMostTad[3].getElementsByRange(svData[5], rightMostTad[2])
				
				rightSideGenes = rightMostTad[3].getGenesByRange(rightMostTad[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
				unaffectedGenesRight = rightMostTad[3].getGenesByRange(svData[5], rightMostTad[2])
				
				#4. Make two new TAD objects. Assign all elements until the inversion start/from the inversion end to the TADs.
				
				#The new end position of the leftmost TAD is the start of the inversion + the length of the inversion part in the right TAD (inversion end - right TAD start)
				
				# leftTadNewEnd = svData[1] + (svData[5] - rightMostTad[1])
				# 
				# #The new start position of the rightmost TAD is the inversion end - (leftmost TAD end - inversion start)
				# rightTadNewStart = svData[5] - (leftMostTad[2] - svData[1])
				# 
				# newLeftTad = TAD(svData[0], leftMostTad[1], leftTadNewEnd)
				# newRightTad = TAD(svData[0], rightTadNewStart, rightMostTad[2])
				
				#5. Assign the gains and losses to the genes
				
				# #All genes that were originally in the left TAD (outisde of the inversion) will gain elements of the right side of the inversion
				# for gene in unaffectedGenesLeft:
				# 	gene.addGainedEQTLs(rightSideElements, svData[7])
				# 
				# #All genes in the right side of the inversion will gain elements from the original left TAD.
				# for gene in rightSideGenes:
				# 	gene.addGainedEQTLs(unaffectedElementsLeft, svData[7])
				# 
				# #vice versa but then for the right TAD and right side of the inversion. 
				# for gene in unaffectedGenesRight:
				# 	gene.addGainedEQTLs(leftSideElements, svData[7])
				# 
				# for gene in leftSideGenes:
				# 	gene.addGainedEQTLs(unaffectedElementsRight, svData[7])
			elif len(leftMostTad) < 1 or len(rightMostTad) < 1:
				print "SV ends in one TAD"
				
				
				
				#Either the left side of the inversion or the right side is within a TAD, but the other part is within a genomic bin. 
				#Here make use of genomic bins instead of TADs.
				
				#1. First get the genomic bin for the non-TAD end
				#Case where the inversion ends in a TAD, but the left side is a genomic bin. 
				if len(leftMostTad) < 1:
					rightMostTad = rightMostTad[0]
					genomicBin = genome.collectGenomicBin(svData[0], svData[1], svData[2])
					
					if genomicBin == None:
						return
					
					#Collect the elements and genes that are gained and lost within the TAD or genomic bin
					
					leftSideElements = genomicBin[3].getElementsByRange(svData[1], genomicBin[2]) #From the start of the inversion until the end of the left most TAD
					unaffectedElementsLeft = genomicBin[3].getElementsByRange(genomicBin[1], svData[1])
			
					#Also get the genes
					leftSideGenes = genomicBin[3].getGenesByRange(svData[1], genomicBin[2]) #From the start of the inversion until the end of the left most TAD
					unaffectedGenesLeft = genomicBin[3].getGenesByRange(genomicBin[1], svData[1])
					
					#3. Collect all elements from the left TAD boundary until the end of the inversion.
					
					rightSideElements = rightMostTad[3].getElementsByRange(rightMostTad[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
					unaffectedElementsRight = rightMostTad[3].getElementsByRange(svData[5], rightMostTad[2])
					
					rightSideGenes = rightMostTad[3].getGenesByRange(rightMostTad[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
					unaffectedGenesRight = rightMostTad[3].getGenesByRange(svData[5], rightMostTad[2])
				#Case where the inversion starts in a TAD, but the right side is in a genomic bin.
				elif len(rightMostTad) < 1:
					leftMostTad = leftMostTad[0]
					genomicBin = genome.collectGenomicBin(svData[0], svData[4], svData[5])
					if genomicBin == None:
						return
					
					#2. Collect all elements until the right TAD boundary inside the inversion.
					
					leftSideElements = leftMostTad[3].getElementsByRange(svData[1], leftMostTad[2]) #From the start of the inversion until the end of the left most TAD
					unaffectedElementsLeft = leftMostTad[3].getElementsByRange(leftMostTad[1], svData[1])
			
					#Also get the genes
					leftSideGenes = leftMostTad[3].getGenesByRange(svData[1], leftMostTad[2]) #From the start of the inversion until the end of the left most TAD
					unaffectedGenesLeft = leftMostTad[3].getGenesByRange(leftMostTad[1], svData[1])
					
					#3. Collect all elements from the left TAD boundary until the end of the inversion.
					
					rightSideElements = genomicBin[3].getElementsByRange(genomicBin[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
					unaffectedElementsRight = genomicBin[3].getElementsByRange(svData[5], genomicBin[2])
					
					rightSideGenes = genomicBin[3].getGenesByRange(genomicBin[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
					unaffectedGenesRight = genomicBin[3].getGenesByRange(svData[5], genomicBin[2])
						
			
			
			#Assigning the gains and losses to the genes is independent of the type of inversion
			print "Copying genes and elements after SV"		
			#All genes that were originally in the left TAD (outisde of the inversion) will gain elements of the right side of the inversion
			#All unaffected genes on the left will lose the eQTLs that are in the left side of the inversion
			for gene in unaffectedGenesLeft:
				
				if gene.name == "PDE4DIP":
					print svData[0], svData[1], svData[5]
					exit()
				
				gene.addGainedEQTLs(rightSideElements, svData[7])
				gene.addLostEQTLs(leftSideElements, svData[7])
				#print "Number of gained right side elements: ", len(rightSideElements), " for genes ", len(unaffectedGenesLeft)
				
				#if svData[7] in gene.gainedEQTLs:
				#print "r: ", len(gene.gainedEQTLs[svData[7]])
				#print "l: ", len(gene.lostEQTLs[svData[7]])
			
			#All genes in the right side of the inversion will gain elements from the original left TAD.
			#All genes in the right side will lose interactions with eQTLs in the unaffected right TAD. 
			for gene in rightSideGenes:
				if gene.name == "ARID1A":
					print svData[0], svData[1], svData[5]
					exit()
				gene.addGainedEQTLs(unaffectedElementsLeft, svData[7])
				#print "Number of unaffected elements right: ", len(unaffectedElementsRight), " for genes ", len(rightSideGenes)
				gene.addLostEQTLs(unaffectedElementsRight, svData[7])
			
			#vice versa but then for the right TAD and right side of the inversion.
			#The lost eQTLs are the ones that are in the right side of the inversion
			for gene in unaffectedGenesRight:
				if gene.name == "ARID1A":
					print svData[0], svData[1], svData[5]
					exit()
				gene.addGainedEQTLs(leftSideElements, svData[7])
				#print "Number of gained right side elements 2: ", len(rightSideElements), " for genes ", len(unaffectedGenesRight)
				gene.addLostEQTLs(rightSideElements, svData[7])
			
			#The lost eQTLs are the ones that are in the unaffected original left TAD
			for gene in leftSideGenes:
				if gene.name == "ARID1A":
					print svData[0], svData[1], svData[5]
					exit()
				gene.addGainedEQTLs(unaffectedElementsRight, svData[7])
				#print "Number of unaffected left elements: ", len(unaffectedElementsLeft), " for genes ", len(leftSideGenes)
				gene.addLostEQTLs(unaffectedElementsLeft, svData[7])
			
			
				
			return
		
		### DUPLICATION ###
		if svType == "dup":
			
			
			#1. Determine which TADs are involved in the duplication (only the outmost 2 are affected, the rest can be kept in tact)
			tadChrSubsetInd = svData[0] == tadData[:,0]
			tadChrSubset = tadData[tadChrSubsetInd]
			tadChrSubset = tadChrSubset[tadChrSubset[:,1].argsort()]
			
			startMatches = svData[1] < tadChrSubset[:,2]
			endMatches = svData[5] > tadChrSubset[:,1]  ##? Is this correct?
			
			matches = startMatches * endMatches
	
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
			
			#Possible cases:
			#- The duplication ends in a TAD on both sides.
			#- The duplication ends in a genomic bin on either side (and can span multiple boundaries)
			#- The duplication is in a genomic bin on both sides (skip this case)
			
			
			if len(filteredTads) > 1: #The SV spans multiple TADs
				
				#First get all TADs overlapped by the start of the inversion
				startMatches = svData[1] <= tadChrSubset[:,2]
				endMatches = svData[1] >= tadChrSubset[:,1]
				
				dupStartMatches = startMatches * endMatches #either the start or end needs to match
				
				#Then get the TADs overlapped by the end of the inversion
				startMatches = svData[5] >= tadChrSubset[:,1]
				endMatches = svData[5] <= tadChrSubset[:,2]
				
				dupEndMatches = startMatches * endMatches
				
				leftMostTad = tadChrSubset[dupStartMatches]
				rightMostTad = tadChrSubset[dupEndMatches]
					
				if len(leftMostTad) < 1 or len(rightMostTad) < 1:
					return #For now here only use cases where the duplication ends in 2 TADs. 
				
				#Here we should also distinguish between cases where the duplications are either in a genomic bin or not.
				#Can this code be simplified? In principle, the ideas are the same, but then with either a TAD or bin. 
				
				
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
				
				if len(followingTads) > 0: #We now do not have the case where the dup crosses only 1 boundary but does end in two TADs
					
					
					
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
					# 
					# if len(firstTadInteractions) > 0:
					# 	print "First TAD has interactions"
					# if len(lastTadInteractions) > 0:
					# 	print "Last TAD has interactions"
					# 
					
					#Assign the elements to the new TADs in the right order.
					#The first TAD gets the eQTLs within the SV of the last TAD.
					#The last TAD gets the eQTLs within the SV of the last TAD.
					
					svInteractionsFirstTad = filteredTads[0][3].getElementsByRange(svData[1], filteredTads[0][2])
					svInteractionsLastTad = filteredTads[len(filteredTads)-1][3].getElementsByRange(filteredTads[len(filteredTads)-1][2], svData[5])
					
					# if len(svInteractionsFirstTad) > 0:
					# 	print "SV in first TAD has interactions"
					# if len(svInteractionsLastTad) > 0:
					# 	print "SV in last TAD has interactions"
					# 
					# 
					# newTad1.eQTLInteractions = svInteractionsFirstTad + firstTadInteractions
					# newLastTad.eQTLInteractions = svInteractionsLastTad + lastTadInteractions #Actually, this is the original TAD!
					# 
					#Determine the gains for every gene. Also for the copied TADs, there are now multiple of these genes. 
					
					#For the new TADs, this is the same principle as for the eQTLs.
					#For the duplicated TADs, we can do * 2 of the elements
					
					#For TAD 1, the first part of C can interact with the second half of A.
					svGenesFirstTad = filteredTads[0][3].getGenesByRange(svData[1], filteredTads[0][2])
					svGenesLastTad = filteredTads[len(filteredTads)-1][3].getGenesByRange(filteredTads[len(filteredTads)-1][2], svData[5])
					# 
					# 
					# if len(svGenesFirstTad) > 0:
					# 	print "SV in first TAD has genes"
					# if len(svGenesLastTad) > 0:
					# 	print "SV in last TAD has genes"
					#
					print "Number of genes to add gains: ", len(svGenesFirstTad)
					for gene in svGenesFirstTad:
						#print "adding gains from right TAD: ", len(svInteractionsLastTad)
						gene.addGainedEQTLs(svInteractionsLastTad, svData[7])
					print "(2) Number of genes to add gains: ", len(svGenesLastTad)
					for gene in svGenesLastTad:
						#print "adding gains from left TAD: ", len(svInteractionsFirstTad)
						gene.addGainedEQTLs(svInteractionsFirstTad, svData[7])
					
					#The last TAD remains the same overall.
					#Only the TADs in the middle are duplicated.
					
					#Each gene gains the eQTLs that are within the TAD that the gene is located in. 
					
					for tad in followingTads:
						for gene in tad[3].genes:
							
							#1. Get all eQTLs within this TAD
							tadEQTLs = tad[3].eQTLInteractions
							
								
							#2. Filter these for the eQTLs of the gene
							gainedEQTLs = []
							for eQTL in tadEQTLs:
								if gene in eQTL.genes:
									gainedEQTLs.append(eQTL)
							#3. Add the eQTLs to the gene for the current sample
							gene.addGainedEQTLs(gainedEQTLs, svData[7])
							if len(gainedEQTLs) > 0:
								print "adding gains from unaffected region: ", len(gainedEQTLs)
				else: #Case where the duplication crosses 1 boundary
				
					#The left TAD (A) stays the same
					#A new TAD originates in the middle, where the dup part of B can interact with the dup part of A
					#The right TAD is a dup part of B with the original B, and here everything remains the same.
					#So we only need to model the gained interactions in the middle TAD (the new TAD)
					
					
					newTad1Start = filteredTads[0][2]
					#The end of the new TAD1 is the leftmost TAD (A) end + leftmostTAD (A) end - duplication start
					newTad1End = filteredTads[0][2] + (filteredTads[0][2] - svData[1])
					
					#Get all interactions from the start of the SV until the end of the first TAD (A).
					svInteractionsFirstTad = filteredTads[0][3].getElementsByRange(svData[1], filteredTads[0][2])
					svGenesFirstTad = filteredTads[0][3].getGenesByRange(svData[1], filteredTads[0][2])
					
					#Get the second TAD (B) and the interactions inside the SV part specifically
					
					svInteractionsSecondTad = filteredTads[1][3].getElementsByRange(filteredTads[0][2], svData[5])
					svGenesSecondTad = filteredTads[1][3].getGenesByRange(filteredTads[0][2], svData[5])
					
					#Genes in the sv part of the first TAD will interact with the eQTLs in the SV of the second TAD
					
					
					for gene in svGenesFirstTad:
						#Each gene in this bin gets all eQTLs that are within the SV.
						gene.addGainedEQTLs(svInteractionsSecondTad, svData[7])
						
					for gene in svGenesSecondTad:
						#Each gene here gains eQTLs from outside of the SV in the bin.
						gene.addGainedEQTLs(svInteractionsFirstTad, svData[7])
			
					
					
				
			# else: #There is only 1 overlapped boundary.
			# 	print "1 overlapped boundary"
			# 	#If the dup overlaps with the start of the TAD, we need different coordinates than at the end.
			# 	#For the start, the start of the new TAD is the start of the original TAD
			# 	#The end of the first TAD is the start of the new TAD +  + (start of TAD - duplication start)
			# 	#Then the final TAD is the end of the first TAD until the original TAD end + duplication size - original TAD start
			# 	
			# 	if svData[1] < filteredTads[0][1] and svData[5] > filteredTads[0][1]:
			# 		newTad1Start = filteredTads[0][1]
			# 		newTad1End = filteredTads[0][1] + (filteredTads[0][1] - svData[1])
			# 		
			# 		newLastTadStart = newTad1End
			# 		#The TAD end is the SV end - original TAD start + original TAD end - SV end + the new TAD end
			# 		newLastTadEnd = (svData[5] - filteredTads[0][1]) + (filteredTads[0][2] - svData[5]) + newTad1End
			# 		
			# 		#The first part of the TAD will gain interactions from the SV in the bin on the left.
			# 		#The genes in the SV from the bin will gain interactions from the first part of the TAD.
			# 		#The genes in the first part of the TAD will gainn interactions from the SV from the bin. 
			# 		
			# 		#The rest of the TAD remains in tact and does not lose anything. 
			# 		
			# 		genomicBin = genome.collectGenomicBin(svData[0], svData[5], filteredTads[0][1])
			# 		
			# 		svInteractionsBin = genomicBin[3].getElementsByRange(svData[1], filteredTads[0][2])
			# 		svGenesBin = genomicBin[3].getGenesByRange(svData[1], filteredTads[0][2])
			# 		
			# 		svInteractionsFirstTad = filteredTads[0][3].getElementsByRange(filteredTads[0][1], svData[5])
			# 		svGenesFirstTad = filteredTads[0][3].getGenesByRange(filteredTads[0][1], svData[5])
			# 		
			# 		for gene in svGenesBin:
			# 			#These genes gain interactions from the TAD.
			# 		#	print "adding gains from left TAD: ", len(svInteractionsFirstTad)
			# 			gene.addGainedEQTLs(svInteractionsFirstTad, svData[7])
			# 			
			# 		
			# 		for gene in svGenesFirstTad:
			# 			#These genes gain interactions from the SV part in the bin.
			# 		#	print "adding gains from bin: ", len(svInteractionsBin)
			# 			gene.addGainedEQTLs(svInteractionsBin, svData[7])
			# 			
			# 		
			# 		
			# 	else:	
			# 		#In this first case, we only have a new part after the original TAD, where genes can gain. there are no losses.
			# 		
			# 		#First part of DUP until TAD boundary is in an extra TAD. All genes in there get eQTLs in there.
			# 		#The second part of the DUP is duplicated into a bin, so all genes in the SV and in the bin gain eQTLs in the bin and in the SV. 
			# 		
			# 		newTad1Start = filteredTads[0][2]
			# 		#The end of the new TAD1 is the leftmost TAD end + leftmostTAD end - duplication start
			# 		newTad1End = filteredTads[0][2] + (filteredTads[0][2] - svData[1])
			# 		
			# 		#Get all interactions from the start of the SV until the end of the first TAD.
			# 		svInteractionsFirstTad = filteredTads[0][3].getElementsByRange(svData[1], filteredTads[0][2])
			# 		svGenesFirstTad = filteredTads[0][3].getGenesByRange(svData[1], filteredTads[0][2])
			# 		
			# 		#Get the bin and the interactions inside the SV part specifically
			# 		genomicBin = genome.collectGenomicBin(svData[0], filteredTads[0][2], svData[5])
			# 		
			# 		svInteractionsBin = genomicBin[3].getElementsByRange(svData[1], filteredTads[0][2])
			# 		svGenesBin = genomicBin[3].getGenesByRange(svData[1], filteredTads[0][2])
			# 		
			# 		#Get the genes in the bin outside of the SV
			# 		#Get the interactions in the bin outside of the SV
			# 		interactionsBin = genomicBin[3].getElementsByRange(svData[5], genomicBin[2])
			# 		genesBin = genomicBin[3].getGenesByRange(svData[5], genomicBin[2])
			# 		
			# 		for gene in svGenesFirstTad:
			# 			#Each gene in this bin gets all eQTLs that are within the SV.
			# 		#	print "adding gains from left TAD: ", len(svInteractionsFirstTad)
			# 			gene.addGainedEQTLs(svInteractionsFirstTad, svData[7])
			# 			
			# 		
			# 		for gene in svGenesBin:
			# 			#Each gene here gains eQTLs from outside of the SV in the bin.
			# 		#	print "adding gains from bin: ", len(interactionsBin)
			# 			gene.addGainedEQTLs(interactionsBin, svData[7])
			# 
			# 			
			# 		for gene in genesBin:
			# 			#Each gene here gains eQTLs from inside the SV.
			# 		#	print "adding gains from bin within SV: ", len(svInteractionsBin)
			# 			gene.addGainedEQTLs(svInteractionsBin, svData[7])
			# 			
					
			#print len(gene.gainedEQTLs[svData[7]])		
				
				
			

