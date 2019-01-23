"""
	The goal of this script is to take a set of SVs as input, and then for each of these determine what the derivate TADs are on the part of the genome that they disrupt.
	
	For duplications, inversions, and translocations. 


"""

import re
import numpy as np
from tad import TAD
from sv import SV
from gene import Gene
from eQTL import EQTL

class DerivativeTADMaker:
	
	
	def __init__(self, svData, genes, tadData, genome):
		
		# allBins = genome.getAllBinsInTADFormat()
		# #Combine the TADs with the bins
		# combinedData = np.concatenate((tadData, allBins))
		# #Should the combined data be sorted?
		# print combinedData.shape
		# print combinedData
		
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
				
		#For the translocations separately
		#1. For each SV, determine which TAD the SVs are in
		# tadsPerSV = self.matchTADsWithTranslocations(svData, tadData)
		# 
		# #2. Group the SVs that form a 'chain' and have at least 1 TAD in overlap
		# svGroups = self.defineGroupsOfTranslocations(tadsPerSV)
		# 
		# #3. Call the derivative TAD maker on this group of SVs and let it assign the gains/losses to the genes
		# self.determineDerivativeTADs([svGroups, tadsPerSV], tadData, genome, "trans")
		# 
		
		print "done making derivative TADs"
		
		return 0
	
	def defineGroupsOfTranslocations(self, tadsPerSV):
		"""
			Loop through the SVs and determine which SVs form a 'chain' and affect the same TADs
		"""
		
		### TO DO:
		# - Make sure that the SVs affecting only 1 position are also captured
		
		svGroups = dict()
		tadsPerSVKeys = tadsPerSV.keys()
		for svInd in range(0, len(tadsPerSV)):
			sv1 = tadsPerSVKeys[svInd]
			
			#Make sure that the breakpoint/sv is not taking palce entirely within one tad
			if tadsPerSV[sv1][0][0][3] == tadsPerSV[sv1][1][0][3]:
				continue
			
			currentGroup = []
			currentGroupReversed = []
			uniqueSVs = []
			for sv2Ind in range(svInd+1, len(tadsPerSV)):
				
				sv2 = tadsPerSVKeys[sv2Ind]
				
				if sv1.sampleName != sv2.sampleName:
					continue
				
				#If the 2nd SV starts in the same TAD as the 1st ends in, and these coordinates are withi 100 bp of each other, cluster them together
				#Here we then assume that the first SV in the list is actually the smallest SV on the smaller chromosme (the data is sorted so this should work fine)
				
				#Make sure that neither of the SVs starts and ends entirely within the same TAD.
				if tadsPerSV[sv2][0][0][3] == tadsPerSV[sv2][1][0][3]:
					continue
					
				
				firstTad = tadsPerSV[sv1][1][0][3]
				secondTad = tadsPerSV[sv2][0][0][3]
				if firstTad == secondTad:
					
					if sv1 not in uniqueSVs:
						#print "adding SV: ", sv1, sv1.chr1, sv1.s1, sv1.chr2, sv1.e2, sv1.sampleName
						currentGroup.append([sv1.s1, sv1])
						currentGroupReversed.append(sv1)
						uniqueSVs.append(sv1)
					if sv2 not in uniqueSVs:
						#print "adding SV: ", sv2, sv2.chr1, sv2.s1, sv2.chr2, sv2.e2, sv2.sampleName
						currentGroup.append([sv2.s1, sv2])
						currentGroupReversed.append(sv2)
						uniqueSVs.append(sv2)
					
			
			if len(currentGroup) > 0:
				#sort the current group by position
				currentGroup = np.array(currentGroup)
				
				currentGroupSorted = currentGroup[currentGroup[:,0].argsort()]
				
				groupList = []
				for sv in currentGroupSorted:
					groupList.append(sv[1])
				
				svGroups[sv1] = groupList
				
		# for group in svGroups:
		# 	print "new group: "
		# 	
		# 	for sv in svGroups[group]:
		# 		print sv.chr1, sv.s1, sv.chr2, sv.e2, sv.sampleName
		return svGroups

	def matchTADsWithTranslocations(self, svData, tadData):
		
		tadsPerSV = dict()
		
		for sv in svData:
			
			#Focus on translocations only
			typeMatch = re.search("chr", sv[8].svType, re.IGNORECASE)
			if typeMatch is None:
				continue

			#1. Check which TADs are affected by the breakpoint (there should be only 1 on each side of the SV)
			
			#Check if the SV is intrachromosomal or interchromosomal
			if sv[0] == sv[3]:
				tadChrSubset = tadData[tadData[:,0] == sv[0]]
				
				startMatches = (sv[1] > tadChrSubset[:,1]) * (sv[1] < tadChrSubset[:,2])
				matchingTadStart = tadChrSubset[startMatches]
				
				if len(matchingTadStart) < 1:
					continue #if there is no TAD affected, we skip it for now
			
				endMatches = (sv[5] > tadChrSubset[:,1]) * (sv[5] < tadChrSubset[:,2])
				matchingTadEnd = tadChrSubset[endMatches]
			
				if len(matchingTadEnd) < 1:
					continue
			
				tadsPerSV[sv[8]] = [matchingTadStart, matchingTadEnd]
		
		return tadsPerSV
		
		
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
		
		## translocation cases to test:
		#1. 1 translocation involving 2 TADs
		#2. 2 translocations involving 3 TADs
		#3. 3 translocations involving 3 TADs, where the 3rd translocation is in the same TAD as the 2nd
		#4. 
		
		#Try a very simple translocation that affects only 2 TADs.
		
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr1", 500, 700, TAD("chr1", 500, 700)]
		# 
		# sv1 = ["chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# svData = [svGroups, tadsPerSV]
		# 
		# #Assign test genes and eQTLs to the TADs
		# #A genes eQTLB. A loses eQTL A. 
		# geneA = Gene("A", "chr1", 100, 105)
		# geneB = Gene("B", "chr1", 670, 680)
		# eQTLA = EQTL("chr1", 111, 111)
		# eQTLB = EQTL("chr1", 610, 610)
		# tad1[3].genes = [geneA]
		# tad1[3].eQTLInteractions = [eQTLA]
		# tad2[3].genes = [geneB]
		# tad2[3].eQTLInteractions = [eQTLB]
		# geneA.eQTLs = [eQTLA]
		# geneB.eQTLs = [eQTLB]
		
		
		# 2 translocations involving 3 TADs
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr1", 500, 700, TAD("chr1", 500, 700)]
		# tad3 = ["chr1", 800, 900, TAD("chr1", 800, 900)]
		# 
		# sv1 = ["chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 650, 670, "chr1", 820, 830, "sample1", "cancerTypeA", SV("chr1", 650, 660, "chr1", 820, 830, "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8], sv2[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# tadsPerSV[sv2[8]] = [[tad2], [tad3]]
		# svData = [svGroups, tadsPerSV]
		# 
		# #TAD A and C have a gene before the SV.
		# #Tad B has 2 eQTLs on either side of the breakpoint. The first will be gained by gene A, the second by gene C. 
		# geneA = Gene("A", "chr1", 100, 105)
		# geneC = Gene("C", "chr1", 840, 845)
		# eQTLA = EQTL("chr1", 601, 601)
		# eQTLB = EQTL("chr1", 680, 680)
		# tad1[3].genes = [geneA]
		# tad3[3].genes = [geneC]
		# tad2[3].eQTLInteractions = [eQTLA, eQTLB]
		
		
		
		# # 2 translocations involving 3 TADs, where the 2nd translocation is in the same TAD as the 2nd
		#With the first SV, geneA gains an eQTL, but then loses it again after the second SV. 
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr1", 500, 700, TAD("chr1", 500, 700)]
		# tad3 = ["chr1", 800, 900, TAD("chr1", 800, 900)]
		# 
		# sv1 = ["chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 130, 140, "chr1", 820, 830, "sample1", "cancerTypeA", SV("chr1", 130, 140, "chr1", 820, 830, "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8], sv2[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# tadsPerSV[sv2[8]] = [[tad1], [tad3]]
		# svData = [svGroups, tadsPerSV]
		# 
		# geneA = Gene("A", "chr1", 100, 105)
		# eQTLB = EQTL("chr1", 610, 610)
		# tad1[3].genes = [geneA]
		# tad2[3].eQTLInteractions = [eQTLB]
		# 
		# # # 2 translocations involving 3 TADs, where the 2nd translocation is in the same TAD as the 2nd
		# #Here in the remaining part, the gene of TAD B can get into contact with the eQTL in TAD C
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr1", 500, 700, TAD("chr1", 500, 700)]
		# tad3 = ["chr1", 800, 900, TAD("chr1", 800, 900)]
		# 
		# sv1 = ["chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 130, 140, "chr1", 820, 830, "sample1", "cancerTypeA", SV("chr1", 130, 140, "chr1", 820, 830, "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8], sv2[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# tadsPerSV[sv2[8]] = [[tad1], [tad3]]
		# svData = [svGroups, tadsPerSV]
		# 
		# geneB = Gene("B", "chr1", 650, 660)
		# eQTLC = EQTL("chr1", 810, 810)
		# tad2[3].genes = [geneB]
		# tad3[3].eQTLInteractions = [eQTLC]
		
		# 
		# # 2 translocations involving 2 TADs, but the 2nd translocation ends before the first translocation
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr2", 500, 700, TAD("chr2", 500, 700)]
		# tad3 = ["chr1", 800, 900, TAD("chr1", 800, 900)]
		# 
		# sv1 = ["chr1", 110, 120, "chr2", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "chr2", 550, 600, "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 130, 140, "chr2", 510, 520, "sample1", "cancerTypeA", SV("chr1", 130, 140, "chr2", 510, 520, "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8], sv2[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# tadsPerSV[sv2[8]] = [[tad1], [tad2]]
		# svData = [svGroups, tadsPerSV]
		# 
		### TRANSLOCATIONS ###
		if svType == "trans":
			
			svGroups = svData[0]
			tadsPerSV = svData[1]
			
			#1. For each group of SVs, determine a derivative genome
			#1.1 Keep a list of elements that are gained/lost for each translocation individually
			#1.2 For the next translocation, update that list again. If there are elements that are now lost, remove them from the list
			
			updatedTadPos = dict() #keep the TADs and the new start/ends after 
			
			for group in svGroups:
				print "new group:"
				
				gains = dict()
				losses = dict()
				
				#First we need to keep track of all healthy TADs as well that the sv group affects
				#From this list we then need to later on at the first iteration determine which ones are affected
				#Then in the code itself we will copy the TADs if these are not affected and will remain as-is. 
				allInvolvedTads = []
				for sv in svGroups[group]:
					for tad in range(0, len(tadsPerSV[sv])):
						if tadsPerSV[sv][tad][0][3] not in allInvolvedTads:
							allInvolvedTads.append(tadsPerSV[sv][tad][0][3])
				filteredInvolvedTads = [] #after we filter out the first TADs that are affected, keep these here
				
				updatedTads = []
				
				for sv in svGroups[group]:
					print "SV: ", sv.chr1, sv.s1, sv.e1, sv.chr2, sv.s2, sv.e2, sv.sampleName
					
					
					
					#Failsafe to make sure that we never use the TADs from the previous iteration here
					rightTad = None
					leftTad = None
					
					#1. Get the affected TADs
					if len(updatedTads) < 1:
						leftTad = tadsPerSV[sv][0][0][3]
						rightTad = tadsPerSV[sv][1][0][3]
						
						#Every TAD that is not affected should be removed from the involved TADs. This is to make sure that we later on 
						for tadInd in range(0, len(allInvolvedTads)):
							tad = allInvolvedTads[tadInd]
							if tad != leftTad or tad != rightTad:
								filteredInvolvedTads.append([tad.chromosome, tad.start, tad.end, tad])
						
					else: #get the affected TADs from the updated TADs.
						#print "SV pos is: ", sv.s1, sv.e2
						for tad in updatedTads: ###Here this is going wrong because the original TAD is never here, and so the right TAD that we use is from the previous iteratio. 
							
							if sv.s1 < tad[2] and sv.s1 > tad[1]:
								leftTad = tad[3]
							if sv.e2 < tad[2] and sv.e2 > tad[1]:
								rightTad = tad[3]
						#This is the extra check to make sure that we also include TADs that are not affected in the first iteration
						for tad in filteredInvolvedTads:
							if sv.s1 < tad[2] and sv.s1 > tad[1]:
								leftTad = tad[3]
							if sv.e2 < tad[2] and sv.e2 > tad[1]:
								rightTad = tad[3]
					
						
					
					#2. Get the elements on the left of the breakpoint and the right of the breakpoint
					leftSideElements = leftTad.getElementsByRange(leftTad.start, sv.s1)
					leftSideGenes = leftTad.getGenesByRange(leftTad.start, sv.s1)
					
					rightSideElements = rightTad.getElementsByRange(sv.e2, rightTad.end)
					rightSideGenes = rightTad.getGenesByRange(sv.e2, rightTad.end)
					
					#Also get the elements in the remaining part
					remainingElementsLeft = leftTad.getElementsByRange(sv.s1, leftTad.end)
					remainingGenesLeft = leftTad.getGenesByRange(sv.s1, leftTad.end)
					
					remainingElementsRight = rightTad.getElementsByRange(rightTad.start, sv.e2)
					remainingGenesRight = rightTad.getGenesByRange(rightTad.start, sv.e2)
					
					#3. Make derivative TADs from the SV and add elements, keeping the reference positions.
					
					#The new TADs will be the left part together with the right part, and everything between that remains a separate part (to be determined by another SV)
					newTad = [leftTad.chromosome, leftTad.start, rightTad.end, TAD(leftTad.chromosome, leftTad.start, rightTad.end)]
					
					
					#add the elements to the TAD
					newTad[3].setEQTLInteractions(leftSideElements + rightSideElements)
					newTad[3].setGenes(leftSideGenes + rightSideGenes)
					
					#Also set the elements to the remaining parts
					remainingPart = [leftTad.chromosome, sv.s1, sv.e2, TAD(leftTad.chromosome, sv.s1, sv.e2)]
					
					remainingPart[3].setEQTLInteractions(remainingElementsLeft + remainingElementsRight)
					remainingPart[3].setGenes(remainingGenesLeft + remainingGenesRight)
					
					#For interchromosomal translocations, this would mean that we have 2 open ends, or intrachromosomal translocations that does not happen.
					#How about pieces that are moved to another location? How are these encoded? E.g. s1 - e1 moves to s2 - e2?
					
					#4. Store the new TADs.
					
					#Copy all TADs but the left and right affected one to a new set of TADs.
					updatedTads.append(newTad)
					updatedTads.append(remainingPart)
					for tad in tadsPerSV[sv]:
						if tad[0][3] != leftTad and tad[0][3] != rightTad:
							updatedTads.append(tad[0])

					#5. For the next SV, use the updated TADs to get the elements/positions and continue until the last SV
					
				#6. Eventually, check for the genes in the original TADs which elements thesd has, and compare it to the derivative. Use the delta informatin to collect gains and losses
				
				#For each of the TADs, see which genes are in them and are affected, and determine their eQTL interactions.
				
					print "new SV:"
				
					for tad in updatedTads:
						print "tad: ", tad
						for gene in tad[3].genes:
							print gene.name
							eQTLsInTad = []
							for eQTL in gene.eQTLs:
								if eQTL.start > tad[3].start and eQTL.start < tad[3].end:
									eQTLsInTad.append([eQTL.chromosome, eQTL.start, eQTL])
							print "original eQTLs:"
							print eQTLsInTad
							print "new eQTLs:"
							for eQTL in tad[3].eQTLInteractions:
								print eQTL.chromosome, eQTL.start, eQTL
							#print tad[3].eQTLInteractions
				exit()
						
				
				#Compare this to the original case and make the derivative. 
					
					
			# 		#Get all elements and genes on the left side of the translocation
			# 		#Get the TAD that the sv is in on the left side
			# 		leftTad = tadsPerSV[sv][0][0][3]
			# 		if leftTad not in updatedTadPos:
			# 			print "Left side gains right side"
			# 			leftSideElements = leftTad.getElementsByRange(leftTad.start, sv.s1)
			# 			leftSideGenes = leftTad.getGenesByRange(leftTad.start, sv.s1)
			# 			updatedTadPos[leftTad] = [[leftTad.start, sv.s1],[sv.s1, leftTad.end]]
			# 			
			# 			
			# 			
			# 			#Determine the elements that the left side loses (which is everything in the TAD after the SV breakpoint), reverse for the lost side
			# 			lostElements = leftTad.getElementsByRange(sv.s1, leftTad.end)
			# 			lostGenes = leftTad.getGenesByRange(sv.s1, leftTad.end)
			# 			for gene in leftSideGenes:
			# 				if gene not in losses:
			# 					losses[gene] = dict()
			# 				for element in lostElements:
			# 					losses[gene][element] = sv.sampleName
			# 			for gene in lostGenes:
			# 				if gene not in losses:
			# 					losses[gene] = dict()
			# 				for element in leftSideElements:
			# 					losses[gene][element] = sv.sampleName
			# 
			# 		else: #if the TAD was previously affected by an SV, we should use the start position of the SV, which is the remaining part of the TAD. 
			# 			print "Left side not intact"
			# 			#In these specific cases, if the left side combines with the right side, we need to make sure that the right side also gains the updated part, not the entire left side again. 
			# 			
			# 			
			# 			#Check if the SV takes place in the remaining left part or in the right part.
			# 			if sv.s1 < updatedTadPos[leftTad][0][1] and sv.s1 > updatedTadPos[leftTad][0][0]:
			# 				print "Left side of previously affected left side gains right side"
			# 				#if the SV is in the left part of the TAD, update the positions in the left TAD following this SV.
			# 				
			# 				#Get the elements using the positions in the already updated TAD
			# 				leftSideElements = leftTad.getElementsByRange(updatedTadPos[leftTad][0][0], sv.s1)
			# 				leftSideGenes = leftTad.getGenesByRange(updatedTadPos[leftTad][0][0], sv.s1)
			# 				
			# 				#In this case, the updated positions should be: the previous start, with the new SV end. The right side will not get updated
			# 				updatedTadPos[leftTad] = [[updatedTadPos[leftTad][0][0], sv.s1],[updatedTadPos[leftTad][1][0], updatedTadPos[leftTad][1][1]]]
			# 				
			# 				#Determine the elements that the left side loses (which is everything in the TAD after the SV breakpoint), reverse for the lost side
			# 				
			# 				lostElements = leftTad.getElementsByRange(sv.s1, updatedTadPos[leftTad][0][1])
			# 				lostGenes = leftTad.getGenesByRange(sv.s1, updatedTadPos[leftTad][0][1])
			# 				for gene in leftSideGenes:
			# 					if gene not in losses:
			# 						losses[gene] = dict()
			# 					for element in lostElements:
			# 						losses[gene][element] = sv.sampleName
			# 				
			# 
			# 				for gene in lostGenes:
			# 					if gene not in losses:
			# 						losses[gene] = dict()
			# 					for element in leftSideElements:
			# 						losses[gene][element] = sv.sampleName
			# 
			# 		
			# 			if sv.s1 < updatedTadPos[leftTad][1][1] and sv.s1 > updatedTadPos[leftTad][1][0]:
			# 				#if the SV is in the left part of the TAD, update the positions in the left TAD following this SV.
			# 				print "Right side of previously affected left side gains right side"
			# 				
			# 				#Get the eleemnts using the positions in the already updated TAD
			# 				leftSideElements = leftTad.getElementsByRange(sv.s1, updatedTadPos[leftTad][1][1])
			# 				leftSideGenes = leftTad.getGenesByRange(sv.s1, updatedTadPos[leftTad][1][1])
			# 				
			# 				#Here update the right side
			# 				updatedTadPos[leftTad] = [[updatedTadPos[leftTad][0][0], updatedTadPos[leftTad][0][1]],[sv.s1, updatedTadPos[leftTad][1][1]]]
			# 		
			# 				lostElements = leftTad.getElementsByRange(updatedTadPos[leftTad][1][0], sv.s1)
			# 				lostGenes = leftTad.getGenesByRange(updatedTadPos[leftTad][1][0], sv.s1)
			# 				for gene in leftSideGenes:
			# 					if gene not in losses:
			# 						losses[gene] = dict()
			# 					for element in lostElements:
			# 						losses[gene][element] = sv.sampleName
			# 				
			# 
			# 				for gene in lostGenes:
			# 					if gene not in losses:
			# 						losses[gene] = dict()
			# 					for element in leftSideElements:
			# 						losses[gene][element] = sv.sampleName
			# 		
			# 		
			# 		#Get all elements and genes on the right side of the translocation
			# 		rightTad = tadsPerSV[sv][1][0][3]
			# 		print rightTad
			# 		if rightTad not in updatedTadPos:
			# 			print "Right side gains left side"
			# 			rightSideElements = rightTad.getElementsByRange(sv.e2, rightTad.end)
			# 			rightSideGenes = rightTad.getGenesByRange(sv.e2, rightTad.end)
			# 			
			# 			print "ements"
			# 			print rightSideElements
			# 			print rightSideGenes
			# 			
			# 			updatedTadPos[rightTad] = [[rightTad.start, sv.e2], [sv.e2, rightTad.end]]
			# 			
			# 			#Determine the elements that the left side loses (which is everything in the TAD after the SV breakpoint), reverse for the lost side
			# 			lostElements = rightTad.getElementsByRange(rightTad.start, sv.e2)
			# 			lostGenes = rightTad.getGenesByRange(rightTad.end, sv.e2)
			# 			for gene in rightSideGenes:
			# 				if gene not in losses:
			# 					losses[gene] = dict()
			# 				for element in lostElements:
			# 					losses[gene][element] = sv.sampleName
			# 			for gene in lostGenes:
			# 				if gene not in losses:
			# 					losses[gene] = dict()
			# 				for element in rightSideElements:
			# 					losses[gene][element] = sv.sampleName
			# 			
			# 		else: #if the TAD was previously affected by an SV, we should use the start position of the SV, which is the remaining part of the TAD. 
			# 			
			# 			#Check if the SV takes place in the remaining left part or in the right part.
			# 			if sv.e2 < updatedTadPos[rightTad][0][1] and sv.e2 > updatedTadPos[rightTad][0][0]:
			# 				print "Left side of previously affected right side gains left side"
			# 				#if the SV is in the left part of the TAD, update the positions in the left TAD following this SV.
			# 				
			# 				#Get the eleemnts using the positions in the already updated TAD
			# 				rightSideElements = rightTad.getElementsByRange(updatedTadPos[rightTad][0][0], sv.e2)
			# 				rightSideGenes = rightTad.getGenesByRange(updatedTadPos[rightTad][0][0], sv.e2)
			# 				
			# 				#Update the left side
			# 				updatedTadPos[rightTad] = [[updatedTadPos[rightTad][0][0], sv.e2],[updatedTadPos[rightTad][1][0], updatedTadPos[rightTad][1][1]]]
			# 				
			# 				lostElements = rightTad.getElementsByRange(sv.e2, updatedTadPos[rightTad][0][1])
			# 				lostGenes = rightTad.getGenesByRange(sv.e2, updatedTadPos[rightTad][0][1])
			# 				for gene in rightSideGenes:
			# 					if gene not in losses:
			# 						losses[gene] = dict()
			# 					for element in lostElements:
			# 						losses[gene][element] = sv.sampleName
			# 				
			# 
			# 				for gene in lostGenes:
			# 					if gene not in losses:
			# 						losses[gene] = dict()
			# 					for element in rightSideElements:
			# 						losses[gene][element] = sv.sampleName
			# 				
			# 		
			# 			if sv.e2 < updatedTadPos[rightTad][1][1] and sv.e2 > updatedTadPos[rightTad][1][0]:
			# 				#if the SV is in the left part of the TAD, update the positions in the left TAD following this SV.
			# 				print "Right side of previously affected right side gains left side"
			# 				#Get the eleemnts using the positions in the already updated TAD
			# 				rightSideElements = rightTad.getElementsByRange(sv.e2, updatedTadPos[rightTad][1][0])
			# 				rightSideGenes = rightTad.getGenesByRange(sv.e2, updatedTadPos[rightTad][1][0])
			# 		
			# 				#Update the right side
			# 				updatedTadPos[rightTad] = [[updatedTadPos[rightTad][0][0], updatedTadPos[rightTad][0][1]],[sv.e2, updatedTadPos[rightTad][1][1]]]
			# 				
			# 				lostElements = rightTad.getElementsByRange(updatedTadPos[rightTad][1][0], sv.e2)
			# 				lostGenes = rightTad.getGenesByRange(updatedTadPos[rightTad][1][0], sv.e2)
			# 				for gene in rightSideGenes:
			# 					if gene not in losses:
			# 						losses[gene] = dict()
			# 					for element in lostElements:
			# 						losses[gene][element] = sv.sampleName
			# 				
			# 
			# 				for gene in lostGenes:
			# 					if gene not in losses:
			# 						losses[gene] = dict()
			# 					for element in rightSideElements:
			# 						losses[gene][element] = sv.sampleName
			# 		
			# 		
			# 		### Needs a better check for gains/losses that are updated in the next iteration
			# 		
			# 		#updatedTadPos[leftTad] = [[leftTad.start, sv.s1],[sv.s1, leftTad.end]]
			# 		#updatedTadPos[rightTad] = [[rightTad.start, sv.e2], [sv.e2, rightTad.end]]
			# 		
			# 		
			# 		#Now every gene on either side gains the elements from the other side.
			# 		for gene in rightSideGenes:
			# 			if gene not in gains:
			# 			 	gains[gene] = dict()
			# 			for element in leftSideElements:
			# 				
			# 				gains[gene][element] = sv.sampleName
			# 			#gene.addGainedEQTLs(leftSideElements, sv.sampleName)
			# 		
			# 		for gene in leftSideGenes:
			# 			#gene.addGainedEQTLs(rightSideElements, sv.sampleName)
			# 			if gene not in gains:
			# 			 	gains[gene] = dict()
			# 			for element in rightSideElements:
			# 				
			# 				gains[gene][element] = sv.sampleName
			# 		
			# 		#Keep the updated positions in the TADs after each SV
			# 		#Use these positions to determine the gains/losses for each new SV
			# 
			# print 'Gains:'
			# for gene in gains:
			# 	print "gene: ", gene.name
			# 	for element in gains[gene]:
			# 		print "element: ", element.start
			# 
			# print 'Losses:'
			# for gene in losses:
			# 	print "gene: ", gene.name
			# 	for element in losses[gene]:
			# 		print "element: ", element.start

			
			
		
		
		
		
			
		### INVERSION ###
		if svType == "inv":
			
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
			# elif len(leftMostTad) < 1 or len(rightMostTad) < 1:
			# 	print "SV ends in one TAD"
			# 	
			# 	
			# 	
			# 	#Either the left side of the inversion or the right side is within a TAD, but the other part is within a genomic bin. 
			# 	#Here make use of genomic bins instead of TADs.
			# 	
			# 	#1. First get the genomic bin for the non-TAD end
			# 	#Case where the inversion ends in a TAD, but the left side is a genomic bin. 
			# 	if len(leftMostTad) < 1:
			# 		rightMostTad = rightMostTad[0]
			# 		genomicBin = genome.collectGenomicBin(svData[0], svData[1], svData[2])
			# 		
			# 		if genomicBin == None:
			# 			return
			# 		
			# 		#Collect the elements and genes that are gained and lost within the TAD or genomic bin
			# 		
			# 		leftSideElements = genomicBin[3].getElementsByRange(svData[1], genomicBin[2]) #From the start of the inversion until the end of the left most TAD
			# 		unaffectedElementsLeft = genomicBin[3].getElementsByRange(genomicBin[1], svData[1])
			# 
			# 		#Also get the genes
			# 		leftSideGenes = genomicBin[3].getGenesByRange(svData[1], genomicBin[2]) #From the start of the inversion until the end of the left most TAD
			# 		unaffectedGenesLeft = genomicBin[3].getGenesByRange(genomicBin[1], svData[1])
			# 		
			# 		#3. Collect all elements from the left TAD boundary until the end of the inversion.
			# 		
			# 		rightSideElements = rightMostTad[3].getElementsByRange(rightMostTad[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
			# 		unaffectedElementsRight = rightMostTad[3].getElementsByRange(svData[5], rightMostTad[2])
			# 		
			# 		rightSideGenes = rightMostTad[3].getGenesByRange(rightMostTad[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
			# 		unaffectedGenesRight = rightMostTad[3].getGenesByRange(svData[5], rightMostTad[2])
			# 	#Case where the inversion starts in a TAD, but the right side is in a genomic bin.
			# 	elif len(rightMostTad) < 1:
			# 		leftMostTad = leftMostTad[0]
			# 		genomicBin = genome.collectGenomicBin(svData[0], svData[4], svData[5])
			# 		if genomicBin == None:
			# 			return
			# 		
			# 		#2. Collect all elements until the right TAD boundary inside the inversion.
			# 		
			# 		leftSideElements = leftMostTad[3].getElementsByRange(svData[1], leftMostTad[2]) #From the start of the inversion until the end of the left most TAD
			# 		unaffectedElementsLeft = leftMostTad[3].getElementsByRange(leftMostTad[1], svData[1])
			# 
			# 		#Also get the genes
			# 		leftSideGenes = leftMostTad[3].getGenesByRange(svData[1], leftMostTad[2]) #From the start of the inversion until the end of the left most TAD
			# 		unaffectedGenesLeft = leftMostTad[3].getGenesByRange(leftMostTad[1], svData[1])
			# 		
			# 		#3. Collect all elements from the left TAD boundary until the end of the inversion.
			# 		
			# 		rightSideElements = genomicBin[3].getElementsByRange(genomicBin[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
			# 		unaffectedElementsRight = genomicBin[3].getElementsByRange(svData[5], genomicBin[2])
			# 		
			# 		rightSideGenes = genomicBin[3].getGenesByRange(genomicBin[1], svData[5]) #from the start of the rightmost TAD until the end of the inversion
			# 		unaffectedGenesRight = genomicBin[3].getGenesByRange(svData[5], genomicBin[2])
			# 			
			
			
			#Assigning the gains and losses to the genes is independent of the type of inversion
			#print "Copying genes and elements after SV"		
			#All genes that were originally in the left TAD (outisde of the inversion) will gain elements of the right side of the inversion
			#All unaffected genes on the left will lose the eQTLs that are in the left side of the inversion
			for gene in unaffectedGenesLeft:
				
				gene.addGainedEQTLs(rightSideElements, svData[7])
				gene.addLostEQTLs(leftSideElements, svData[7])
				#print "Number of gained right side elements: ", len(rightSideElements), " for genes ", len(unaffectedGenesLeft)
				
				#if svData[7] in gene.gainedEQTLs:
				#print "r: ", len(gene.gainedEQTLs[svData[7]])
				#print "l: ", len(gene.lostEQTLs[svData[7]])
			
			#All genes in the right side of the inversion will gain elements from the original left TAD.
			#All genes in the right side will lose interactions with eQTLs in the unaffected right TAD. 
			for gene in rightSideGenes:
				
				gene.addGainedEQTLs(unaffectedElementsLeft, svData[7])
				#print "Number of unaffected elements right: ", len(unaffectedElementsRight), " for genes ", len(rightSideGenes)
				gene.addLostEQTLs(unaffectedElementsRight, svData[7])
			
			#vice versa but then for the right TAD and right side of the inversion.
			#The lost eQTLs are the ones that are in the right side of the inversion
			for gene in unaffectedGenesRight:
				
				gene.addGainedEQTLs(leftSideElements, svData[7])
				#print "Number of gained right side elements 2: ", len(rightSideElements), " for genes ", len(unaffectedGenesRight)
				gene.addLostEQTLs(rightSideElements, svData[7])
			
			#The lost eQTLs are the ones that are in the unaffected original left TAD
			for gene in leftSideGenes:
				
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
						if gene.name == "ARHGEF10L":
							print svData
						#	exit()
							
						gene.addGainedEQTLs(svInteractionsLastTad, svData[7])
					print "(2) Number of genes to add gains: ", len(svGenesLastTad)
					for gene in svGenesLastTad:
						#print "adding gains from left TAD: ", len(svInteractionsFirstTad)
						if gene.name == "ARHGEF10L":
							print svData
							#exit()
						gene.addGainedEQTLs(svInteractionsFirstTad, svData[7])
					
					#The last TAD remains the same overall.
					#Only the TADs in the middle are duplicated.
					
					#Each gene gains the eQTLs that are within the TAD that the gene is located in. 
					
					for tad in followingTads:
						for gene in tad[3].genes:
							if gene.name == "ARHGEF10L":
								print svData
								#exit()
							
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
						if gene.name == "ARHGEF10L":
							print svData
							#exit()
						#Each gene in this bin gets all eQTLs that are within the SV.
						gene.addGainedEQTLs(svInteractionsSecondTad, svData[7])
						
					for gene in svGenesSecondTad:
						if gene.name == "ARHGEF10L":
							print svData
							#exit()
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
				
				
			

