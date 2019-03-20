"""
	The goal of this class is to take a set of SVs as input, and then for each of these determine what the derivate TADs are on the part of the genome that they disrupt.
	In the derivative TADs, we re-assign gains and losses of genomic elements to the genes that are affected when the derivative TADs are made. 
	
	For duplications, inversions, and (intra & interchromosomal) translocations.
	
	TO DO:
	- Make sure that the effects of SVs are not counted more than once. Currently, since translocations and e.g. inversions are considered separately, we do not consider effects caused by e.g. first a translocation, and then an inversion
	that will un-do the gains/losses of the translocation.
	- finish documentation


"""

import re
import numpy as np
from tad import TAD
from sv import SV
from gene import Gene
from element import Element

class DerivativeTADMaker:
	
	
	def __init__(self, svData, genes, tadData):
		
		self.linkSVEffectsToGenes(svData, genes, tadData)
		
	
	def linkSVEffectsToGenes(self, svData, genes, tadData):
		
		"""
			For every SV, determine the type. If this is an inversion or duplication, we only need this particular SV.
			If the SV is a translocation, we collect a set of SVs that 'work together'. Based on a window around their start/end positions.
			For the collected SV or SVs, we first make derivative TADs. 
			Then after these derivative TADs have been made, we go through the genes that are present within these new TADs, and add the newly gained or remove newly lost genomic elements for these genes.
			
			genes: (numpy array) array with the genes and their information. chr, start, end, geneObject
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject.
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
			
			
		"""
		print "Linking SV effects to genes"
		invCount = 0
		dupCount = 0
		delCount = 0
		
		#Inversions
		for sv in svData:
			typeMatch = re.search("inv", sv[8].svType, re.IGNORECASE)
			if typeMatch is not None:
				self.determineDerivativeTADs(sv, tadData, "inv")
				invCount += 1
				print "inversion count: ", invCount
		
		
		#Duplications
		for sv in svData:
			typeMatch = re.search("dup", sv[8].svType, re.IGNORECASE)
			if typeMatch is not None:
				
				self.determineDerivativeTADs(sv, tadData, "dup")
				dupCount += 1
				print "duplication count: ", dupCount
				
		#Deletions		
		for sv in svData:
			typeMatch = re.search("del", sv[8].svType, re.IGNORECASE)
			if typeMatch is not None:
				
				self.determineDerivativeTADs(sv, tadData, "del")
				delCount += 1
				print "deletion count: ", delCount		
				
				
		#For the translocations separately
		# 1. For each SV, determine which TAD the SVs are in
		print "matching TADs with translocations"
		tadsPerSV = self.matchTADsWithTranslocations(svData, tadData)
		
		#2. Group the SVs that form a 'chain' and have at least 1 TAD in overlap
		print "making SV groups"
		svGroups = self.defineGroupsOfTranslocations(tadsPerSV)
		
		print "determine derivative TADs"
		#3. Call the derivative TAD maker on this group of SVs and let it assign the gains/losses to the genes
		import time
		startTime = time.time()
		self.determineDerivativeTADs([svGroups, tadsPerSV], tadData, "trans")
		endTime = time.time()
		print "Took ", endTime - startTime, " to determine the derivative TADs"
		
		print "done making derivative TADs"
	
	def defineGroupsOfTranslocations(self, tadsPerSV):
		"""
			Loop through the SVs and determine which SVs form a 'chain' and affect the same TADs.
			For every SV, we first get the TAD that it ends in and that the next SV in the group starts in. If this is true, we group them together on the 'chain'. 
			
			TO DO:
			- Make sure that the SVs affecting only 1 position are also captured
			
			tadsPerSV: (dictionary) the keys of this dictionary are an SV object, the values are lists with the TADs that are affected by this SV. Output from matchTADsWithTranslocations()
			
			return
			svGroups: (dictionary) the keys are the first SV object that is part of the chain of SVs. The values are the SV objects that are part of this group. 
		"""

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
				
				#If the 2nd SV starts in the same TAD as the 1st ends in, cluster them together
				#Here we then assume that the first SV in the list is actually the smallest SV on the smaller chromosme (the data is sorted so this should work fine)
				
				#Make sure that neither of the SVs starts and ends entirely within the same TAD.
				if tadsPerSV[sv2][0][0][3] == tadsPerSV[sv2][1][0][3]:
					continue
					
				
				firstTad = tadsPerSV[sv1][1][0][3]
				secondTad = tadsPerSV[sv2][0][0][3]
				if firstTad == secondTad:
					
					if sv1 not in uniqueSVs:
						currentGroup.append([sv1.s1, sv1])
						currentGroupReversed.append(sv1)
						uniqueSVs.append(sv1)
					if sv2 not in uniqueSVs:
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
				
		return svGroups

	def matchTADsWithTranslocations(self, svData, tadData):
		"""
			For every SV, find the TAD that it ends in on the left and on the right. This can be the same TAD.
			
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject. Can be empty in SNV mode.
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
			
			return:
			tadsPerSV: (dictionary) dictionary with the SV object as the key, and the left TAD (in np array format, like tadData) as first array element,
									right TAD as the second. 
			
		"""
		tadsPerSV = dict()
		
		for sv in svData:
			
			#Focus on translocations only
			interChrTypeMatch = re.search("chr", sv[8].svType, re.IGNORECASE)
			transTypeMatch = re.search("trans", sv[8].svType, re.IGNORECASE)
			rangeTypeMatch = re.search("range", sv[8].svType, re.IGNORECASE)
			if interChrTypeMatch is None and transTypeMatch is None and rangeTypeMatch is None:
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
			else: #interchromosomal SV
				tadChr1Subset = tadData[tadData[:,0] == sv[0]]
				
				startMatches = (sv[1] > tadChr1Subset[:,1]) * (sv[1] < tadChr1Subset[:,2])
				matchingTadStart = tadChr1Subset[startMatches]
				
				if len(matchingTadStart) < 1:
					continue #if there is no TAD affected, we skip it for now
				
				tadChr2Subset = tadData[tadData[:,0] == sv[3]]
			
				endMatches = (sv[5] > tadChr2Subset[:,1]) * (sv[5] < tadChr2Subset[:,2])
				matchingTadEnd = tadChr2Subset[endMatches]
			
				if len(matchingTadEnd) < 1:
					continue
			
				tadsPerSV[sv[8]] = [matchingTadStart, matchingTadEnd]
	
		return tadsPerSV
		
		
	def determineDerivativeTADs(self, svData, tadData, svType):	
	
		"""
			Given an SV or a set of SVs, depending on the type of SVs, we compute how the affected region of the genome will look after the SV.
			We then make a set of new TAD objects that are located next to each other, and update all elements that are now inside these new/affected TADs. 
		
			svData: (numpy array) array with the SVs and their information. chr1, s1, e1, chr2, s2, e2, cancerType, sampleName, svObject. Can be empty in SNV mode.
			tadData: (numpy array) array with the TADs and their information. chr, start, end, tadObject
			svType: (string) type of SV that we should determine the derivative TAD for. Either del, inv, dup or trans.
			
		
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
		# sv1 = ["chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "+", "chr1", 550, 600, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "-", "chr1", 550, 600, "+", "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8], sv2[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# tadsPerSV[sv2[8]] = [[tad1], [tad2]]
		# svData = [svGroups, tadsPerSV]
		# 
		# #Assign test genes and eQTLs to the TADs
		# #A gains eQTLB. A loses eQTL A. 
		# geneA = Gene("A", "chr1", 100, 105)
		# geneB = Gene("B", "chr1", 580, 590)
		# eQTLA = EQTL("chr1", 128, 128)
		# eQTLB = EQTL("chr1", 591, 592)
		# tad1[3].genes = [geneA]
		# tad1[3].eQTLInteractions = [eQTLA]
		# tad2[3].genes = [geneB]
		# tad2[3].eQTLInteractions = [eQTLB]
		# geneA.eQTLs = [eQTLA]
		# geneB.eQTLs = [eQTLB]
		
		
		#Interchromosomal case
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr2", 500, 700, TAD("chr2", 500, 700)]
		# 
		# sv1 = ["chr1", 110, 120, "chr2", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "+", "chr2", 550, 600, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 110, 120, "chr2", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "-", "chr2", 550, 600, "+", "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8], sv2[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# tadsPerSV[sv2[8]] = [[tad1], [tad2]]
		# svData = [svGroups, tadsPerSV]
		# 
		# #Assign test genes and eQTLs to the TADs
		# #A gains eQTLB. A loses eQTL A. 
		# geneA = Gene("A", "chr1", 100, 105)
		# geneB = Gene("B", "chr2", 580, 590)
		# eQTLA = EQTL("chr1", 128, 128)
		# eQTLB = EQTL("chr2", 591, 592)
		# tad1[3].genes = [geneA]
		# tad1[3].eQTLInteractions = [eQTLA]
		# tad2[3].genes = [geneB]
		# tad2[3].eQTLInteractions = [eQTLB]
		# geneA.eQTLs = [eQTLA]
		# geneB.eQTLs = [eQTLB]
		
		#Case where one end remains open. B gains A, A loses A. 
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr2", 500, 700, TAD("chr2", 500, 700)]
		# 
		# sv1 = ["chr1", 110, 120, "chr2", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "-", "chr2", 550, 600, "+", "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# svData = [svGroups, tadsPerSV]
		# 
		# #Assign test genes and eQTLs to the TADs
		# #A gains eQTLB. A loses eQTL A. 
		# geneA = Gene("A", "chr1", 100, 105)
		# geneB = Gene("B", "chr2", 580, 590)
		# eQTLA = EQTL("chr1", 128, 128)
		# eQTLB = EQTL("chr2", 591, 592)
		# tad1[3].genes = [geneA]
		# tad1[3].eQTLInteractions = [eQTLA]
		# tad2[3].genes = [geneB]
		# tad2[3].eQTLInteractions = [eQTLB]
		# geneA.eQTLs = [eQTLA]
		# geneB.eQTLs = [eQTLB]
		# 
		# Case where the DNA is inverted when re-joining the translocations. Gene A gains eQTL B and loses eQTL A. Gene B does not lose eQTL B but does not gain anything. 
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr2", 500, 700, TAD("chr2", 500, 700)]
		# 
		# sv1 = ["chr1", 110, 120, "chr2", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "+", "chr2", 550, 600, "+", "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 110, 120, "chr2", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "-", "chr2", 550, 600, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8], sv2[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# tadsPerSV[sv2[8]] = [[tad1], [tad2]]
		# svData = [svGroups, tadsPerSV]
		# 
		# #Assign test genes and eQTLs to the TADs
		# #A gains eQTLB. A loses eQTL A. 
		# geneA = Gene("A", "chr1", 100, 105)
		# geneB = Gene("B", "chr2", 580, 590)
		# eQTLA = EQTL("chr1", 128, 128)
		# eQTLB = EQTL("chr2", 591, 592)
		# tad1[3].genes = [geneA]
		# tad1[3].eQTLInteractions = [eQTLA]
		# tad2[3].genes = [geneB]
		# tad2[3].eQTLInteractions = [eQTLB]
		# geneA.eQTLs = [eQTLA]
		# geneB.eQTLs = [eQTLB]
		
		#Case where only the right 2 parts are joined and inverted, the other ends remain open
		#Nothing happens in eQTL gains/losses
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr2", 500, 700, TAD("chr2", 500, 700)]
		# 
		# sv1 = ["chr1", 110, 120, "chr2", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "-", "chr2", 550, 600, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# 
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# svData = [svGroups, tadsPerSV]
		# 
		# #Assign test genes and eQTLs to the TADs
		# #A gains eQTLB. A loses eQTL A. 
		# geneA = Gene("A", "chr1", 100, 105)
		# geneB = Gene("B", "chr2", 580, 590)
		# eQTLA = EQTL("chr1", 128, 128)
		# eQTLB = EQTL("chr2", 591, 592)
		# tad1[3].genes = [geneA]
		# tad1[3].eQTLInteractions = [eQTLA]
		# tad2[3].genes = [geneB]
		# tad2[3].eQTLInteractions = [eQTLB]
		# geneA.eQTLs = [eQTLA]
		# geneB.eQTLs = [eQTLB]
		# 
		
		
		
		# 2 translocations involving 3 TADs
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr1", 500, 700, TAD("chr1", 500, 700)]
		# tad3 = ["chr1", 800, 900, TAD("chr1", 800, 900)]
		# 
		# sv1 = ["chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "+", "chr1", 550, 600, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 650, 670, "chr1", 820, 830, "sample1", "cancerTypeA", SV("chr1", 650, 660, "+", "chr1", 820, 830, "-", "sample1", "cancerTypeA", "intraChromosomal")]
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
		# 
		
		
		# # 2 translocations involving 3 TADs, where the 2nd translocation is in the same TAD as the 2nd
		# A loses eQTL A and gains eQTL B
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr1", 500, 700, TAD("chr1", 500, 700)]
		# tad3 = ["chr1", 800, 900, TAD("chr1", 800, 900)]
		# 
		# sv1 = ["chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "+", "chr1", 550, 600, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 130, 140, "chr1", 820, 830, "sample1", "cancerTypeA", SV("chr1", 130, 140, "+", "chr1", 820, 830, "-", "sample1", "cancerTypeA", "intraChromosomal")]
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
		
		
		# # # 2 translocations involving 3 TADs, where the 2nd translocation is in the same TAD as the 2nd
		#Here in the remaining part, the gene of TAD B can get into contact with the eQTL in TAD C
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr1", 500, 700, TAD("chr1", 500, 700)]
		# tad3 = ["chr1", 800, 900, TAD("chr1", 800, 900)]
		# 
		# sv1 = ["chr1", 110, 120, "chr1", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "+", "chr1", 550, 600, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 640, 640, "chr1", 820, 830, "sample1", "cancerTypeA", SV("chr1", 640, 640, "-", "chr1", 820, 830, "+", "sample1", "cancerTypeA", "intraChromosomal")]
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
		
		
		# 2 translocations involving 2 TADs, but the 2nd translocation ends before the first translocation
		# tad1 = ["chr1", 100, 300, TAD("chr1", 100, 300)]
		# tad2 = ["chr2", 500, 700, TAD("chr2", 500, 700)]
		# 
		# sv1 = ["chr1", 110, 120, "chr2", 550, 600, "sample1", "cancerTypeA", SV("chr1", 110, 120, "+", "chr2", 550, 600, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr1", 130, 140, "chr2", 510, 520, "sample1", "cancerTypeA", SV("chr1", 130, 140, "+", "chr2", 510, 520, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8], sv2[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad2]]
		# tadsPerSV[sv2[8]] = [[tad1], [tad2]]
		# svData = [svGroups, tadsPerSV]
		
		
		# VCF example case with 3 translocations
		
		# tad1 = ["chr2", 320000, 322000, TAD("chr2", 320000, 322000)]
		# tad2 = ["chr13", 123400, 123500, TAD("chr13", 123400, 123500)]
		# tad3 = ["chr17", 198900, 199000, TAD("chr17", 198900, 199000)]
		# 
		# sv1 = ["chr2", 321681, 321681, "chr17", 198982, 198982, "sample1", "cancerTypeA", SV("chr2", 321681, 321681, "+", "chr17", 198982, 198982, "+", "sample1", "cancerTypeA", "intraChromosomal")]
		# sv2 = ["chr2", 321681, 321681, "chr13", 123456, 123456, "sample1", "cancerTypeA", SV("chr2", 321681, 321681, "-", "chr13", 123456, 123456, "+", "sample1", "cancerTypeA", "intraChromosomal")]
		# sv3 = ["chr13", 123456, 123456, "chr17", 198982, 198982, "sample1", "cancerTypeA", SV("chr13", 123456, 123456, "-", "chr17", 198982, 198982, "-", "sample1", "cancerTypeA", "intraChromosomal")]
		# 
		# svGroups = dict()
		# tadsPerSV = dict()
		# svGroups[sv1[8]] = [sv1[8], sv2[8], sv3[8]]
		# tadsPerSV[sv1[8]] = [[tad1], [tad3]]
		# tadsPerSV[sv2[8]] = [[tad1], [tad2]]
		# tadsPerSV[sv3[8]] = [[tad2], [tad3]]
		# svData = [svGroups, tadsPerSV]
		
		
		
		
		### TRANSLOCATIONS ###
		if svType == "trans":
			
			svGroups = svData[0]
			tadsPerSV = svData[1]
			
			#1. For each group of SVs, determine a derivative genome
			#1.1 Keep a list of elements that are gained/lost for each translocation individually
			#1.2 For the next translocation, update that list again. If there are elements that are now lost, remove them from the list
			
			updatedTadPos = dict() #keep the TADs and the new start/ends after 
			
			for group in svGroups:
				
				
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
				fullNewTads = [] #Keep the TADs as a result of one translocation separate. This is a test for now, assuming that there are not more than 1 translocations re-affecting these TADs
				previousPieces = [] #Keep all previously remaining parts, to check if a previous SV already put these back together or not. 
				for sv in svGroups[group]:
					#print "SV: ", sv.chr1, sv.s1, sv.e1, sv.chr2, sv.s2, sv.e2, sv.sampleName, sv.o1, sv.o2
					
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
							
							if tad != leftTad and tad != rightTad:
							
								filteredInvolvedTads.append([tad.chromosome, tad.start, tad.end, tad])
						
					else: #get the affected TADs from the updated TADs.
						#print "SV pos is: ", sv.s1, sv.e2
						for tad in updatedTads: 
							
							if sv.s1 <= tad[2] and sv.s1 >= tad[1]:
								leftTad = tad[3]
							if sv.e2 <= tad[2] and sv.e2 >= tad[1]:
								rightTad = tad[3]
						#This is the extra check to make sure that we also include TADs that are not affected in the first iteration
						if leftTad is None:
							for tad in filteredInvolvedTads:
								if sv.s1 < tad[2] and sv.s1 > tad[1]:
									leftTad = tad[3]
						if rightTad is None:
							for tad in filteredInvolvedTads:
								if sv.e2 < tad[2] and sv.e2 > tad[1]:
									rightTad = tad[3]
						
						
						#Only check in the full TADs if the left/right TADs are still None. This is to prioritize looking for smaller pieces that are joined by SVs rather than re-updating an already existing TAD. 
						if leftTad is None:
							for tad in fullNewTads:
								if sv.s1 < tad[2] and sv.s1 > tad[1]:
									leftTad = tad[3]
						if rightTad is None:
							for tad in fullNewTads:
								if sv.s1 < tad[2] and sv.s1 > tad[1]:
									rightTad = tad[3]
						
						#If then the TADs are still None, continue. It may not be in a TAD
						if rightTad is None or leftTad is None:
							continue
					
					#These are the scenarios that we need to incorporate for translocations:
					#For +-, genes in the left side of chr1 gain eQTLs from the right side of chr2, and vice versa. 
					#For -+, genes in the left side of chr2 gain eQTLs from the right side of chr1, and vice versa. 
					#For ++, genes in the left side of chr1 gain eQTLs from the inverted left side of chr2, so starting from the breakpoint until the TAD start. Vice versa for the genes in the other part. 
					#For --, genes in the right side of chr1 gain eQTLs from the inverted right side of chr2, so starting from the breakpoint until the TAD end. Vice versa for the genes in the other part.

					if sv.o1 == "+" and sv.o2 == "-": #These are the same scenarios, but with the cases above we have already determined the correct TADs. 

						#2. Get the elements on the left of the breakpoint and the right of the breakpoint
						leftSideElements = leftTad.getElementsByRange(leftTad.start, sv.s1)
						leftSideGenes = leftTad.getGenesByRange(leftTad.start, sv.s1)
						
						rightSideElements = rightTad.getElementsByRange(sv.e2, rightTad.end)
						rightSideGenes = rightTad.getGenesByRange(sv.e2, rightTad.end)
						
						# #Also get the elements in the remaining part
						remainingElementsLeft = leftTad.getElementsByRange(sv.s1, leftTad.end)
						remainingGenesLeft = leftTad.getGenesByRange(sv.s1, leftTad.end)
						
						remainingElementsRight = rightTad.getElementsByRange(rightTad.start, sv.e2)
						remainingGenesRight = rightTad.getGenesByRange(rightTad.start, sv.e2)
						
						#3. Make derivative TADs from the SV and add elements, keeping the reference positions.
					
						#The new TADs will be the left part together with the right part, and everything between that remains a separate part (to be determined by another SV)
						newTad = [leftTad.chromosome, leftTad.start, rightTad.end, TAD(leftTad.chromosome, leftTad.start, rightTad.end)]
						
						
						#add the elements to the TAD
						newTad[3].setElements(leftSideElements + rightSideElements)
						newTad[3].setGenes(leftSideGenes + rightSideGenes)
						
						#Also set the elements to the remaining parts
						
						
						
						
						#The right part of the first chromosome
						remainingPart1 = [leftTad.chromosome, sv.s1, leftTad.end, TAD(leftTad.chromosome, sv.s1, leftTad.end)]
						remainingPart1[3].setElements(remainingElementsLeft)
						remainingPart1[3].setGenes(remainingGenesLeft)
						
						if [leftTad.chromosome, leftTad.start, leftTad.end, leftTad] in updatedTads:
							
							#Remove the remaining part, it was now added to another part of DNA.
							
							ind = updatedTads.index([leftTad.chromosome, leftTad.start, leftTad.end, leftTad])
							del updatedTads[ind]
						#But only re-add the remaining part if it is not exactly the same as a previous TAD (untouched)
						else:
						#if remainingPart1 not in updatedTads and remainingPart1 not in previousPieces:
							updatedTads.append(remainingPart1)
							
						
						
						remainingPart2 = [rightTad.chromosome, rightTad.start, sv.e2, TAD(rightTad.chromosome, rightTad.start, sv.e2)]
						remainingPart2[3].setElements(remainingElementsRight)
						remainingPart2[3].setGenes(remainingGenesRight)
						
						#Also remove the old TADs if necessary
						if [rightTad.chromosome, rightTad.start, rightTad.end, rightTad] in updatedTads:
							
							#Remove the remaining part, it was now added to another part of DNA.
							
							ind = updatedTads.index([rightTad.chromosome, rightTad.start, rightTad.end, rightTad])
							del updatedTads[ind]
						else:
						#if remainingPart2 not in updatedTads and remainingPart2 not in previousPieces:
							updatedTads.append(remainingPart2)

						#Also make sure that we remove the previously made TADs if these are updated again. 
						if [leftTad.chromosome, leftTad.start, leftTad.end, leftTad] in fullNewTads:
							ind = fullNewTads.index([leftTad.chromosome, leftTad.start, leftTad.end, leftTad])
							del fullNewTads[ind]
						if [rightTad.chromosome, rightTad.start, rightTad.end, rightTad] in fullNewTads:
							ind = fullNewTads.index([rightTad.chromosome, rightTad.start, rightTad.end, rightTad])
							del fullNewTads[ind]
							
						#4. Store the new TADs.
					
						#Copy all TADs but the left and right affected one to a new set of TADs.
						fullNewTads.append(newTad)
					
					if sv.o1 == "-" and sv.o2 == "+":
						#2. Get the elements on the left of the breakpoint and the right of the breakpoint
						#Left side is the first chromosome, right side the second chromosome. 
						leftSideElements = leftTad.getElementsByRange(sv.s1, leftTad.end)
						leftSideGenes = leftTad.getGenesByRange(sv.s1, leftTad.end)
						
						rightSideElements = rightTad.getElementsByRange(rightTad.start, sv.e2)
						rightSideGenes = rightTad.getGenesByRange(rightTad.start, sv.e2)

						# #Also get the elements in the remaining part
						remainingElementsLeft = leftTad.getElementsByRange(leftTad.start, sv.s1)
						remainingGenesLeft = leftTad.getGenesByRange(leftTad.start, sv.s1)
						
						remainingElementsRight = rightTad.getElementsByRange(sv.e2, rightTad.end)
						remainingGenesRight = rightTad.getGenesByRange(sv.e2, rightTad.end)

						#3. Make derivative TADs from the SV and add elements, keeping the reference positions.
						
						#The new TADs will be the left part together with the right part, and everything between that remains a separate part (to be determined by another SV)
						newTad = [rightTad.chromosome, rightTad.start, leftTad.end, TAD(rightTad.chromosome, rightTad.start, leftTad.end)]
						
						
						#add the elements to the TAD
						newTad[3].setElements(leftSideElements + rightSideElements)
						newTad[3].setGenes(leftSideGenes + rightSideGenes)
						
						#Also set the elements to the remaining parts
						
						#There are always 2 remaining parts, on both ends (in case of interchromosomal translocations)
						
						#But since we were already using smaller pieces, there are no remaining pieces, since the other one was already added.
						#So, check if the TADs that we are attaching are part of the previous remaining pieces. (this may go wrong, not check for all translocation sccenarios)
						
						#Only add the remaining parts if the current TAD is not already a previous remaining part. In that case, it would be in the TADs. 
						
						#The left part of the first chromosome
						remainingPart1 = [leftTad.chromosome, leftTad.start, sv.s1, TAD(leftTad.chromosome, leftTad.start, sv.s1)]
						remainingPart1[3].setElements(remainingElementsLeft)
						remainingPart1[3].setGenes(remainingGenesLeft)
						
						if [leftTad.chromosome, leftTad.start, leftTad.end, leftTad] in updatedTads:
							
							#Remove the remaining part, it was now added to another part of DNA.
							
							ind = updatedTads.index([leftTad.chromosome, leftTad.start, leftTad.end, leftTad])
							del updatedTads[ind]
						else:
						#if remainingPart1 not in updatedTads and remainingPart1 not in previousPieces:
							updatedTads.append(remainingPart1)
							
						remainingPart2 = [rightTad.chromosome, sv.e2, rightTad.end, TAD(rightTad.chromosome, sv.e2, rightTad.end)]
						remainingPart2[3].setElements(remainingElementsRight)
						remainingPart2[3].setGenes(remainingGenesRight)
						
						#Also remove the old TADs if necessary
						if [rightTad.chromosome, rightTad.start, rightTad.end, rightTad] in updatedTads:
							
							#Remove the remaining part, it was now added to another part of DNA.
							
							ind = updatedTads.index([rightTad.chromosome, rightTad.start, rightTad.end, rightTad])
							del updatedTads[ind]
						else:
						#if remainingPart2 not in updatedTads and remainingPart2 not in previousPieces:
							updatedTads.append(remainingPart2)

						#Also make sure that we remove the previously made TADs if these are updated again. 
						if [leftTad.chromosome, leftTad.start, leftTad.end, leftTad] in fullNewTads:
							ind = fullNewTads.index([leftTad.chromosome, leftTad.start, leftTad.end, leftTad])
							del fullNewTads[ind]
						if [rightTad.chromosome, rightTad.start, rightTad.end, rightTad] in fullNewTads:
							ind = fullNewTads.index([rightTad.chromosome, rightTad.start, rightTad.end, rightTad])
							del fullNewTads[ind]

						#4. Store the new TADs.
					
						#Copy all TADs but the left and right affected one to a new set of TADs.
						fullNewTads.append(newTad)
						
					if sv.o1 == "+" and sv.o2 == "+": 
						leftSideElements = leftTad.getElementsByRange(leftTad.start, sv.s1)
						leftSideGenes = leftTad.getGenesByRange(leftTad.start, sv.s1)
						
						#This part is inverted, so we start from the SV until the TAD start
						rightSideElements = rightTad.getElementsByRange(rightTad.start, sv.e2)
						rightSideGenes = rightTad.getGenesByRange(rightTad.start, sv.e2)
						
						# #Also get the elements in the remaining part
						#Left side is the first chromosome, right side the second chromosome.
						remainingElementsLeft = leftTad.getElementsByRange(sv.s1, leftTad.end)
						remainingGenesLeft = leftTad.getGenesByRange(sv.s1, leftTad.end)
						
						remainingElementsRight = rightTad.getElementsByRange(rightTad.start, sv.e2)
						remainingGenesRight = rightTad.getGenesByRange(rightTad.start, sv.e2)

						#3. Make derivative TADs from the SV and add elements, keeping the reference positions.
					
						#The new TAD is here the left part until the SV, and the right TAD from the start until the SV.
						
						newTad1 = [leftTad.chromosome, leftTad.start, rightTad.start, TAD(leftTad.chromosome, leftTad.start, rightTad.start)]
						
						#add the elements to the TAD
						newTad1[3].setElements(leftSideElements + rightSideElements)
						newTad1[3].setGenes(leftSideGenes + rightSideGenes)
						
						#Also set the elements to the remaining parts
						
						#There are always 2 remaining parts, on both ends (in case of interchromosomal translocations)
						
						remainingPart1 = [leftTad.chromosome, sv.s1, leftTad.end, TAD(leftTad.chromosome, sv.s1, leftTad.end)]
						remainingPart1[3].setElements(remainingElementsLeft)
						remainingPart1[3].setGenes(remainingGenesLeft)
						
						if [leftTad.chromosome, leftTad.start, leftTad.end, leftTad] in updatedTads:
							
							#Remove the remaining part, it was now added to another part of DNA.
							
							ind = updatedTads.index([leftTad.chromosome, leftTad.start, leftTad.end, leftTad])
							del updatedTads[ind]
						else:
						#if remainingPart1 not in updatedTads and remainingPart1 not in previousPieces:
							updatedTads.append(remainingPart1)
							
						remainingPart2 = [rightTad.chromosome, sv.e2, rightTad.end, TAD(rightTad.chromosome, sv.e2, rightTad.end)]
						remainingPart2[3].setElements(remainingElementsRight)
						remainingPart2[3].setGenes(remainingGenesRight)
						
						#Also remove the old TADs if necessary
						if [rightTad.chromosome, rightTad.start, rightTad.end, rightTad] in updatedTads:
							
							#Remove the remaining part, it was now added to another part of DNA.
							
							ind = updatedTads.index([rightTad.chromosome, rightTad.start, rightTad.end, rightTad])
							del updatedTads[ind]
						else:
						#if remainingPart2 not in updatedTads and remainingPart2 not in previousPieces:
							updatedTads.append(remainingPart2)
						#Also make sure that we remove the previously made TADs if these are updated again. 
						if [leftTad.chromosome, leftTad.start, leftTad.end, leftTad] in fullNewTads:
							ind = fullNewTads.index([leftTad.chromosome, leftTad.start, leftTad.end, leftTad])
							del fullNewTads[ind]
						if [rightTad.chromosome, rightTad.start, rightTad.end, rightTad] in fullNewTads:
							ind = fullNewTads.index([rightTad.chromosome, rightTad.start, rightTad.end, rightTad])
							del fullNewTads[ind]
						
						
						#4. Store the new TADs.
					
						#Copy all TADs but the left and right affected one to a new set of TADs.
						fullNewTads.append(newTad1)
						
						
					if sv.o1 == "-" and sv.o2 == "-":
						#This is the left part of chr1
						leftSideElements = leftTad.getElementsByRange(sv.s1, leftTad.end)
						leftSideGenes = leftTad.getGenesByRange(sv.s1, leftTad.end)
						
						#This part is inverted, so we start SV until the TAD end
						rightSideElements = rightTad.getElementsByRange(sv.e2, rightTad.end)
						rightSideGenes = rightTad.getGenesByRange(sv.e2, rightTad.end)
						
						
						# #Also get the elements in the remaining part
						#Left side is the first chromosome, right side the second chromosome.
						remainingElementsLeft = leftTad.getElementsByRange(leftTad.start, sv.s1)
						remainingGenesLeft = leftTad.getGenesByRange(leftTad.start, sv.s1)
						
						remainingElementsRight = rightTad.getElementsByRange(rightTad.start, sv.e2)
						remainingGenesRight = rightTad.getGenesByRange(rightTad.start, sv.e2)
						
						
						#3. Make derivative TADs from the SV and add elements, keeping the reference positions.
					
						#The new TAD is here the left part until the SV, and the right TAD from the start until the SV.
						
						newTad1 = [leftTad.chromosome, rightTad.end, leftTad.end, TAD(leftTad.chromosome, rightTad.end, leftTad.end)]
						
						#add the elements to the TAD
						newTad1[3].setElements(leftSideElements + rightSideElements)
						newTad1[3].setGenes(leftSideGenes + rightSideGenes)
						
						#There are always 2 remaining parts, on both ends (in case of interchromosomal translocations)
						
						remainingPart1 = [leftTad.chromosome, leftTad.start, sv.s1, TAD(leftTad.chromosome, leftTad.start, sv.s1)]
						remainingPart1[3].setElements(remainingElementsLeft)
						remainingPart1[3].setGenes(remainingGenesLeft)
						
						if [leftTad.chromosome, leftTad.start, leftTad.end, leftTad] in updatedTads:
							
							#Remove the remaining part, it was now added to another part of DNA.
							
							ind = updatedTads.index([leftTad.chromosome, leftTad.start, leftTad.end, leftTad])
							del updatedTads[ind]
						else:
							updatedTads.append(remainingPart1)

							
						remainingPart2 = [rightTad.chromosome, rightTad.start, sv.e2, TAD(rightTad.chromosome, rightTad.start, sv.e2)]
						remainingPart2[3].setElements(remainingElementsRight)
						remainingPart2[3].setGenes(remainingGenesRight)
						
						#Also remove the old TADs if necessary
						
						if [rightTad.chromosome, rightTad.start, rightTad.end, rightTad] in updatedTads:
							
							#Remove the remaining part, it was now added to another part of DNA.
						
							ind = updatedTads.index([rightTad.chromosome, rightTad.start, rightTad.end, rightTad])
							del updatedTads[ind]
						else:
							updatedTads.append(remainingPart2)

						
						#Also make sure that we remove the previously made TADs if these are updated again. 
						if [leftTad.chromosome, leftTad.start, leftTad.end, leftTad] in fullNewTads:
							ind = fullNewTads.index([leftTad.chromosome, leftTad.start, leftTad.end, leftTad])
							del fullNewTads[ind]
						if [rightTad.chromosome, rightTad.start, rightTad.end, rightTad] in fullNewTads:
							ind = fullNewTads.index([rightTad.chromosome, rightTad.start, rightTad.end, rightTad])
							del fullNewTads[ind]
						
							
						#4. Store the new TADs.
					
						#Copy all TADs but the left and right affected one to a new set of TADs.
						fullNewTads.append(newTad1)
					
				#5. For the next SV, use the updated TADs to get the elements/positions and continue until the last SV
					
				#6. Eventually, check for the genes in the original TADs which elements thesd has, and compare it to the derivative. Use the delta informatin to collect gains and losses
				
				#For each of the TADs, see which genes are in them and are affected, and determine their eQTL interactions.
				
				######	
				
					allTads = updatedTads + fullNewTads
					for tad in allTads:
						
						if tad[1] == tad[2]: #skip TADs that have no length
							continue
						
						for gene in tad[3].genes:
							
							#Get the elements that are inside this TAD for this gene	
							elementsInTad = []
							for element in gene.elements:
								elementStr = element[0] + "_" + str(element[1]) + "_" + str(element[2]) + "_" + str(element[3]) + "_" + str(element[4])
								
								#elementsInTad.append([element.chromosome, element.start, element])
								elementsInTad.append(elementStr)
							
							originalTadElements = []
							for element in gene.leftTAD.elements: #assume left TAD is the same as the right TAD, because the gene is within a TAD
								elementStr = element[0] + "_" + str(element[1]) + "_" + str(element[2]) + "_" + str(element[3]) + "_" + str(element[4])
								originalTadElements.append(elementStr)
							
							#Get all elements that are now in this TAD
							newElements = []
							for element in tad[3].elements:
								#if eQTL not in gene.leftTAD.eQTLInteractions:
								elementStr = element[0] + "_" + str(element[1]) + "_" + str(element[2]) + "_" + str(element[3]) + "_" + str(element[4])
								
								newElements.append(elementStr)
							
							
							filteredNewElements = np.setdiff1d(newElements, originalTadElements)
							#Make the derivative, determine which eQTLs are gained and lost for this gene.
							lostElements = []
							for element in elementsInTad:
								if element not in newElements: #this eQTL has been lost. 
									lostElements.append(element.split("_"))
							
							gene.addLostElements(lostElements, sv.sampleName)
							gene.addLostElementsSVs(lostElements, sv.chr1 + "_" + str(sv.s1) + "_" + str(sv.e1) + "_" + sv.chr2 + "_" + str(sv.s2) + "_" + str(sv.e2) + "_" + sv.sampleName)
							
							gainedElements = []
							for element in filteredNewElements:
								if element not in elementsInTad:
									#if eQTL[2] not in gene.leftTAD.eQTLInteractions: #Also do this check to make sure that we do not add eQTLs that were in the oroginal TAD of the gene, which are simply not associated to the gene. 
									gainedElements.append(element.split("_"))
							
							if len(gainedElements) < 1:
								continue
							
							gene.addGainedElements(gainedElements, sv.sampleName)
		
		### DELETIONS ###
		if svType == "del":
			
			#1. For every SV, determine the TADs that it occurs in.
			#2. Collect all elements within the SV region of these TADs
			#3. Assign the elements as lost for these genes (ignore if the gene itself is also deleted for now, that could be a different feature). If the elements were linked to genes and are lost, or if these are TAD-wide, is determined
			#by the gene object itself. 
			
			#Determine all overlapping TADs.
			
			tadChrSubsetInd = svData[0] == tadData[:,0]
			tadChrSubset = tadData[tadChrSubsetInd]
			
			#If the SV start is before the end of the TAD, and the SV end after the start of the TAD, the TAD is overlapped.
			startMatches = svData[1] <= tadChrSubset[:,2]
			endMatches = svData[5] >= tadChrSubset[:,1]
			
			tadMatches = tadChrSubset[startMatches * endMatches]
			
			if tadMatches.shape[0] < 1: #no matches
				return
								
			#Filter for TADs that are entirely overlapped and TADs that only contain part of the deletion.
			
			for tad in tadMatches:
				
				lostElements = []
				remainingGenes = []
				if svData[1] > tad[1] or svData[5] < tad[2]: #if the SV overlaps the TAD entirely, this will never be true.
					#Determine which part of the TAD is disrupted by the SV
					if svData[1] > tad[1] and svData[5] > tad[2]: #If the SV starts after the TAD start, but the TAD ends before the SV end, the SV is in the leftmost TAD.
						lostElements = tad[3].getElementsByRange(svData[1], tad[2]) #Elements in the deletion
						remainingGenes = tad[3].getGenesByRange(tad[1], svData[1])
						#For now, all genes will lose these elements, not just the genes that remain in the TAD. At a later point we could perhaps exclude the genes that are also deleted.  
					if svData[5] > tad[1] and svData[5] < tad[2]: #If the SV ends after the start of the TAD, and also ends before the end of the TAD, the SV is in the rightmost TAD.
						lostElements = tad[3].getElementsByRange(tad[1], svData[5])
						remainingGenes = tad[3].getGenesByRange(svData[5], tad[2])
				# else: #do not do this, if the deletion covers the entire TAD, it also removes the gene itself, so that is not intersting as a 'loss' 
				# 	lostElements = tad[3].elements
				for gene in remainingGenes:
					gene.addLostElements(lostElements, svData[8].sampleName)
					gene.addLostElementsSVs(lostElements, svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName)
					
			
					
								
							
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
				
			#Assigning the gains and losses to the genes is independent of the type of inversion
			#print "Copying genes and elements after SV"		
			#All genes that were originally in the left TAD (outisde of the inversion) will gain elements of the right side of the inversion
			#All unaffected genes on the left will lose the eQTLs that are in the left side of the inversion
			for gene in unaffectedGenesLeft:
				
				gene.addGainedElements(rightSideElements, svData[7])
				
				gene.addLostElements(leftSideElements, svData[7])
				gene.addLostElementsSVs(lostElements, svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName)
				
			#All genes in the right side of the inversion will gain elements from the original left TAD.
			#All genes in the right side will lose interactions with eQTLs in the unaffected right TAD. 
			for gene in rightSideGenes:
				
				gene.addGainedElements(unaffectedElementsLeft, svData[7])
				#print "Number of unaffected elements right: ", len(unaffectedElementsRight), " for genes ", len(rightSideGenes)
				gene.addLostElements(unaffectedElementsRight, svData[7])
				gene.addLostElementsSVs(lostElements, svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName)
			
			#vice versa but then for the right TAD and right side of the inversion.
			#The lost eQTLs are the ones that are in the right side of the inversion
			for gene in unaffectedGenesRight:
				
				gene.addGainedElements(leftSideElements, svData[7])
				gene.addLostElements(rightSideElements, svData[7])
				gene.addLostElementsSVs(lostElements, svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName)
			
			#The lost eQTLs are the ones that are in the unaffected original left TAD
			for gene in leftSideGenes:
				
				gene.addGainedElements(unaffectedElementsRight, svData[7])
				gene.addLostElements(unaffectedElementsLeft, svData[7])
				gene.addLostElementsSVs(lostElements, svData[0] + "_" + str(svData[1]) + "_" + str(svData[2]) + "_" + svData[3] + "_" + str(svData[4]) + "_" + str(svData[5]) + "_" + svData[8].sampleName)
			
			
				
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
			#- We skip cases where the duplication does not end in a TAD on either side. 
			
			
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
				
					#Assign the elements to the new TADs in the right order.
					#The first TAD gets the eQTLs within the SV of the last TAD.
					#The last TAD gets the eQTLs within the SV of the last TAD.
					
					svInteractionsFirstTad = filteredTads[0][3].getElementsByRange(svData[1], filteredTads[0][2])
					svInteractionsLastTad = filteredTads[len(filteredTads)-1][3].getElementsByRange(filteredTads[len(filteredTads)-1][2], svData[5])
					
					#Determine the gains for every gene. Also for the copied TADs, there are now multiple of these genes. 
					
					#For the new TADs, this is the same principle as for the eQTLs.
					#For the duplicated TADs, we can do * 2 of the elements
					
					#For TAD 1, the first part of C can interact with the second half of A.
					svGenesFirstTad = filteredTads[0][3].getGenesByRange(svData[1], filteredTads[0][2])
					svGenesLastTad = filteredTads[len(filteredTads)-1][3].getGenesByRange(filteredTads[len(filteredTads)-1][2], svData[5])
					
					for gene in svGenesFirstTad:
						
						gene.addGainedElements(svInteractionsLastTad, svData[7])
					
					for gene in svGenesLastTad:
					
						gene.addGainedElements(svInteractionsFirstTad, svData[7])
					
					#The last TAD remains the same overall.
					#Only the TADs in the middle are duplicated.
					
					#Each gene gains the eQTLs that are within the TAD that the gene is located in. 
					
					for tad in followingTads:
						for gene in tad[3].genes:
							
							#1. Get all eQTLs within this TAD
							tadEQTLs = tad[3].elements
							
								
							#2. Filter these for the eQTLs of the gene
							gainedEQTLs = []
							for eQTL in tadEQTLs:
								if eQTL[4] == gene.name:
									gainedEQTLs.append(eQTL)
								# if eQTL in gene.elements:
								# 	gainedEQTLs.append(eQTL)
								# if gene in eQTL.genes:
								# 	gainedEQTLs.append(eQTL)
							#3. Add the eQTLs to the gene for the current sample
							gene.addGainedElements(gainedEQTLs, svData[7])
							
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
						gene.addGainedElements(svInteractionsSecondTad, svData[7])
						
					for gene in svGenesSecondTad:
						
						#Each gene here gains eQTLs from outside of the SV in the bin.
						gene.addGainedElements(svInteractionsFirstTad, svData[7])
		