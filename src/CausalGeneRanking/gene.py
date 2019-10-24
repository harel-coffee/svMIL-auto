from __future__ import absolute_import
from copy import deepcopy
import numpy as np

import settings

class Gene:
	"""
		Class to describe a gene. Holds all other information related to the neighborhood of the gene as well, like the TADs, eQTLs and SVs.
	"""
	def __init__(self, name, chromosome, start, end):
		
		self.name = name
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.SVs = dict()
		self.SNVs = None
		self.leftTAD = None
		self.rightTAD = None
		self.elements = []
		self.gainedElements = dict()
		self.lostElements = dict()
		self.lostElementsSVs = dict() #lost elements per SV, not per sample
		self.gainedElementsSVs = dict()
		self.alteredElements = dict()
		self.elementsNotLinkedToGenes = ['cpg', 'tf', 'hic', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
									'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin',
									'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed']
	def setTADs(self, leftTAD, rightTAD):
		
		self.leftTAD = leftTAD
		self.rightTAD = rightTAD
	
	#I'm not sure if the left and right TADs are really used anymore, needs to be checked
	def setLeftTAD(self, leftTAD):
		self.leftTAD = leftTAD
		
	def setRightTAD(self, rightTAD):
		self.rightTAD = rightTAD
		
	
	def setElements(self, elements):
		
		self.elements = elements
		
	def setSVs(self, SVs):
		
		self.SVs = SVs
		
	def setSNVs(self, SNVs):
		self.SNVs = SNVs
		
	def addElement(self, element):
		self.elements.append(element)
		
	def setGainedElements(self, gainedElements, sample):
		self.gainedElements[sample] = gainedElements #keep the gained eQTLs separate per patient to later do mutual exclusivity.
		
	def addGainedElements(self, gainedElements, sample):
		
		if len(gainedElements) > 0:
			if sample not in self.gainedElements:
				self.gainedElements[sample] = dict()
				
		#Have a dictionary where we count the number of elements of a specific type that are gained per sample.
		#This is much faster than storing the actual elements that are gained, and we do not use that information in the ranking, so it can be discarded here. 
		for gainedElement in gainedElements:
			if gainedElement[3] not in self.gainedElements[sample]:
				self.gainedElements[sample][gainedElement[3]] = 0
			self.gainedElements[sample][gainedElement[3]] += 1

		
	def addGainedElementsSVs(self, gainedElements, sv):
		
		if len(gainedElements) > 0:
			if sv not in self.gainedElements:
				self.gainedElementsSVs[sv] = dict()
		
		#Have a dictionary where we count the number of elements of a specific type that are gained per sample.
		#This is much faster than storing the actual elements that are gained, and we do not use that information in the ranking, so it can be discarded here. 
		for gainedElement in gainedElements:
			if gainedElement[3] not in self.gainedElementsSVs[sv]:
				self.gainedElementsSVs[sv][gainedElement[3]] = 0
			self.gainedElementsSVs[sv][gainedElement[3]] += 1
		
		self.addAlteredElements(gainedElements, sv, 'gain')	
	
	def addLostElements(self, lostElements, sample):
		
		if len(lostElements) > 0:
			if sample not in self.lostElements:
				self.lostElements[sample] = dict()
		
		#Have a dictionary where we count the number of elements of a specific type that are lost per sample.
		#This is much faster than storing the actual elements that are lost, and we do not use that information in the ranking, so it can be discarded here.
	

		for lostElement in lostElements:
			if lostElement[3] in self.elementsNotLinkedToGenes:
				if lostElement[3] not in self.lostElements[sample]:
					self.lostElements[sample][lostElement[3]] = 0
				self.lostElements[sample][lostElement[3]] +=1
			else:
				
				if lostElement[4] == self.name:#filter by elements that are linked to genes in the data, exclude these as losses if not linked to the gene
					if lostElement[3] not in self.lostElements[sample]: 
						self.lostElements[sample][lostElement[3]] = 0
					self.lostElements[sample][lostElement[3]] +=1
	
		
	
	def addLostElementsSVs(self, lostElements, sv):
		
		if len(lostElements) > 0:
			if sv not in self.lostElements:
				self.lostElementsSVs[sv] = dict()
		
		#Have a dictionary where we count the number of elements of a specific type that are lost per sample.
		#This is much faster than storing the actual elements that are lost, and we do not use that information in the ranking, so it can be discarded here.
		
		
		for lostElement in lostElements:
			if lostElement[3] in self.elementsNotLinkedToGenes:
				if lostElement[3] not in self.lostElementsSVs[sv]:
					self.lostElementsSVs[sv][lostElement[3]] = 0
				self.lostElementsSVs[sv][lostElement[3]] +=1
			else:
				
				if lostElement[4] == self.name:#filter by elements that are linked to genes in the data, exclude these as losses if not linked to the gene
					if lostElement[3] not in self.lostElementsSVs[sv]: 
						self.lostElementsSVs[sv][lostElement[3]] = 0
					self.lostElementsSVs[sv][lostElement[3]] +=1
		
		self.addAlteredElements(lostElements, sv, 'loss')
		
	def addAlteredElements(self, elements, sv, alterationType):
		"""
			For this gene, make a dictionary where we can look up by SV which elements were altered by that SV for this gene.
			The values for this element are the feature vector that we will use to describe that element. 
		"""
		
		allowedElements = ['enhancer']
		
		elementsNotLinkedToGenes = ['cpg', 'tf', 'hic', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
									'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin',
									'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol',
									'enhancer'] #add enhancer temporarily to work with unlinked data

		if len(elements) > 0:
			if sv not in self.alteredElements:
				self.alteredElements[sv] = dict()
			else:
				return #don't add the SV another time if we see it again for some reason, we already saw it once. 
		
		#For methylation marks, gather all relevant marks here for easy lookup.
		#For now, just focus on what is relevant for enhancers
		methylationMarks = []
		annotationElements = ['dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3', 'CTCF', 'CTCF+Enhancer', 'Enhancer',
							  'Heterochromatin', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol']
		for element in elements:
			
			#if element[3] in ['h3k27ac', 'h3k4me1']:
			#if element[3] in ['h3k27ac', 'h3k4me1', 'CTCF', 'CTCF+Enhancer', 'Enhancer', 'Heterochromatin', 'Repeat', 'Repressed', 'Transcribed', 'dnaseI', 'rnaPol']:
			if element[3] in annotationElements:
				methylationMarks.append(element)
		
		methylationMarks = np.array(methylationMarks, dtype='object')	
		
		for element in elements:
			
			if element[3] not in allowedElements:
				continue
			
			elementStr = element[0] + "_" + str(element[1]) + "_" + str(element[2]) + "_" + element[3]

			#Check if the element is methylated or not
			#first, just set the things for enhancers only, this should later be element un-specific
			#For all the other elements, determine if the enhancer has a specific mark or not (at the same position)
			
			#Find overlap with the methylated elements
			#Any overlap is accepted
			methylationMatches = []
			if len(methylationMarks) > 0:
				methylationMatches = methylationMarks[(methylationMarks[:,0] == element[0]) * (methylationMarks[:,2] >= element[1]) * (methylationMarks[:,1] <= element[2])]
			
			#Fix this later and make it not-so-hardcoded
			#elementMethylation = [0, 0] #keep order of elements
			elementMethylation = [0]*len(annotationElements) #keep order of elements
			for match in methylationMatches:
				
				for elementInd in range(0, len(annotationElements)):
					annotationElement = annotationElements[elementInd]
					if match[3] == annotationElement:
						elementMethylation[elementInd] = 1
				
			lossGains = [0,0]
			if alterationType == 'loss':
				if element[3] in self.elementsNotLinkedToGenes:
					lossGains[0] = 1
				else: #make sure that elements that belong to the gene are only lost. 
					if element[4] == self.name:
						lossGains[0] = 1
					else: #if the loss is from an element that was not interacting with this gene, it is not a true loss. 
						lossGains[0] = 0

			if alterationType == "gain":
				lossGains[1] = 1
			
			#if we get here, we passed all checks and there is a valid gain OR loss
			if elementStr not in self.alteredElements[sv]:
				self.alteredElements[sv][elementStr] = lossGains + elementMethylation
				#self.alteredElements[sv][elementStr] = lossGains
		#something with methylation for the affected genes only
		#first make sure that all elements are gathered, then afterwards, add the methylation specifically for each of them. 
		#methylationData = InputParser().getMethylationFromFile('../../data/methylation/BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt')
		
		
		
		