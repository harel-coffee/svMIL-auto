from __future__ import absolute_import
from copy import deepcopy
import numpy as np

import settings

class Gene:
	"""
		Class to describe a gene. Holds all other information related to the neighborhood of the gene as well, like the TADs, regulatory elements and SVs affecting the gene.
	"""
	def __init__(self, name, chromosome, start, end):
		
		self.name = name
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.SVs = dict()
		self.SNVs = []
		self.CNVs = []
		self.leftTAD = None
		self.rightTAD = None
		self.elements = []
		self.gainedElements = dict() #gained elements per sample
		self.lostElements = dict()
		self.lostElementsSVs = dict() #lost elements per SV, not per sample
		self.gainedElementsSVs = dict()
		self.gainedElementsStrengthsSVs = dict()
		self.lostElementsStrengthsSVs = dict()
		self.alteredElements = dict()
		self.elementsNotLinkedToGenes = ['ctcf', 'rnaPol', 'cpg', 'tf', 'hic', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
									'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin',
									'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'superEnhancer']
		self.strengthElements = ['enhancer', 'ctcf', 'rnaPol', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
		self.cosmic = 0
		
		
	def setTADs(self, leftTAD, rightTAD):
		
		self.leftTAD = leftTAD
		self.rightTAD = rightTAD

	#Left/right TADs are not used anymore, but this would be the TADs on the left/right of the gene, especially for those that are in multiple TADs.
	def setLeftTAD(self, leftTAD):
		self.leftTAD = leftTAD
		
	def setRightTAD(self, rightTAD):
		self.rightTAD = rightTAD
		
	
	def setElements(self, elements):
		
		self.elements = elements
		
	def setSVs(self, SVs):
		
		self.SVs = SVs
		
	def addElement(self, element):
		self.elements.append(element)
		
	def setGainedElements(self, gainedElements, sample): #set the gained elements per sample. 
		self.gainedElements[sample] = gainedElements
		
	def addGainedElements(self, gainedElements, sample): #add elements one by one. 
		
		if len(gainedElements) > 0:
			if sample not in self.gainedElements:
				self.gainedElements[sample] = dict()
				
		#Have a dictionary where we count the number of elements of a specific type that are gained per sample.
		#This is much faster than storing the actual elements that are gained, and we do not use that information in the ranking, so it can be discarded here. 
		for gainedElement in gainedElements:
			if gainedElement[3] not in self.gainedElements[sample]:
				self.gainedElements[sample][gainedElement[3]] = 0
			self.gainedElements[sample][gainedElement[3]] += 1


	#this is what we use in the gene-SV pair output file. The provided SV (object) causes this gene to gain these elements.
	def addGainedElementsSVs(self, gainedElements, sv):
		
		if len(gainedElements) > 0:
			if sv not in self.gainedElementsSVs:
				self.gainedElementsSVs[sv] = dict()
		
		#Have a dictionary where we count the number of elements of a specific type that are gained per sample.
		#This is much faster than storing the actual elements that are gained, and we do not use that information in the ranking, so it can be discarded here. 
		for gainedElement in gainedElements:
			
			if gainedElement[3] not in self.gainedElementsSVs[sv]:
				self.gainedElementsSVs[sv][gainedElement[3]] = 0
			self.gainedElementsSVs[sv][gainedElement[3]] += 1
		
		
		self.addAlteredElements(gainedElements, sv, 'gain') #add to bags
		self.addGainedElementsStrengthsSVs(gainedElements, sv) #get the strengths of the elements that are gained

	#Add the strengths of the gained elements. Because there could be multiple, take the max
	#This was only intended for the SV-gene pairs .txt output file, but the strentghs are there not used anymore.
	#For MIL the strengths are implemented differently! (see addAlteredElements)
	def addGainedElementsStrengthsSVs(self, gainedElements, sv):
		
		#for each element type, get the strongest one, and keep it as a feature.
		if len(gainedElements) > 0:
			if sv not in self.gainedElementsStrengthsSVs:
				self.gainedElementsStrengthsSVs[sv] = dict()

		#Find the element with the maximum strength for this element type. 
		strengthList = dict()
		for gainedElement in gainedElements:
			
			if gainedElement[3] in self.strengthElements:
				
				if gainedElement[3] not in self.gainedElementsStrengthsSVs[sv]:
					self.gainedElementsStrengthsSVs[sv][gainedElement[3]] = 0
				
				if gainedElement[5] > self.gainedElementsStrengthsSVs[sv][gainedElement[3]]:
					self.gainedElementsStrengthsSVs[sv][gainedElement[3]] = gainedElement[5]
				
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
			if sv not in self.lostElementsSVs:
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
		
		self.addAlteredElements(lostElements, sv, 'loss') #for MIL, add lost elements to the bags. 
		self.addLostElementsStrengthsSVs(lostElements, sv)
	
	#Add the strengths of the lost elements. Because there could be multiple, take the max
	#This was only intended for the SV-gene pairs .txt output file, but the strentghs are there not used anymore.
	#For MIL the strengths are implemented differently! (see addAlteredElements)
	def addLostElementsStrengthsSVs(self, lostElements, sv):
		
		#for each element type, get the strongest one, and keep it as a feature.
		if len(lostElements) > 0:
			if sv not in self.lostElementsStrengthsSVs:
				self.lostElementsStrengthsSVs[sv] = dict()

		#list for each type what the max value is so far
		strengthList = dict()
		for lostElement in lostElements:

			if lostElement[3] in self.strengthElements:
				
				
				if lostElement[3] not in self.lostElementsStrengthsSVs[sv]:
					self.lostElementsStrengthsSVs[sv][lostElement[3]] = 0
				
				if lostElement[5] > self.lostElementsStrengthsSVs[sv][lostElement[3]]:
					self.lostElementsStrengthsSVs[sv][lostElement[3]] = lostElement[5]
					
			
	def addAlteredElements(self, elements, sv, alterationType):
		"""
			This is where we get what we need to set up MIL.
			This gene is affected by the provided SV, and gains or loses the provided element (alterationType specifies which).
			
			To gather this information, for this gene, the altered elements are added to a dictionary. The keys are first the SVs that disrupt the gene.
			Then, within that, there is a dictionary per element (str format), which then holds the features of the gains and losses (instances).

			Later on, we go through this gene, and make bags for each SV affecting this gene (SV-gene pair as bag label).

			elements (list): list of all regulatory element disrupted by the SV for this gene
			sv (str): str format of the SV that affects this gene (see DerivativeTADMaker for format)
			alterationType (str): either 'gain' or 'loss'.
			
		"""
		
		allowedElements = ['enhancer', 'eQTL', 'superEnhancer']
		#allowedElements = ['enhancer', 'eQTL']

		if len(elements) > 0:
			if sv not in self.alteredElements:
				self.alteredElements[sv] = dict()

		#For methylation marks, gather all relevant marks here for easy lookup.
		methylationMarks = []
		annotationElements = ['cpg', 'tf', 'hic', 'ctcf', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3',
							  'CTCF', 'CTCF+Enhancer', 'CTCF+Promoter', 'Enhancer', 'Heterochromatin', 'Poised_Promoter', 'Promoter', 'Repeat', 'Repressed', 'Transcribed', 'rnaPol']

		strengthElements = ['enhancer', 'ctcf', 'rnaPol', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
		for element in elements:

			if element[3] in annotationElements:
				methylationMarks.append(element)

		methylationMarks = np.array(methylationMarks, dtype='object')

		for element in elements:

			if element[3] not in allowedElements:
				continue

			#skip these cases here already to improve run time. These element losses are never
			#accepted because the elements are not linked to the current gene.
			if alterationType == 'loss':
				if element[3] not in self.elementsNotLinkedToGenes and element[4] != self.name:
					continue

			elementStr = element[0] + "_" + str(element[1]) + "_" + str(element[2]) + "_" + element[3]

			#Check if the element is methylated or not

			#Find overlap with the methylated elements
			#Any overlap is accepted
			methylationMatches = []
			if len(methylationMarks) > 0:
				#methylationMatches = methylationMarks[(methylationMarks[:,0] == element[0]) * (methylationMarks[:,2] >= element[1]) * (methylationMarks[:,1] <= element[2])]

				methylationMatches = methylationMarks[(methylationMarks[:,2] >= element[1]) * (methylationMarks[:,1] <= element[2])]



			#Fix this later and make it not-so-hardcoded
			#elementMethylation = [0, 0] #keep order of elements
			elementMethylation = [0]*len(annotationElements) #keep order of elements
			elementStrength = [0]*len(strengthElements) #keep order of elements
			for match in methylationMatches:

				for elementInd in range(0, len(annotationElements)):
					annotationElement = annotationElements[elementInd]
					if match[3] == annotationElement:
						elementMethylation[elementInd] = 1

				for elementInd in range(0, len(strengthElements)):
					strengthElement = strengthElements[elementInd]
					if match[3] == strengthElement:
						elementStrength[elementInd] = match[5]



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

			if lossGains[0] == 0 and lossGains[1] == 0:
				continue #this pair ended up with a bad loss, so we need to skip it then if it also has no gains.

			enhancerType = 0
			eQTLType = 0
			superEnhancerType = 0

			if element[3] == 'enhancer':
				enhancerType = 1
			elif element[3] == 'eQTL':
				eQTLType = 1
			elif element[3] == 'superEnhancer':
				superEnhancerType = 1


			#if we get here, we passed all checks and there is a valid gain OR loss
			if elementStr not in self.alteredElements[sv]:
				self.alteredElements[sv][elementStr] = lossGains + elementMethylation + elementStrength + [enhancerType, eQTLType, superEnhancerType]
