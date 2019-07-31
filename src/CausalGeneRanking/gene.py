from copy import deepcopy

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
			
	
	def addLostElements(self, lostElements, sample):
		
		if len(lostElements) > 0:
			if sample not in self.lostElements:
				self.lostElements[sample] = dict()
		
		#Have a dictionary where we count the number of elements of a specific type that are lost per sample.
		#This is much faster than storing the actual elements that are lost, and we do not use that information in the ranking, so it can be discarded here.
		elementsNotLinkedToGenes = ['cpg', 'tf', 'hic', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
		
		for lostElement in lostElements:
			if lostElement[3] in elementsNotLinkedToGenes:
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
		elementsNotLinkedToGenes = ['cpg', 'tf', 'hic', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
		
		for lostElement in lostElements:
			if lostElement[3] in elementsNotLinkedToGenes:
				if lostElement[3] not in self.lostElementsSVs[sv]:
					self.lostElementsSVs[sv][lostElement[3]] = 0
				self.lostElementsSVs[sv][lostElement[3]] +=1
			else:
				
				if lostElement[4] == self.name:#filter by elements that are linked to genes in the data, exclude these as losses if not linked to the gene
					if lostElement[3] not in self.lostElementsSVs[sv]: 
						self.lostElementsSVs[sv][lostElement[3]] = 0
					self.lostElementsSVs[sv][lostElement[3]] +=1
		
	# def setLostElements(self, lostElements, sample):
	# 	
	# 	for lostElement in lostElements:
	# 		if lostElement in self.elements:
	# 			self.addLostElement(lostElement, sample)
	# 	#self.lostEQTLs[sample] = lostEQTLs
	# 
	# def addLostElement(self, lostElement, sample, types):
	# 	
	# 	#Treat losses differently for elements that we cannot link to the gene
	# 	elementsNotLinkedToGenes = ['cpg', 'tf', 'hic', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
	# 	
	# 	if sample not in self.lostElements:
	# 		self.lostElements[sample] = dict()
	# 
	# 	for elementType in types:
	# 		if elementType not in self.lostElements[sample]:
	# 			self.lostElements[sample][elementType] = 0
	# 		self.lostElements[sample][elementType] += 1
	# 
	# 	