from copy import deepcopy

import settings

class Gene:
	"""
		Class to describe a gene. Will hold all other information related to the neighborhood of the gene as well, like the TADs, eQTLs and SVs. 
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
		self.interactions = []
		self.gainedElements = dict()
		self.lostElements = dict()
		
	def setTADs(self, leftTAD, rightTAD):
		
		self.leftTAD = leftTAD
		self.rightTAD = rightTAD
	
	def setLeftTAD(self, leftTAD):
		self.leftTAD = leftTAD
		
	def setRightTAD(self, rightTAD):
		self.rightTAD = rightTAD
		
	
	def setElements(self, elements):
		
		self.elements = elements
		
	def setInteractions(self, interactions):
		self.interactions = interactions
		
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
		
		for gainedElement in gainedElements:
			if gainedElement[3] not in self.gainedElements[sample]:
				self.gainedElements[sample][gainedElement[3]] = 0
			self.gainedElements[sample][gainedElement[3]] += 1
			#print "final no of gains: ", len(self.gainedEQTLs[sample])
		# types = dict()
		# for gainedElement in gainedElements:
		# 	types[gainedElement[3]] = 0
		# 	#if gainedElement.type not in types:
		# 	
		# 	#	types.append(gainedElement.type)
		# 	# if gainedElement[3] not in types:
		# 	# 	types.append(gainedElement[3])
		# 
		# #For the recurrence ranking, speed up code by adding only 1 per sample. More is not necessary.
		# 
		# if len(gainedElements) > 0:
		# 	if sample not in self.gainedElements:
		# 		self.gainedElements[sample] = types.keys()
		# 	else:
		# 		for elementType in types:
		# 			
		# 			self.gainedElements[sample].append(elementType)
		# 	
		#Use this part for when we need all eQTLs in a list
		# if len(gainedEQTLs) > 0:
		# 	if sample not in self.gainedEQTLs:
		# 		self.gainedEQTLs[sample] = gainedEQTLs
		# 	else:
		# 		self.gainedEQTLs[sample] += gainedEQTLs
		# 
		# 
		# if sample in self.gainedEQTLs:
		# 	print "no of gains: ", len(gainedEQTLs)
		# 	print "final no of gains: ", len(self.gainedEQTLs[sample])
		# 	exit()
		
	
	def addLostElements(self, lostElements, sample):
		
		if sample not in self.lostElements:
			self.lostElements[sample] = dict()
		
		for lostElement in lostElements:
			if lostElement[3] not in self.lostElements[sample]:
				self.lostElements[sample][lostElement[3]] = 0
			self.lostElements[sample][lostElement[3]] +=1
			
			
			
			
		
		# for lostElement in lostElements:
		# 	self.addLostElement(lostElement, sample, types)
		# 
		
	def setLostElements(self, lostElements, sample):
		for lostElement in lostElements:
			if lostElement in self.elements:
				self.addLostElement(lostElement, sample)
		#self.lostEQTLs[sample] = lostEQTLs
	
	def addLostElement(self, lostElement, sample, types):

		#Treat losses differently for elements that we cannot link to the gene
		elementsNotLinkedToGenes = ['cpg', 'tf', 'hic', 'dnaseI', 'h3k9me3', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k36me3']
		
		# if lostElement[3] in elementsNotLinkedToGenes:
		# 	
		# 	if sample not in self.lostElements:
		# 		self.lostElements[sample] = types.keys()
		# 	else:
		# 		for elementType in types:
		# 			self.lostElements[sample].append(elementType)
		if sample not in self.lostElements:
			self.lostElements[sample] = dict()

		for elementType in types:
			if elementType not in self.lostElements[sample]:
				self.lostElements[sample][elementType] = 0
			self.lostElements[sample][elementType] += 1
		# else:	
		# 	if lostElement in self.elements:
		# 		if sample not in self.lostElements:
		# 			self.lostElements[sample] = types.keys()
		# 		else:
		# 			for elementType in types:
		# 				self.lostElements[sample].append(elementType)
		