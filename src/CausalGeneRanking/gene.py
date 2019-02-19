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
		self.SVs = None
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

			#print "final no of gains: ", len(self.gainedEQTLs[sample])
		types = []
		for gainedElement in gainedElements:
			if gainedElement.type not in types:
				
				types.append(gainedElement.type)
		
		#For the recurrence ranking, speed up code by adding only 1 per sample. More is not necessary. 
		if len(gainedElements) > 0:
			if sample not in self.gainedElements:
				self.gainedElements[sample] = types
			else:
				for elementType in types:
					self.gainedElements[sample].append(elementType)
		
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

		#An eQTL can only be lost if it is associated with that specific gene. 	
		for lostElement in lostElements:
			if lostElement in self.elements:
				self.addLostElement(lostElement, sample)
		#self.lostEQTLs[sample] += lostEQTLs
		
	def setLostElements(self, lostElements, sample):
		for lostElement in lostElements:
			if lostElement in self.elements:
				self.addLostElement(lostElement, sample)
		#self.lostEQTLs[sample] = lostEQTLs
	
	def addLostElement(self, lostElement, sample):
	
		if settings.general['lncRNA'] == True:
			
			if sample not in self.lostElements:
				self.lostElements[sample] = []
			
			self.lostElements[sample].append(lostElement.type)
	
		else:	
			if lostElement in self.elements:
				if sample not in self.lostElements:
					self.lostElements[sample] = []
				
				self.lostElements[sample].append(lostElement.type)
		