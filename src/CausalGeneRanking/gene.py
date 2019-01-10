from copy import deepcopy

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
		self.eQTLs = []
		self.interactions = []
		self.gainedEQTLs = dict()
		self.lostEQTLs = dict()
		
	def setTADs(self, leftTAD, rightTAD):
		
		self.leftTAD = leftTAD
		self.rightTAD = rightTAD
	
	def setLeftTAD(self, leftTAD):
		self.leftTAD = leftTAD
		
	def setRightTAD(self, rightTAD):
		self.rightTAD = rightTAD
		
	
	def setEQTLS(self, eQTLS):
		
		self.eQTLS = eQTLS
		
	def setInteractions(self, interactions):
		self.interactions = interactions
		
	def setSVs(self, SVs):
		
		self.SVs = SVs
		
	def setSNVs(self, SNVs):
		self.SNVs = SNVs
		
	def addEQTL(self, eQTL):
		self.eQTLs.append(eQTL)
		
	def setGainedEQTLs(self, gainedEQTLs, sample):
		self.gainedEQTLs[sample] = gainedEQTLs #keep the gained eQTLs separate per patient to later do mutual exclusivity.
		
	def addGainedEQTLs(self, gainedEQTLs, sample):
		if sample not in self.gainedEQTLs:
			self.gainedEQTLs[sample] = deepcopy(gainedEQTLs)
		else:
			self.gainedEQTLs[sample] += deepcopy(gainedEQTLs)
		if len(gainedEQTLs) > 0:
			print "no of gains: ", len(gainedEQTLs)
			print "final no of gains: ", len(self.gainedEQTLs[sample])
		
	
	def addLostEQTLs(self, lostEQTLs, sample):
		if sample not in self.lostEQTLs:
			self.lostEQTLs[sample] = []
		self.lostEQTLs[sample] += lostEQTLs
		
		
		
	def setLostEQTLs(self, lostEQTLs, sample):
		self.lostEQTLs[sample] = lostEQTLs
	
	def addLostEQTL(self, lostEQTL, sample):
		if sample not in self.lostEQTLs:
			self.lostEQTLs[sample] = []	
		self.lostEQTLs[sample].append(lostEQTL)
	