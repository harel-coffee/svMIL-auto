from copy import deepcopy

class TAD:
	
	"""
		Class to describe the location of a TAD and potentially SVs overlapping it.
	"""
	
	def __init__(self, chromosome, start, end):
		
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.SVs = None
		self.SNVs = None
		self.interactions = None
		self.eQTLInteractions = []
		self.genes = []
		
	def setSVs(self, SVs): #All SVs that overlap with this TAD boundary (left or right depending on the gene. Is that safe?)
		
		self.SVs = SVs
	
	def setSNVs(self, SNVs):
		
		self.SNVs = SNVs
		
	def setInteractions(self, interactions): #All interactions that take place within this TAD
		self.interactions = interactions
		
	def setEQTLInteractions(self, eQTLInteractions): #All eQTL interactions that take place within the TAD. 
		self.eQTLInteractions = eQTLInteractions
		
	def addGene(self, gene): #Function to add genes that are within the TAD
		self.genes.append(gene)
		
	def setGenes(self, genes):
		self.genes = genes
	
	def getElementsByRange(self, start, end):
		
		#First do this only for eQTLs
		eQTLsInRange = []
		for eQTL in self.eQTLInteractions:
			if eQTL.start >= start and eQTL.start <= end:
				eQTLsInRange.append(eQTL)
		
		return eQTLsInRange
	
	def getGenesByRange(self, start, end):
		genesInRange = []
		for gene in self.genes:
			if gene.start >= start and gene.start <= end or gene.end <= end and gene.end >= start:
				genesInRange.append(gene)
		
		return genesInRange
		
		