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
		self.elements = []
		self.genes = []
		
	def setSVs(self, SVs): #All SVs that overlap with this TAD boundary (left or right depending on the gene. Is that safe?)
		
		self.SVs = SVs
	
	def setSNVs(self, SNVs):
		
		self.SNVs = SNVs
		
	def setInteractions(self, interactions): #All interactions that take place within this TAD
		self.interactions = interactions
		
	def setElements(self, elements): #All eQTL interactions that take place within the TAD. 
		self.elements = elements
	
	def addElements(self, elements):
		for element in elements:
			self.elements.append(list(element))
		
	def addGene(self, gene): #Function to add genes that are within the TAD
		self.genes.append(gene)
		
	def setGenes(self, genes):
		self.genes = genes
	
	def getElementsByRange(self, start, end):
		
		
		elementsInRange = []
		for element in self.elements:
			if element[1] >= start and element[2] <= end:
				elementsInRange.append(element)
			# 
			# if eQTL.start >= start and eQTL.start <= end:
			# 	elementsInRange.append(eQTL)
		
	
		return elementsInRange
	
	def getGenesByRange(self, start, end):
		genesInRange = []
		for gene in self.genes:
			if gene.start >= start and gene.start <= end or gene.end <= end and gene.end >= start:
				genesInRange.append(gene)
		
		return genesInRange
		
		