class Element:
	
	"""
		Generic class for genomic elements that can be described with a chromosome, start and end position. We assume that these elements can be linked to genes, but that is not necessary. 
		
	"""
	
	
	def __init__(self, chromosome, start, end):
		
		self.chromosome = chromosome
		self.start = start
		self.end = end
		self.enhancerStatus = "N" #y default not an enhancer
		self.SVs = []
		self.SNVs = []
		self.genes = [] #Make sure that we also know which gene(s) the eQTL affects
		self.type = "element"
		
	def setSVs(self, SVs):
		
		self.SVs = SVs
		
	def setSNVs(self, SNVs):
		
		self.SNVs = SNVs
		
	def addGene(self, gene):
		self.genes.append(gene)
